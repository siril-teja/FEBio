//! febio_rs_solver: parallel sparse direct solver for FEBio written in Rust.
//!
//! * symmetric (indefinite allowed) matrices: supernodal LDL^T with METIS nested
//!   dissection ordering, tree + dense-kernel parallelism, PARDISO-style pivot
//!   perturbation and optional iterative refinement.
//! * unsymmetric (structurally symmetric) matrices: supernodal LU with static pivoting,
//!   pivot perturbation, residual-checked iterative refinement and a partial-pivoting
//!   sparse LU (faer) as automatic fallback.
//!
//! C ABI (see `febio_rs_solver.h`).

pub mod lu;
pub mod numeric;
pub mod symbolic;

use std::time::Instant;

pub use symbolic::Ordering;

pub struct Stats {
    pub n: usize,
    pub nnz_a: usize,
    pub nnz_l: usize,
    pub n_supernodes: usize,
    pub n_supervars: usize,
    pub flops: f64,
    pub t_analyze: f64,
    pub t_factor: f64,
    pub t_solve: f64,
    pub n_factor: usize,
    pub n_solve: usize,
    pub n_perturbed: usize,
    pub n_refine: usize,
    pub n_fallback: usize,
    pub n_analyze: usize,
    pub n_reused: usize,
}

impl Default for Stats {
    fn default() -> Self {
        Stats {
            n: 0,
            nnz_a: 0,
            nnz_l: 0,
            n_supernodes: 0,
            n_supervars: 0,
            flops: 0.0,
            t_analyze: 0.0,
            t_factor: 0.0,
            t_solve: 0.0,
            n_factor: 0,
            n_solve: 0,
            n_perturbed: 0,
            n_refine: 0,
            n_fallback: 0,
            n_analyze: 0,
            n_reused: 0,
        }
    }
}

impl Stats {
    /// new sparsity pattern: reset the pattern description, keep cumulative counters/timers
    fn reset_pattern(&mut self) {
        self.n_analyze += 1;
        self.n_supernodes = 0;
        self.n_supervars = 0;
        self.nnz_l = 0;
        self.flops = 0.0;
        self.n_perturbed = 0;
    }
}

enum Kind {
    None,
    Symmetric {
        colptr: Vec<usize>,
        rowidx: Vec<usize>,
        sym: symbolic::Symbolic,
        sched: (Vec<usize>, Vec<usize>),
        num: Option<numeric::Numeric>,
        vals: Vec<f64>,
    },
    Unsymmetric {
        n: usize,
        // CSR of A (sorted column indices) == CSC of A^T
        rowptr: Vec<usize>,
        colidx: Vec<usize>,
        /// map from FEBio's value order to the sorted order
        vorder: Vec<usize>,
        vals: Vec<f64>,
        sym: symbolic::Symbolic,
        lsym: lu::LuSymbolic,
        sched: (Vec<usize>, Vec<usize>),
        num: Option<lu::LuNumeric>,
        /// partial-pivoting fallback (used when static pivoting is not accurate enough)
        fallback: Option<faer::sparse::linalg::solvers::Lu<usize, f64>>,
        use_fallback: bool,
    },
}

pub struct Solver {
    pool: rayon::ThreadPool,
    pub nthreads: usize,
    pub ordering: Ordering,
    pub pivot_eps: f64,
    /// pivot_eps was set explicitly (otherwise LU uses 1e-13, LDL^T 1e-8, as PARDISO)
    pub pivot_eps_set: bool,
    pub par_flops: f64,
    /// 0: refine only if pivots were perturbed (max `max_refine` steps); n>0: always n steps
    pub refine_mode: i32,
    pub max_refine: usize,
    /// reuse the previous fill-reducing ordering when the pattern changes only slightly
    pub reuse_ordering: bool,
    kind: Kind,
    pub stats: Stats,
}

impl Solver {
    pub fn new(nthreads: usize) -> Solver {
        let nthreads = if nthreads == 0 {
            std::thread::available_parallelism()
                .map(|v| v.get())
                .unwrap_or(1)
        } else {
            nthreads
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(nthreads)
            .build()
            .expect("thread pool");
        Solver {
            pool,
            nthreads,
            ordering: Ordering::Metis,
            pivot_eps: 1e-8,
            pivot_eps_set: false,
            par_flops: 4.0e6,
            refine_mode: 0,
            max_refine: 2,
            reuse_ordering: true,
            kind: Kind::None,
            stats: Stats::default(),
        }
    }

    pub fn symbolic(&self) -> Option<&symbolic::Symbolic> {
        match &self.kind {
            Kind::Symmetric { sym, .. } => Some(sym),
            _ => None,
        }
    }

    /// Symmetric matrix, lower triangle in CSC (equivalently upper triangle in CSR), 0- or 1-based.
    pub fn analyze_symmetric(
        &mut self,
        n: usize,
        colptr: &[i32],
        rowidx: &[i32],
        offset: i32,
    ) -> Result<(), String> {
        let t0 = Instant::now();
        let (cp, ri) = to_usize(n, colptr, rowidx, offset)?;
        // validate: lower triangle
        for j in 0..n {
            for &i in &ri[cp[j]..cp[j + 1]] {
                if i < j || i >= n {
                    return Err(format!("entry ({i},{j}) is not in the lower triangle"));
                }
            }
        }
        // FEBio calls PreProcess whenever it rebuilds the matrix profile (e.g. contact
        // updates). Identical pattern -> keep the analysis; similar -> keep the ordering.
        let mut hint: Option<Vec<usize>> = None;
        if let Kind::Symmetric {
            colptr,
            rowidx,
            sym,
            ..
        } = &self.kind
        {
            if colptr == &cp && rowidx == &ri {
                self.stats.n_reused += 1;
                self.stats.t_analyze += t0.elapsed().as_secs_f64();
                return Ok(());
            }
            if self.reuse_ordering && sym.n == n && similar_nnz(rowidx.len(), ri.len()) {
                hint = Some(sym.perm.clone());
            }
        }
        let ordering = self.ordering;
        let sym = self.pool.install(|| {
            let a = symbolic::LowerCsc {
                n,
                colptr: &cp,
                rowidx: &ri,
            };
            symbolic::analyze_with_hint(&a, ordering, hint.as_deref())
        })?;
        self.stats.reset_pattern();
        self.stats.n = n;
        self.stats.nnz_a = ri.len();
        self.stats.nnz_l = sym.nnz_l;
        self.stats.n_supernodes = sym.ns;
        self.stats.n_supervars = sym.n_supervars;
        self.stats.flops = sym.flops;
        let sched = numeric::schedule(&sym, self.nthreads);
        self.kind = Kind::Symmetric {
            colptr: cp,
            rowidx: ri,
            sym,
            sched,
            num: None,
            vals: Vec::new(),
        };
        self.stats.t_analyze += t0.elapsed().as_secs_f64();
        Ok(())
    }

    /// General (structurally symmetric or not) matrix in CSR, 0- or 1-based.
    pub fn analyze_unsymmetric(
        &mut self,
        n: usize,
        rowptr: &[i32],
        colidx: &[i32],
        offset: i32,
    ) -> Result<(), String> {
        let t0 = Instant::now();
        let (rp, ci0) = to_usize(n, rowptr, colidx, offset)?;
        // sort column indices per row, remembering the value permutation
        let mut ci = vec![0usize; ci0.len()];
        let mut vorder = vec![0usize; ci0.len()];
        for r in 0..n {
            let mut idx: Vec<usize> = (rp[r]..rp[r + 1]).collect();
            idx.sort_unstable_by_key(|&t| ci0[t]);
            for (q, &t) in idx.iter().enumerate() {
                ci[rp[r] + q] = ci0[t];
                vorder[rp[r] + q] = t;
            }
        }
        let mut hint: Option<Vec<usize>> = None;
        if let Kind::Unsymmetric {
            rowptr,
            colidx,
            vorder: vo,
            sym,
            ..
        } = &self.kind
        {
            if rowptr == &rp && colidx == &ci && vo == &vorder {
                self.stats.n_reused += 1;
                self.stats.t_analyze += t0.elapsed().as_secs_f64();
                return Ok(());
            }
            if self.reuse_ordering && sym.n == n && similar_nnz(colidx.len(), ci.len()) {
                hint = Some(sym.perm.clone());
            }
        }
        let ordering = self.ordering;
        let (sym, lsym) = self.pool.install(|| -> Result<_, String> {
            let (lcp, lri) = lu::symmetrized_lower(n, &rp, &ci);
            let a = symbolic::LowerCsc {
                n,
                colptr: &lcp,
                rowidx: &lri,
            };
            let sym = symbolic::analyze_with_hint(&a, ordering, hint.as_deref())?;
            let lsym = lu::analyze(&sym, &rp, &ci)?;
            Ok((sym, lsym))
        })?;
        self.stats.reset_pattern();
        self.stats.n = n;
        self.stats.nnz_a = ci.len();
        self.stats.nnz_l = 2 * sym.nnz_l;
        self.stats.n_supernodes = sym.ns;
        self.stats.n_supervars = sym.n_supervars;
        self.stats.flops = 2.0 * sym.flops;
        let sched = numeric::schedule(&sym, self.nthreads);
        let use_fallback = std::env::var("FRS_UNSYM")
            .map(|v| v == "faer")
            .unwrap_or(false);
        self.kind = Kind::Unsymmetric {
            n,
            rowptr: rp,
            colidx: ci,
            vorder,
            vals: Vec::new(),
            sym,
            lsym,
            sched,
            num: None,
            fallback: None,
            use_fallback,
        };
        self.stats.t_analyze += t0.elapsed().as_secs_f64();
        Ok(())
    }

    pub fn factor(&mut self, vals: &[f64]) -> Result<(), String> {
        let t0 = Instant::now();
        let opt = numeric::FactorOptions {
            nthreads: self.nthreads,
            pivot_eps: self.pivot_eps,
            par_flops: self.par_flops,
        };
        match &mut self.kind {
            Kind::None => return Err("factor called before analyze".into()),
            Kind::Symmetric {
                colptr,
                rowidx,
                sym,
                sched,
                num,
                vals: v,
            } => {
                if vals.len() != rowidx.len() {
                    return Err("value array has wrong length".into());
                }
                v.clear();
                v.extend_from_slice(vals);
                if let Ok(path) = std::env::var("FRS_DUMP") {
                    dump_lower(&path, colptr, rowidx, vals);
                }
                let old = num.take();
                let r = self
                    .pool
                    .install(|| numeric::factor(sym, vals, &opt, old, sched))?;
                self.stats.n_perturbed = r.n_perturbed;
                *num = Some(r);
            }
            Kind::Unsymmetric {
                n,
                rowptr,
                colidx,
                vorder,
                vals: v,
                sym,
                lsym,
                sched,
                num,
                fallback,
                use_fallback,
            } => {
                if vals.len() != colidx.len() {
                    return Err("value array has wrong length".into());
                }
                v.clear();
                v.extend(vorder.iter().map(|&t| vals[t]));
                *fallback = None;
                if *use_fallback {
                    *fallback = Some(self.pool.install(|| faer_lu(*n, rowptr, colidx, v))?);
                } else {
                    let old = num.take();
                    let nt = self.nthreads;
                    // PARDISO's default perturbation for unsymmetric matrices is 1e-13 (iparm[9]=13)
                    let eps = if self.pivot_eps_set {
                        opt.pivot_eps
                    } else {
                        1e-13
                    };
                    let r = self
                        .pool
                        .install(|| lu::factor(sym, lsym, v, sched, nt, eps, opt.par_flops, old));
                    self.stats.n_perturbed = r.n_perturbed;
                    *num = Some(r);
                }
            }
        }
        self.stats.n_factor += 1;
        self.stats.t_factor += t0.elapsed().as_secs_f64();
        Ok(())
    }

    /// solve A x = b
    pub fn solve(&mut self, x: &mut [f64], b: &[f64]) -> Result<(), String> {
        let t0 = Instant::now();
        match &self.kind {
            Kind::None => return Err("solve called before analyze".into()),
            Kind::Symmetric {
                colptr,
                rowidx,
                sym,
                sched,
                num,
                vals,
            } => {
                let num = num.as_ref().ok_or("solve called before factor")?;
                let n = sym.n;
                let solve_once = |rhs: &[f64], out: &mut [f64], w: &mut Vec<f64>| {
                    w.clear();
                    w.extend(sym.perm.iter().map(|&p| rhs[p]));
                    let nt = self.nthreads;
                    self.pool.install(|| {
                        numeric::solve_permuted_par(sym, num, w, &sched.0, &sched.1, nt)
                    });
                    for (i, &p) in sym.perm.iter().enumerate() {
                        out[p] = w[i];
                    }
                };
                let mut w = Vec::with_capacity(n);
                solve_once(b, x, &mut w);
                let steps = if self.refine_mode > 0 {
                    self.refine_mode as usize
                } else if num.n_perturbed > 0 {
                    self.max_refine
                } else {
                    0
                };
                if steps > 0 {
                    let mut r = vec![0.0; n];
                    let mut dx = vec![0.0; n];
                    for _ in 0..steps {
                        sym_spmv_lower(n, colptr, rowidx, vals, x, &mut r);
                        for i in 0..n {
                            r[i] = b[i] - r[i];
                        }
                        solve_once(&r, &mut dx, &mut w);
                        for i in 0..n {
                            x[i] += dx[i];
                        }
                        self.stats.n_refine += 1;
                    }
                }
            }
            Kind::Unsymmetric { .. } => self.solve_unsym(x, b)?,
        }
        self.stats.n_solve += 1;
        self.stats.t_solve += t0.elapsed().as_secs_f64();
        Ok(())
    }
}

impl Solver {
    fn solve_unsym(&mut self, x: &mut [f64], b: &[f64]) -> Result<(), String> {
        use faer::linalg::solvers::SolveCore;
        use faer::MatMut;
        let nt = self.nthreads;
        let max_refine = self.max_refine.max(3);
        let forced = if self.refine_mode > 0 {
            self.refine_mode as usize
        } else {
            0
        };
        let Kind::Unsymmetric {
            n,
            rowptr,
            colidx,
            vals,
            sym,
            lsym,
            num,
            fallback,
            ..
        } = &mut self.kind
        else {
            unreachable!()
        };
        let n = *n;
        let fb_solve = |lu: &faer::sparse::linalg::solvers::Lu<usize, f64>,
                        x: &mut [f64],
                        b: &[f64],
                        pool: &rayon::ThreadPool| {
            x[..n].copy_from_slice(&b[..n]);
            let xm = MatMut::from_column_major_slice_mut(&mut x[..n], n, 1);
            // the CSR of A was handed to faer as CSC, i.e. A^T was factored
            pool.install(|| lu.solve_transpose_in_place_with_conj(faer::Conj::No, xm));
        };
        if let Some(lu) = fallback.as_ref() {
            fb_solve(lu, x, b, &self.pool);
            return Ok(());
        }
        let numr = num.as_ref().ok_or("solve called before factor")?;
        let mut w = Vec::with_capacity(n);
        let solve_once = |rhs: &[f64], out: &mut [f64], w: &mut Vec<f64>| {
            w.clear();
            w.extend(sym.perm.iter().map(|&p| rhs[p]));
            lu::solve_permuted(sym, lsym, numr, w);
            for (i, &p) in sym.perm.iter().enumerate() {
                out[p] = w[i];
            }
        };
        solve_once(b, x, &mut w);
        let bn = b
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
            .max(f64::MIN_POSITIVE);
        let mut r = vec![0.0; n];
        let mut dx = vec![0.0; n];
        let resid = |x: &[f64], r: &mut [f64]| -> f64 {
            lu::csr_spmv(n, rowptr, colidx, vals, x, r);
            let mut s = 0.0;
            for i in 0..n {
                r[i] = b[i] - r[i];
                s += r[i] * r[i];
            }
            s.sqrt() / bn
        };
        let _ = nt;
        let mut rel = resid(x, &mut r);
        let mut steps = 0;
        while (steps < forced) || (rel > 1e-12 && steps < max_refine) {
            solve_once(&r, &mut dx, &mut w);
            for i in 0..n {
                x[i] += dx[i];
            }
            steps += 1;
            self.stats.n_refine += 1;
            rel = resid(x, &mut r);
        }
        if std::env::var_os("FRS_PROFILE").is_some() {
            eprintln!("febio_rs_solver: LU solve rel residual {rel:.2e} after {steps} refinement steps, perturbed pivots {}", numr.n_perturbed);
        }
        if rel > 1e-8 {
            // static pivoting was not accurate enough: partial-pivoting LU for this matrix
            eprintln!("febio_rs_solver: static-pivot LU residual {rel:.2e}; switching to partial pivoting");
            let lu = self.pool.install(|| faer_lu(n, rowptr, colidx, vals))?;
            fb_solve(&lu, x, b, &self.pool);
            *fallback = Some(lu);
            self.stats.n_fallback += 1;
        }
        Ok(())
    }
}

/// partial-pivoting sparse LU (faer) of the CSR matrix, factored as CSC of A^T
fn faer_lu(
    n: usize,
    rp: &[usize],
    ci: &[usize],
    v: &[f64],
) -> Result<faer::sparse::linalg::solvers::Lu<usize, f64>, String> {
    use faer::sparse::linalg::solvers::{Lu, SymbolicLu};
    use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
    let s = SymbolicSparseColMatRef::new_checked(n, n, rp, None, ci);
    let slu = SymbolicLu::try_new(s).map_err(|e| format!("symbolic LU: {e:?}"))?;
    Lu::try_new_with_symbolic(slu, SparseColMatRef::new(s, v)).map_err(|e| format!("LU: {e:?}"))
}

/// patterns whose nonzero counts differ by < 5% are treated as similar (ordering reuse)
fn similar_nnz(a: usize, b: usize) -> bool {
    let (a, b) = (a as f64, b as f64);
    (a - b).abs() <= 0.05 * a.max(b)
}

/// debugging aid: write the matrix as  n, nnz (i64), colptr (i32), rowidx (i32), values (f64)
fn dump_lower(path: &str, cp: &[usize], ri: &[usize], v: &[f64]) {
    use std::io::Write;
    if let Ok(f) = std::fs::File::create(path) {
        let mut f = std::io::BufWriter::new(f);
        let n = cp.len() - 1;
        let _ = f.write_all(&(n as i64).to_le_bytes());
        let _ = f.write_all(&(ri.len() as i64).to_le_bytes());
        for &x in cp {
            let _ = f.write_all(&(x as i32).to_le_bytes());
        }
        for &x in ri {
            let _ = f.write_all(&(x as i32).to_le_bytes());
        }
        for &x in v {
            let _ = f.write_all(&x.to_le_bytes());
        }
    }
}

fn to_usize(
    n: usize,
    ptr: &[i32],
    idx: &[i32],
    offset: i32,
) -> Result<(Vec<usize>, Vec<usize>), String> {
    if ptr.len() != n + 1 {
        return Err("pointer array has wrong length".into());
    }
    let cp: Vec<usize> = ptr.iter().map(|&p| (p - offset) as usize).collect();
    if cp[0] != 0 || cp[n] != idx.len() {
        return Err("inconsistent pointer array".into());
    }
    let ri: Vec<usize> = idx.iter().map(|&i| (i - offset) as usize).collect();
    Ok((cp, ri))
}

/// y = A x for a symmetric matrix stored as lower CSC
pub fn sym_spmv_lower(n: usize, cp: &[usize], ri: &[usize], v: &[f64], x: &[f64], y: &mut [f64]) {
    y[..n].fill(0.0);
    for j in 0..n {
        let xj = x[j];
        let mut acc = 0.0;
        for t in cp[j]..cp[j + 1] {
            let i = ri[t];
            let a = v[t];
            if i == j {
                acc += a * xj;
            } else {
                y[i] += a * xj;
                acc += a * x[i];
            }
        }
        y[j] += acc;
    }
}

// ============================================================================
// C ABI
// ============================================================================

#[repr(C)]
pub struct FrsStats {
    pub n: i64,
    pub nnz_a: i64,
    pub nnz_l: i64,
    pub n_supernodes: i64,
    pub n_supervars: i64,
    pub flops: f64,
    pub t_analyze: f64,
    pub t_factor: f64,
    pub t_solve: f64,
    pub n_factor: i64,
    pub n_solve: i64,
    pub n_perturbed: i64,
    pub n_refine: i64,
    pub n_analyze: i64,
    pub n_reused: i64,
    pub n_fallback: i64,
}

fn report(e: String) -> i32 {
    eprintln!("febio_rs_solver: {e}");
    -1
}

/// Create a solver. nthreads = 0 -> all cores.
#[no_mangle]
pub extern "C" fn frs_create(nthreads: i32) -> *mut Solver {
    Box::into_raw(Box::new(Solver::new(nthreads.max(0) as usize)))
}

#[no_mangle]
pub unsafe extern "C" fn frs_destroy(s: *mut Solver) {
    if !s.is_null() {
        drop(Box::from_raw(s));
    }
}

/// option keys: 1 = ordering (0 metis, 1 amd, 2 natural), 2 = pivot eps,
/// 3 = refinement mode (0 auto, n>0 always n steps), 4 = max refinement steps (auto mode),
/// 5 = parallel dense-kernel flop threshold, 6 = reuse ordering for similar patterns (0/1)
#[no_mangle]
pub unsafe extern "C" fn frs_set_option(s: *mut Solver, key: i32, val: f64) -> i32 {
    let s = &mut *s;
    match key {
        1 => {
            s.ordering = match val as i32 {
                1 => Ordering::Amd,
                2 => Ordering::Natural,
                _ => Ordering::Metis,
            }
        }
        2 => {
            s.pivot_eps = val;
            s.pivot_eps_set = true;
        }
        3 => s.refine_mode = val as i32,
        4 => s.max_refine = val.max(0.0) as usize,
        5 => s.par_flops = val,
        6 => s.reuse_ordering = val != 0.0,
        _ => return -1,
    }
    0
}

/// symmetric: lower-triangular CSC == upper-triangular CSR (the PARDISO mtype=-2 layout)
#[no_mangle]
pub unsafe extern "C" fn frs_analyze_symmetric(
    s: *mut Solver,
    n: i32,
    ptr: *const i32,
    idx: *const i32,
    offset: i32,
) -> i32 {
    let s = &mut *s;
    let n = n as usize;
    let p = std::slice::from_raw_parts(ptr, n + 1);
    let nnz = (p[n] - offset) as usize;
    let i = std::slice::from_raw_parts(idx, nnz);
    match s.analyze_symmetric(n, p, i, offset) {
        Ok(()) => 0,
        Err(e) => report(e),
    }
}

/// unsymmetric: full CSR
#[no_mangle]
pub unsafe extern "C" fn frs_analyze_unsymmetric(
    s: *mut Solver,
    n: i32,
    ptr: *const i32,
    idx: *const i32,
    offset: i32,
) -> i32 {
    let s = &mut *s;
    let n = n as usize;
    let p = std::slice::from_raw_parts(ptr, n + 1);
    let nnz = (p[n] - offset) as usize;
    let i = std::slice::from_raw_parts(idx, nnz);
    match s.analyze_unsymmetric(n, p, i, offset) {
        Ok(()) => 0,
        Err(e) => report(e),
    }
}

#[no_mangle]
pub unsafe extern "C" fn frs_factor(s: *mut Solver, vals: *const f64, nnz: i64) -> i32 {
    let s = &mut *s;
    let v = std::slice::from_raw_parts(vals, nnz as usize);
    match s.factor(v) {
        Ok(()) => 0,
        Err(e) => report(e),
    }
}

#[no_mangle]
pub unsafe extern "C" fn frs_solve(s: *mut Solver, x: *mut f64, b: *const f64) -> i32 {
    let s = &mut *s;
    let n = s.stats.n;
    let x = std::slice::from_raw_parts_mut(x, n);
    let b = std::slice::from_raw_parts(b, n).to_vec();
    match s.solve(x, &b) {
        Ok(()) => 0,
        Err(e) => report(e),
    }
}

/// free the numeric factors (keeps the symbolic analysis for pattern reuse)
#[no_mangle]
pub unsafe extern "C" fn frs_release_numeric(s: *mut Solver) {
    let s = &mut *s;
    match &mut s.kind {
        Kind::Symmetric { num, vals, .. } => {
            *num = None;
            *vals = Vec::new();
        }
        Kind::Unsymmetric {
            num,
            fallback,
            vals,
            ..
        } => {
            *num = None;
            *fallback = None;
            *vals = Vec::new();
        }
        Kind::None => {}
    }
}

#[no_mangle]
pub unsafe extern "C" fn frs_get_stats(s: *const Solver, out: *mut FrsStats) {
    let st = &(*s).stats;
    *out = FrsStats {
        n: st.n as i64,
        nnz_a: st.nnz_a as i64,
        nnz_l: st.nnz_l as i64,
        n_supernodes: st.n_supernodes as i64,
        n_supervars: st.n_supervars as i64,
        flops: st.flops,
        t_analyze: st.t_analyze,
        t_factor: st.t_factor,
        t_solve: st.t_solve,
        n_factor: st.n_factor as i64,
        n_solve: st.n_solve as i64,
        n_perturbed: st.n_perturbed as i64,
        n_refine: st.n_refine as i64,
        n_analyze: st.n_analyze as i64,
        n_reused: st.n_reused as i64,
        n_fallback: st.n_fallback as i64,
    };
}
