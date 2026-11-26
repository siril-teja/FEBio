//! febio_rs_solver: parallel sparse direct solver for FEBio written in Rust.
//!
//! * symmetric (indefinite allowed) matrices: supernodal LDL^T with METIS nested
//!   dissection ordering, tree + dense-kernel parallelism, PARDISO-style pivot
//!   perturbation and optional iterative refinement.
//! * unsymmetric matrices: faer's supernodal sparse LU (partial pivoting).
//!
//! C ABI (see `febio_rs_solver.h`).

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
}

impl Default for Stats {
    fn default() -> Self {
        Stats { n: 0, nnz_a: 0, nnz_l: 0, n_supernodes: 0, n_supervars: 0, flops: 0.0, t_analyze: 0.0, t_factor: 0.0, t_solve: 0.0, n_factor: 0, n_solve: 0, n_perturbed: 0, n_refine: 0 }
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
        // matrix stored as CSC of A^T (= CSR of A)
        colptr: Vec<usize>,
        rowidx: Vec<usize>,
        lu: Option<faer::sparse::linalg::solvers::Lu<usize, f64>>,
        slu: faer::sparse::linalg::solvers::SymbolicLu<usize>,
    },
}

pub struct Solver {
    pool: rayon::ThreadPool,
    pub nthreads: usize,
    pub ordering: Ordering,
    pub pivot_eps: f64,
    pub par_flops: f64,
    /// 0: refine only if pivots were perturbed (max `max_refine` steps); n>0: always n steps
    pub refine_mode: i32,
    pub max_refine: usize,
    kind: Kind,
    pub stats: Stats,
}

impl Solver {
    pub fn new(nthreads: usize) -> Solver {
        let nthreads = if nthreads == 0 {
            std::thread::available_parallelism().map(|v| v.get()).unwrap_or(1)
        } else {
            nthreads
        };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(nthreads).build().expect("thread pool");
        Solver {
            pool,
            nthreads,
            ordering: Ordering::Metis,
            pivot_eps: 1e-8,
            par_flops: 4.0e6,
            refine_mode: 0,
            max_refine: 2,
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
    pub fn analyze_symmetric(&mut self, n: usize, colptr: &[i32], rowidx: &[i32], offset: i32) -> Result<(), String> {
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
        let ordering = self.ordering;
        let sym = self.pool.install(|| {
            let a = symbolic::LowerCsc { n, colptr: &cp, rowidx: &ri };
            symbolic::analyze(&a, ordering)
        })?;
        self.stats = Stats::default();
        self.stats.n = n;
        self.stats.nnz_a = ri.len();
        self.stats.nnz_l = sym.nnz_l;
        self.stats.n_supernodes = sym.ns;
        self.stats.n_supervars = sym.n_supervars;
        self.stats.flops = sym.flops;
        let sched = numeric::schedule(&sym, self.nthreads);
        self.kind = Kind::Symmetric { colptr: cp, rowidx: ri, sym, sched, num: None, vals: Vec::new() };
        self.stats.t_analyze += t0.elapsed().as_secs_f64();
        Ok(())
    }

    /// General matrix in CSR (row pointers / column indices), 0- or 1-based.
    pub fn analyze_unsymmetric(&mut self, n: usize, rowptr: &[i32], colidx: &[i32], offset: i32) -> Result<(), String> {
        use faer::sparse::linalg::solvers::SymbolicLu;
        use faer::sparse::SymbolicSparseColMatRef;
        let t0 = Instant::now();
        let (cp, ri) = to_usize(n, rowptr, colidx, offset)?;
        // faer wants sorted indices per column
        let mut ri = ri;
        for j in 0..n {
            ri[cp[j]..cp[j + 1]].sort_unstable();
        }
        let slu = {
            let a = SymbolicSparseColMatRef::new_checked(n, n, &cp, None, &ri);
            self.pool.install(|| SymbolicLu::try_new(a)).map_err(|e| format!("symbolic LU: {e:?}"))?
        };
        self.stats = Stats::default();
        self.stats.n = n;
        self.stats.nnz_a = ri.len();
        self.kind = Kind::Unsymmetric { n, colptr: cp, rowidx: ri, lu: None, slu };
        self.stats.t_analyze += t0.elapsed().as_secs_f64();
        Ok(())
    }

    pub fn factor(&mut self, vals: &[f64]) -> Result<(), String> {
        let t0 = Instant::now();
        let opt = numeric::FactorOptions { nthreads: self.nthreads, pivot_eps: self.pivot_eps, par_flops: self.par_flops };
        match &mut self.kind {
            Kind::None => return Err("factor called before analyze".into()),
            Kind::Symmetric { colptr, rowidx, sym, sched, num, vals: v } => {
                if vals.len() != rowidx.len() {
                    return Err("value array has wrong length".into());
                }
                v.clear();
                v.extend_from_slice(vals);
                if let Ok(path) = std::env::var("FRS_DUMP") {
                    dump_lower(&path, colptr, rowidx, vals);
                }
                let old = num.take();
                let r = self.pool.install(|| numeric::factor(sym, vals, &opt, old, sched))?;
                self.stats.n_perturbed = r.n_perturbed;
                *num = Some(r);
            }
            Kind::Unsymmetric { n, colptr, rowidx, lu, slu } => {
                use faer::sparse::linalg::solvers::Lu;
                use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
                // the stored matrix is A^T in CSC; fix reordering of values matching sorted indices
                let s = SymbolicSparseColMatRef::new_checked(*n, *n, colptr, None, rowidx);
                if vals.len() != rowidx.len() {
                    return Err("value array has wrong length".into());
                }
                let a = SparseColMatRef::new(s, vals);
                let r = self.pool.install(|| Lu::try_new_with_symbolic(slu.clone(), a)).map_err(|e| format!("LU: {e:?}"))?;
                *lu = Some(r);
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
            Kind::Symmetric { colptr, rowidx, sym, sched, num, vals } => {
                let num = num.as_ref().ok_or("solve called before factor")?;
                let n = sym.n;
                let solve_once = |rhs: &[f64], out: &mut [f64], w: &mut Vec<f64>| {
                    w.clear();
                    w.extend(sym.perm.iter().map(|&p| rhs[p]));
                    let nt = self.nthreads;
                    self.pool.install(|| numeric::solve_permuted_par(sym, num, w, &sched.0, &sched.1, nt));
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
            Kind::Unsymmetric { n, lu, .. } => {
                use faer::linalg::solvers::SolveCore;
                use faer::MatMut;
                let lu = lu.as_ref().ok_or("solve called before factor")?;
                x[..*n].copy_from_slice(&b[..*n]);
                let xm = MatMut::from_column_major_slice_mut(&mut x[..*n], *n, 1);
                // we factored A^T (CSR of A viewed as CSC) -> A x = b is a transpose solve
                self.pool.install(|| lu.solve_transpose_in_place_with_conj(faer::Conj::No, xm));
            }
        }
        self.stats.n_solve += 1;
        self.stats.t_solve += t0.elapsed().as_secs_f64();
        Ok(())
    }
}

/// debugging aid: write the matrix as  n, nnz (i64), colptr (i32), rowidx (i32), values (f64)
fn dump_lower(path: &str, cp: &[usize], ri: &[usize], v: &[f64]) {
    use std::io::Write;
    if let Ok(f) = std::fs::File::create(path) {
        let mut f = std::io::BufWriter::new(f);
        let n = cp.len() - 1;
        let _ = f.write_all(&(n as i64).to_le_bytes());
        let _ = f.write_all(&(ri.len() as i64).to_le_bytes());
        for &x in cp { let _ = f.write_all(&(x as i32).to_le_bytes()); }
        for &x in ri { let _ = f.write_all(&(x as i32).to_le_bytes()); }
        for &x in v { let _ = f.write_all(&x.to_le_bytes()); }
    }
}

fn to_usize(n: usize, ptr: &[i32], idx: &[i32], offset: i32) -> Result<(Vec<usize>, Vec<usize>), String> {
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
/// 5 = parallel dense-kernel flop threshold
#[no_mangle]
pub unsafe extern "C" fn frs_set_option(s: *mut Solver, key: i32, val: f64) -> i32 {
    let s = &mut *s;
    match key {
        1 => s.ordering = match val as i32 { 1 => Ordering::Amd, 2 => Ordering::Natural, _ => Ordering::Metis },
        2 => s.pivot_eps = val,
        3 => s.refine_mode = val as i32,
        4 => s.max_refine = val.max(0.0) as usize,
        5 => s.par_flops = val,
        _ => return -1,
    }
    0
}

/// symmetric: lower-triangular CSC == upper-triangular CSR (the PARDISO mtype=-2 layout)
#[no_mangle]
pub unsafe extern "C" fn frs_analyze_symmetric(s: *mut Solver, n: i32, ptr: *const i32, idx: *const i32, offset: i32) -> i32 {
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
pub unsafe extern "C" fn frs_analyze_unsymmetric(s: *mut Solver, n: i32, ptr: *const i32, idx: *const i32, offset: i32) -> i32 {
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
    };
}
