//! Numeric supernodal LDL^T factorization (left-looking, tree-parallel) and solves.
//!
//! Panel of supernode s: m_s x k_s column-major block holding
//!   rows 0..k   : the dense diagonal block (L11 strictly lower, D on the diagonal)
//!   rows k..m   : L21 (rows given by the supernode's off-diagonal pattern)
//!
//! Parallelism:
//!   * independent subtrees of the supernodal elimination tree are factored concurrently
//!     (rayon work stealing),
//!   * large supernodes (near the root, where tree parallelism runs out) use faer's
//!     multithreaded dense kernels (gemm / trsm / LDL^T).
//! The dependency structure guarantees that a supernode only reads panels of its own
//! (already finished) descendants and only writes its own panel.

use crate::symbolic::Symbolic;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::ldlt::factor::{
    cholesky_in_place, cholesky_in_place_scratch, LdltParams, LdltRegularization,
};
use faer::linalg::matmul::matmul;
use faer::linalg::triangular_solve::{
    solve_unit_lower_triangular_in_place, solve_unit_upper_triangular_in_place,
};
use faer::reborrow::*;
use faer::{Accum, MatMut, MatRef, Par, Spec};
use std::cell::RefCell;
use std::sync::atomic::{AtomicUsize, Ordering as AO};

#[derive(Clone, Copy)]
struct SyncPtr(*mut f64);
unsafe impl Send for SyncPtr {}
unsafe impl Sync for SyncPtr {}
impl SyncPtr {
    #[inline(always)]
    fn get(self) -> *mut f64 {
        self.0
    }
}

pub struct Numeric {
    pub l: Vec<f64>,
    pub d: Vec<f64>,
    pub n_perturbed: usize,
}

pub struct FactorOptions {
    pub nthreads: usize,
    /// pivots with |d| <= pivot_eps * max|A_ij| are replaced by +-pivot_eps * max|A_ij|
    /// (same strategy as PARDISO's pivot perturbation, default 1e-8)
    pub pivot_eps: f64,
    /// supernodes with at least this many flops use parallel dense kernels
    pub par_flops: f64,
}

/// FRS_PROFILE=1 prints per-phase / per-supernode timings (read once)
fn profiling() -> bool {
    static P: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *P.get_or_init(|| std::env::var_os("FRS_PROFILE").is_some())
}

thread_local! {
    static WS: RefCell<Workspace> = RefCell::new(Workspace::default());
}

#[derive(Default)]
struct Workspace {
    map: Vec<u32>,
    t: Vec<f64>,
    c: Vec<f64>,
}

#[inline]
unsafe fn panel_mut<'a>(base: SyncPtr, off: usize, m: usize, k: usize) -> MatMut<'a, f64> {
    MatMut::from_raw_parts_mut(base.0.add(off), m, k, 1, m as isize)
}
#[inline]
unsafe fn panel_ref<'a>(base: *const f64, off: usize, m: usize, k: usize) -> MatRef<'a, f64> {
    MatRef::from_raw_parts(base.add(off), m, k, 1, m as isize)
}

/// Geist-Ng style tree partition: returns (independent subtree roots sorted by decreasing
/// work, remaining top-of-tree supernodes in postorder)
pub fn schedule(sym: &Symbolic, nthreads: usize) -> (Vec<usize>, Vec<usize>) {
    let ns = sym.ns;
    if nthreads <= 1 {
        return (sym.roots.clone(), vec![]);
    }
    let w = &sym.subtree_flops;
    let mut set: Vec<usize> = sym.roots.clone();
    let mut is_top = vec![false; ns];
    let total: f64 = set.iter().map(|&r| w[r]).sum();
    for _ in 0..ns {
        // LPT greedy assignment to estimate balance
        set.sort_by(|a, b| w[*b].partial_cmp(&w[*a]).unwrap());
        let mut bins = vec![0.0f64; nthreads];
        for &r in &set {
            let (imin, _) =
                bins.iter().enumerate().fold(
                    (0, f64::MAX),
                    |acc, (i, &v)| if v < acc.1 { (i, v) } else { acc },
                );
            bins[imin] += w[r];
        }
        let maxb = bins.iter().cloned().fold(0.0, f64::max);
        let sum: f64 = bins.iter().sum();
        let balanced = set.len() >= 2 * nthreads && maxb <= 1.10 * sum / nthreads as f64;
        // stop if balanced, or if the remaining forest is a small fraction of the work
        if balanced || sum < 0.02 * total || set.is_empty() {
            break;
        }
        // split the heaviest subtree
        let r = set.remove(0);
        is_top[r] = true;
        set.extend_from_slice(sym.children(r));
    }
    let top: Vec<usize> = (0..ns).filter(|&s| is_top[s]).collect();
    set.sort_by(|a, b| w[*b].partial_cmp(&w[*a]).unwrap());
    (set, top)
}

/// Allocate a large f64 buffer without touching it, advising the kernel to back it with
/// transparent huge pages (fewer, cheaper page faults; first-touch faults serialize
/// between threads). Each supernode zeroes its own panel inside its task.
fn alloc_uninit_huge(len: usize) -> Vec<f64> {
    let mut v: Vec<f64> = Vec::with_capacity(len);
    #[cfg(target_os = "linux")]
    unsafe {
        let bytes = len * std::mem::size_of::<f64>();
        if bytes >= (4 << 20) {
            let p = v.as_mut_ptr() as usize;
            let page = 2usize << 20;
            let start = (p + page - 1) & !(page - 1);
            let end = (p + bytes) & !(page - 1);
            if end > start {
                libc::madvise(start as *mut libc::c_void, end - start, libc::MADV_HUGEPAGE);
            }
        }
    }
    // SAFETY: every element is written (zeroed then assembled) before it is read
    unsafe { v.set_len(len) };
    v
}

pub fn factor(
    sym: &Symbolic,
    avals: &[f64],
    opt: &FactorOptions,
    reuse: Option<Numeric>,
    sched: &(Vec<usize>, Vec<usize>),
) -> Result<Numeric, String> {
    let n = sym.n;
    let len = sym.len_val();
    let (mut l, mut d) = match reuse {
        Some(old) if old.l.len() == len && old.d.len() == n => (old.l, old.d),
        _ => (alloc_uninit_huge(len), vec![0.0f64; n]),
    };
    let amax = avals.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    let eps = if amax > 0.0 {
        opt.pivot_eps * amax
    } else {
        opt.pivot_eps
    };
    let nperturbed = AtomicUsize::new(0);
    let err = std::sync::Mutex::new(None::<String>);
    let lp = SyncPtr(l.as_mut_ptr());
    let dp = SyncPtr(d.as_mut_ptr());

    let ctx = Ctx {
        sym,
        avals,
        lp,
        dp,
        eps,
        opt,
        nperturbed: &nperturbed,
        err: &err,
    };
    let (subtrees, top) = (&sched.0, &sched.1);
    {
        use rayon::prelude::*;
        // phase 1: independent subtrees, each factored sequentially (largest first)
        let t0 = std::time::Instant::now();
        let timed = |r: usize| {
            let t = std::time::Instant::now();
            ctx.subtree_seq(r);
            if profiling() {
                eprintln!(
                    "  subtree {r} ({:.2e} flops) {:.3}s on {:?}",
                    sym.subtree_flops[r],
                    t.elapsed().as_secs_f64(),
                    rayon::current_thread_index()
                );
            }
        };
        subtrees.par_iter().with_max_len(1).for_each(|&r| timed(r));
        let t1 = t0.elapsed().as_secs_f64();
        // phase 2: the top of the tree, in postorder, with parallel dense kernels
        // siblings in the top part run concurrently (phase 1 is complete, so no long
        // unrelated subtree work can be stolen into a critical-path supernode)
        if top.len() > 1 {
            let mut is_top = vec![false; sym.ns];
            for &s in top.iter() {
                is_top[s] = true;
            }
            let is_top = &is_top;
            fn rec(ctx: &Ctx, s: usize, is_top: &[bool]) {
                let ch: Vec<usize> = ctx
                    .sym
                    .children(s)
                    .iter()
                    .copied()
                    .filter(|&c| is_top[c])
                    .collect();
                if ch.len() > 1 {
                    rayon::scope(|sc| {
                        for &c in &ch {
                            sc.spawn(move |_| rec(ctx, c, is_top));
                        }
                    });
                } else if let Some(&c) = ch.first() {
                    rec(ctx, c, is_top);
                }
                ctx.supernode(s, true);
            }
            let roots: Vec<usize> = top
                .iter()
                .copied()
                .filter(|&s| sym.sn_parent[s] == usize::MAX || !is_top[sym.sn_parent[s]])
                .collect();
            rayon::scope(|sc| {
                for &r in &roots {
                    let ctx = &ctx;
                    sc.spawn(move |_| rec(ctx, r, is_top));
                }
            });
        } else {
            for &s in top.iter() {
                ctx.supernode(s, true);
            }
        }
        if profiling() {
            let w1: f64 = subtrees.iter().map(|&r| sym.subtree_flops[r]).sum();
            eprintln!("schedule: {} subtrees ({:.2e} flops, max {:.2e}) {:.3}s ; top {} supernodes ({:.2e} flops) {:.3}s",
                subtrees.len(), w1, subtrees.first().map(|&r| sym.subtree_flops[r]).unwrap_or(0.0), t1, top.len(), sym.flops - w1, t0.elapsed().as_secs_f64() - t1);
        }
    }
    if let Some(e) = err.into_inner().unwrap() {
        return Err(e);
    }
    Ok(Numeric {
        l,
        d,
        n_perturbed: nperturbed.load(AO::Relaxed),
    })
}

struct Ctx<'a> {
    sym: &'a Symbolic,
    avals: &'a [f64],
    lp: SyncPtr,
    dp: SyncPtr,
    eps: f64,
    opt: &'a FactorOptions,
    nperturbed: &'a AtomicUsize,
    err: &'a std::sync::Mutex<Option<String>>,
}

impl<'a> Ctx<'a> {
    fn subtree_seq(&self, r: usize) {
        // supernodes of the subtree rooted at r are the contiguous postorder range
        // [first_desc(r), r]
        let first = self.sym.first_desc[r];
        for s in first..=r {
            self.supernode(s, false);
        }
    }

    /// apply all descendant updates to columns [c0, c1) of supernode s
    #[allow(clippy::too_many_arguments)]
    fn apply_updates(
        &self,
        s: usize,
        c0: usize,
        c1: usize,
        map: &[u32],
        t: &mut Vec<f64>,
        cbuf: &mut Vec<f64>,
        par: Par,
    ) {
        let sym = self.sym;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let base = self.lp;
        let raw = unsafe { std::slice::from_raw_parts_mut(base.0.add(sym.val_ptr[s]), m * k) };
        for u in sym.updates(s) {
            let dsn = u.d as usize;
            let kd = sym.k(dsn);
            let md = sym.m(dsn);
            let pd = sym.pattern(dsn);
            let (r0, r1) = (u.r0 as usize, u.r1 as usize);
            // restrict to target columns [b+c0, b+c1)
            let p1 = &pd[r0..r1];
            let lo = r0 + p1.partition_point(|&r| r < b + c0);
            let hi = r0 + p1.partition_point(|&r| r < b + c1);
            if lo >= hi {
                continue;
            }
            let len1 = hi - lo;
            // rows always start at r0 (not at the block's first column lo): the gemm below then
            // has the same row alignment however the target columns are split into blocks
            let nrow = pd.len() - r0;
            let ld = unsafe { panel_ref(base.0, sym.val_ptr[dsn], md, kd) };
            let dd = unsafe { std::slice::from_raw_parts(self.dp.0.add(sym.sn_begin[dsn]), kd) };
            let ld1 = ld.subrows(kd + lo, len1);
            let ldr = ld.subrows(kd + r0, nrow);
            // T = Ld1 * D   (len1 x kd)
            t.clear();
            t.resize(len1 * kd, 0.0);
            for j in 0..kd {
                let dj = dd[j];
                let src = ld1.col(j);
                let dst = &mut t[j * len1..(j + 1) * len1];
                for i in 0..len1 {
                    dst[i] = src[i] * dj;
                }
            }
            let tm = MatRef::from_column_major_slice(t, len1, kd);
            // C = L_d[rows r0..] * T^T  (nrow x len1). One plain gemm: every entry is produced by
            // the same kernel with the same alignment no matter how the target columns are
            // split, which keeps the factorization bitwise independent of the thread count.
            // Only the lower trapezoid (rows >= the target column) is used.
            cbuf.clear();
            cbuf.resize(nrow * len1, 0.0);
            {
                let cm = MatMut::from_column_major_slice_mut(cbuf, nrow, len1);
                matmul(cm, Accum::Replace, ldr, tm.transpose(), 1.0, par);
            }
            // scatter-subtract into the panel of s
            let prow = &pd[r0..];
            let off = lo - r0;
            for cc in 0..len1 {
                let tcol = prow[off + cc] - b;
                let dst = &mut raw[tcol * m..(tcol + 1) * m];
                let src = &cbuf[cc * nrow..(cc + 1) * nrow];
                for rr in (off + cc)..nrow {
                    let tr = map[prow[rr]] as usize;
                    dst[tr] -= src[rr];
                }
            }
        }
    }

    fn supernode(&self, s: usize, top: bool) {
        let t_start = std::time::Instant::now();
        self.supernode_inner(s, top);
        if profiling() {
            let el = t_start.elapsed().as_secs_f64();
            if el > 0.01 {
                eprintln!(
                    "sn {s}: k={} m={} nupd={} sub_flops={:.2e} t={:.4}s",
                    self.sym.k(s),
                    self.sym.m(s),
                    self.sym.updates(s).len(),
                    self.sym.subtree_flops[s],
                    el
                );
            }
        }
    }

    fn supernode_inner(&self, s: usize, top: bool) {
        let t_in = std::time::Instant::now();
        let sym = self.sym;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        let off = sym.val_ptr[s];
        let sflops = {
            let mut f = 0.0;
            for c in 0..k {
                let below = (m - c - 1) as f64;
                f += below * below;
            }
            f
        };
        // intra-node parallelism only in the top part of the tree (no tree parallelism left)
        let par = if top && sflops >= self.opt.par_flops && self.opt.nthreads > 1 {
            Par::rayon(self.opt.nthreads)
        } else {
            Par::Seq
        };
        // SAFETY: panel of s is written only here; descendants are complete (tree order)
        let base = self.lp;
        let mut ls = unsafe { panel_mut(base, off, m, k) };
        // 1. zero the panel and assemble A
        {
            let raw = unsafe { std::slice::from_raw_parts_mut(base.0.add(off), m * k) };
            raw.fill(0.0);
            for &(t, pos) in &sym.amap[sym.amap_ptr[s]..sym.amap_ptr[s + 1]] {
                raw[pos] += self.avals[t];
            }
        }
        // 2. left-looking updates from descendants
        let upds = sym.updates(s);
        if !upds.is_empty() {
            let mut uflops = 0.0f64;
            for u in upds {
                let kd = sym.k(u.d as usize);
                let len1 = (u.r1 - u.r0) as f64;
                let len2 = (sym.pattern(u.d as usize).len() - u.r0 as usize) as f64;
                uflops += len1 * len2 * kd as f64;
            }
            if top && uflops >= 4.0 * self.opt.par_flops && self.opt.nthreads > 1 && k >= 64 {
                // parallel over column blocks of the target panel: race free since
                // every block only writes its own columns
                let mut map = vec![0u32; sym.n];
                for c in 0..k {
                    map[b + c] = c as u32;
                }
                for (r, &row) in pat.iter().enumerate() {
                    map[row] = (k + r) as u32;
                }
                let bw = 32usize.max(k / (8 * self.opt.nthreads)).min(256);
                let nb = (k + bw - 1) / bw;
                use rayon::prelude::*;
                let map = &map;
                (0..nb).into_par_iter().for_each(|blk| {
                    let c0 = blk * bw;
                    let c1 = (c0 + bw).min(k);
                    let mut t = Vec::new();
                    let mut c = Vec::new();
                    self.apply_updates(s, c0, c1, map, &mut t, &mut c, Par::Seq);
                });
            } else {
                WS.with(|ws| {
                    let ws = &mut *ws.borrow_mut();
                    if ws.map.len() < sym.n {
                        ws.map.resize(sym.n, 0);
                    }
                    for c in 0..k {
                        ws.map[b + c] = c as u32;
                    }
                    for (r, &row) in pat.iter().enumerate() {
                        ws.map[row] = (k + r) as u32;
                    }
                    let Workspace { map, t, c } = ws;
                    self.apply_updates(s, 0, k, map, t, c, Par::Seq);
                });
            }
        }
        let t_upd = t_in.elapsed().as_secs_f64();
        // 3. dense LDL^T of the diagonal block with pivot perturbation
        {
            let (l11, _) = ls.rb_mut().split_at_row_mut(k);
            let params: Spec<LdltParams, f64> = Default::default();
            let req = cholesky_in_place_scratch::<f64>(k, par, params);
            let mut mem = MemBuffer::new(req);
            let reg = LdltRegularization {
                dynamic_regularization_signs: None,
                dynamic_regularization_delta: self.eps,
                dynamic_regularization_epsilon: self.eps,
            };
            // count pivots that will be perturbed (approximate: checked after the fact)
            match cholesky_in_place(l11, reg, par, MemStack::new(&mut mem), params) {
                Ok(_) => {}
                Err(e) => {
                    *self.err.lock().unwrap() =
                        Some(format!("LDLT failed in supernode {s}: {e:?}"));
                    return;
                }
            }
        }
        let dd = unsafe { std::slice::from_raw_parts_mut(self.dp.0.add(b), k) };
        let mut np = 0;
        for c in 0..k {
            let v = ls[(c, c)];
            if v.abs() <= self.eps * (1.0 + 1e-12) {
                np += 1;
            }
            dd[c] = v;
        }
        if np > 0 {
            self.nperturbed.fetch_add(np, AO::Relaxed);
        }
        // 4. L21 = A21 * L11^{-T} * D^{-1}
        if m > k {
            let (l11, l21) = ls.rb_mut().split_at_row_mut(k);
            let l11 = l11.rb();
            solve_unit_lower_triangular_in_place(l11, l21.transpose_mut(), par);
            let mut l21 = ls.rb_mut().subrows_mut(k, m - k);
            for c in 0..k {
                let inv = 1.0 / dd[c];
                for x in l21.rb_mut().col_mut(c).iter_mut() {
                    *x *= inv;
                }
            }
        }
        if profiling() && t_in.elapsed().as_secs_f64() > 0.05 {
            let mut uflops = 0.0f64;
            let mut skd = 0.0;
            for u in sym.updates(s) {
                let kd = sym.k(u.d as usize);
                let len1 = (u.r1 - u.r0) as f64;
                let len2 = (sym.pattern(u.d as usize).len() - u.r0 as usize) as f64;
                uflops += len1 * len2 * kd as f64;
                skd += len1 * len2;
            }
            eprintln!(
                "   sn {s}: upd {:.4}s ({:.2e} flop, {:.1} GF/s, avg kd {:.1}), dense {:.4}s",
                t_upd,
                uflops,
                uflops / t_upd * 1e-9,
                uflops / skd,
                t_in.elapsed().as_secs_f64() - t_upd
            );
        }
    }
}

/// x <- L^{-1} D^{-1} L^{-T} x  in the permuted ordering (in place)
pub fn solve_permuted(sym: &Symbolic, num: &Numeric, x: &mut [f64]) {
    let ns = sym.ns;
    let lp = num.l.as_ptr();
    // forward: L y = b
    for s in 0..ns {
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        let ls = unsafe { panel_ref(lp, sym.val_ptr[s], m, k) };
        {
            let xs = MatMut::from_column_major_slice_mut(&mut x[b..b + k], k, 1);
            solve_unit_lower_triangular_in_place(ls.subrows(0, k), xs, Par::Seq);
        }
        if m > k {
            let l21 = ls.subrows(k, m - k);
            for c in 0..k {
                let xc = x[b + c];
                if xc == 0.0 {
                    continue;
                }
                let col = l21.col(c);
                for (r, &row) in pat.iter().enumerate() {
                    x[row] -= col[r] * xc;
                }
            }
        }
    }
    // diagonal
    for i in 0..sym.n {
        x[i] /= num.d[i];
    }
    // backward: L^T x = y
    for s in (0..ns).rev() {
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        let ls = unsafe { panel_ref(lp, sym.val_ptr[s], m, k) };
        if m > k {
            let l21 = ls.subrows(k, m - k);
            for c in 0..k {
                let col = l21.col(c);
                let mut acc = 0.0;
                for (r, &row) in pat.iter().enumerate() {
                    acc += col[r] * x[row];
                }
                x[b + c] -= acc;
            }
        }
        let xs = MatMut::from_column_major_slice_mut(&mut x[b..b + k], k, 1);
        solve_unit_upper_triangular_in_place(ls.subrows(0, k).transpose(), xs, Par::Seq);
    }
}

/// Tree-parallel version of `solve_permuted`.
///
/// forward  (L y = b): left-looking; supernode s gathers the contributions of its
///          descendants into its own entries only -> independent subtrees run in parallel,
///          then the top of the tree (with parallel gathers for large supernodes);
/// backward (L^T x = y): top of the tree first (parallel over columns of large supernodes),
///          then the independent subtrees in parallel; each supernode only writes its own
///          entries and reads ancestors, which are final.
pub fn solve_permuted_par(
    sym: &Symbolic,
    num: &Numeric,
    x: &mut [f64],
    subtrees: &[usize],
    top: &[usize],
    nthreads: usize,
) {
    use rayon::prelude::*;
    // Always the same (left-looking, chunked) algorithm, also for one thread, so that the
    // solution is bitwise identical for any thread count.
    let _ = nthreads;
    let xp = SyncPtr(x.as_mut_ptr());
    let lp = num.l.as_ptr() as usize;
    let n = sym.n;

    // All accesses to x go through the raw pointer: concurrently running supernodes write
    // disjoint index ranges (their own columns) and only read entries that are final.
    let fwd = |s: usize, par: bool| {
        let xq = xp;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let upds = sym.updates(s);
        let gather = |u: &crate::symbolic::Update, acc: &mut [f64]| {
            let dsn = u.d as usize;
            let kd = sym.k(dsn);
            let md = sym.m(dsn);
            let bd = sym.sn_begin[dsn];
            let pd = sym.pattern(dsn);
            let (r0, r1) = (u.r0 as usize, u.r1 as usize);
            let ld = unsafe { panel_ref(lp as *const f64, sym.val_ptr[dsn], md, kd) };
            for j in 0..kd {
                let yj = unsafe { *xq.get().add(bd + j) };
                if yj == 0.0 {
                    continue;
                }
                let col = ld.col(j);
                for t in r0..r1 {
                    acc[pd[t] - b] -= col[kd + t] * yj;
                }
            }
        };
        // Deterministic gather: the updates are split into fixed-size chunks (independent of
        // the thread count), each chunk is summed in update order, and the chunk sums are
        // added in chunk order. The serial and parallel paths perform exactly the same
        // floating-point operations, so the result is bitwise reproducible.
        const CHUNK: usize = 32;
        let nchunks = upds.len().div_ceil(CHUNK);
        let chunk_sum = |ci: usize| -> Vec<f64> {
            let mut acc = vec![0.0f64; k];
            for u in &upds[ci * CHUNK..((ci + 1) * CHUNK).min(upds.len())] {
                gather(u, &mut acc);
            }
            acc
        };
        let partial: Vec<Vec<f64>> = if par && nchunks > 1 {
            (0..nchunks).into_par_iter().map(chunk_sum).collect()
        } else {
            (0..nchunks).map(chunk_sum).collect()
        };
        let mut acc = vec![0.0f64; k];
        for pc in &partial {
            for (a, v) in acc.iter_mut().zip(pc) {
                *a += v;
            }
        }
        // SAFETY: x[b..b+k] is owned by supernode s during this call
        let xs = unsafe { std::slice::from_raw_parts_mut(xq.get().add(b), k) };
        for c in 0..k {
            xs[c] += acc[c];
        }
        let ls = unsafe { panel_ref(lp as *const f64, sym.val_ptr[s], sym.m(s), k) };
        let xm = MatMut::from_column_major_slice_mut(xs, k, 1);
        solve_unit_lower_triangular_in_place(ls.subrows(0, k), xm, Par::Seq);
    };

    let bwd = |s: usize, par: bool| {
        let xq = xp;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        let ls = unsafe { panel_ref(lp as *const f64, sym.val_ptr[s], m, k) };
        let mut v = vec![0.0f64; k];
        if m > k {
            let l21 = ls.subrows(k, m - k);
            let dot = |c: usize| -> f64 {
                let col = l21.col(c);
                let mut acc = 0.0;
                for (r, &row) in pat.iter().enumerate() {
                    // rows of the pattern belong to ancestors (final values)
                    acc += col[r] * unsafe { *xq.get().add(row) };
                }
                acc
            };
            if par && (m - k) * k >= 200_000 {
                v.par_iter_mut()
                    .enumerate()
                    .for_each(|(c, out)| *out = dot(c));
            } else {
                for c in 0..k {
                    v[c] = dot(c);
                }
            }
        }
        // SAFETY: x[b..b+k] is owned by supernode s during this call
        let xs = unsafe { std::slice::from_raw_parts_mut(xq.get().add(b), k) };
        for c in 0..k {
            xs[c] -= v[c];
        }
        let xm = MatMut::from_column_major_slice_mut(xs, k, 1);
        solve_unit_upper_triangular_in_place(ls.subrows(0, k).transpose(), xm, Par::Seq);
    };

    // forward
    subtrees.par_iter().with_max_len(1).for_each(|&r| {
        for s in sym.first_desc[r]..=r {
            fwd(s, false);
        }
    });
    for &s in top {
        fwd(s, true);
    }
    // diagonal
    {
        let x = unsafe { std::slice::from_raw_parts_mut(xp.0, n) };
        x.par_chunks_mut(16384).enumerate().for_each(|(ci, ch)| {
            let o = ci * 16384;
            for (i, v) in ch.iter_mut().enumerate() {
                *v /= num.d[o + i];
            }
        });
    }
    // backward
    for &s in top.iter().rev() {
        bwd(s, true);
    }
    subtrees.par_iter().with_max_len(1).for_each(|&r| {
        for s in (sym.first_desc[r]..=r).rev() {
            bwd(s, false);
        }
    });
}
