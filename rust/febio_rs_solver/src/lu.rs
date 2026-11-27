//! Supernodal LU for structurally symmetric, numerically unsymmetric matrices
//! (FEBio's REAL_UNSYMMETRIC / REAL_SYMM_STRUCTURE systems: biphasic, multiphasic,
//! fluid, non-symmetric contact, ...).
//!
//! The symbolic structure of |A| + |A|^T is shared with the LDL^T solver (METIS ordering,
//! supernodes, tree schedule). Factorization A = L U uses static pivoting (no row
//! interchanges) with PARDISO-style perturbation of tiny pivots; accuracy is safeguarded by
//! iterative refinement and a residual check in the caller (which falls back to a
//! partial-pivoting LU if needed).
//!
//! Storage per supernode s (k columns, m = k + |pattern| rows):
//!   L panel  m x k (col-major): top k x k = L11 (strict lower, unit) and U11 (upper incl.
//!            diagonal) in place; rows k..m = L21
//!   Ut panel (m-k) x k (col-major): U12^T (row r = pattern row, column c = U row b+c)
//!
//! Left-looking update of target s from descendant d (rows pd[r0..r1) in s's columns):
//!   L_s[rows pd[r0..], cols pd[lo..hi)] -= L_d[pd[r0..]] * U_d[:, pd[lo..hi)]
//!   U_s[rows pd[lo..hi), cols pd[r1..]] -= L_d[pd[lo..hi)] * U_d[:, pd[r1..]]

use crate::symbolic::Symbolic;
use faer::linalg::matmul::matmul;
use faer::linalg::triangular_solve::{
    solve_lower_triangular_in_place, solve_unit_lower_triangular_in_place,
    solve_upper_triangular_in_place,
};
use faer::reborrow::*;
use faer::{Accum, MatMut, MatRef, Par};
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

/// assembly target of an entry of A
#[derive(Clone, Copy)]
pub struct Target {
    /// index into the value array of A (CSR order)
    pub t: usize,
    /// absolute position: < len_l -> L storage, otherwise Ut storage at pos - len_l
    pub pos: usize,
}

pub struct LuSymbolic {
    /// offsets of Ut panels
    pub ut_ptr: Vec<usize>,
    pub amap_ptr: Vec<usize>,
    pub amap: Vec<Target>,
}

pub struct LuNumeric {
    pub l: Vec<f64>,
    pub ut: Vec<f64>,
    pub n_perturbed: usize,
}

/// build the lower-CSC pattern of |A| + |A|^T (incl. full diagonal) from a CSR pattern
pub fn symmetrized_lower(n: usize, rp: &[usize], ci: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let mut cnt = vec![1usize; n]; // diagonal
    for r in 0..n {
        for &c in &ci[rp[r]..rp[r + 1]] {
            if r != c {
                cnt[r.min(c)] += 1;
            }
        }
    }
    let mut cp = vec![0usize; n + 1];
    for j in 0..n {
        cp[j + 1] = cp[j] + cnt[j];
    }
    let mut pos = cp.clone();
    let mut ri = vec![0usize; cp[n]];
    for j in 0..n {
        ri[pos[j]] = j;
        pos[j] += 1;
    }
    for r in 0..n {
        for &c in &ci[rp[r]..rp[r + 1]] {
            if r != c {
                let (lo, hi) = (r.min(c), r.max(c));
                ri[pos[lo]] = hi;
                pos[lo] += 1;
            }
        }
    }
    // sort + dedupe each column
    let mut out_cp = vec![0usize; n + 1];
    let mut out_ri = Vec::with_capacity(ri.len());
    for j in 0..n {
        let col = &mut ri[cp[j]..cp[j + 1]];
        col.sort_unstable();
        let start = out_ri.len();
        for &i in col.iter() {
            if out_ri.len() == start || *out_ri.last().unwrap() != i {
                out_ri.push(i);
            }
        }
        out_cp[j + 1] = out_ri.len();
    }
    (out_cp, out_ri)
}

pub fn analyze(sym: &Symbolic, rp: &[usize], ci: &[usize]) -> Result<LuSymbolic, String> {
    let ns = sym.ns;
    let n = sym.n;
    let mut sn_of = vec![0usize; n];
    for s in 0..ns {
        for j in sym.sn_begin[s]..sym.sn_begin[s + 1] {
            sn_of[j] = s;
        }
    }
    let mut ut_ptr = vec![0usize; ns + 1];
    for s in 0..ns {
        ut_ptr[s + 1] = ut_ptr[s] + (sym.m(s) - sym.k(s)) * sym.k(s);
    }
    let len_l = sym.len_val();
    let mut cnt = vec![0usize; ns + 1];
    for r in 0..n {
        for &c in &ci[rp[r]..rp[r + 1]] {
            cnt[sn_of[sym.iperm[r].min(sym.iperm[c])]] += 1;
        }
    }
    let mut amap_ptr = vec![0usize; ns + 1];
    for s in 0..ns {
        amap_ptr[s + 1] = amap_ptr[s] + cnt[s];
    }
    let mut fill = amap_ptr.clone();
    let mut amap = vec![Target { t: 0, pos: 0 }; amap_ptr[ns]];
    let local_row = |s: usize, row: usize| -> Option<usize> {
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        if row >= b && row < b + k {
            Some(row - b)
        } else {
            sym.pattern(s).binary_search(&row).ok().map(|x| k + x)
        }
    };
    for r in 0..n {
        for t in rp[r]..rp[r + 1] {
            let c = ci[t];
            let (pr, pc) = (sym.iperm[r], sym.iperm[c]);
            let s = sn_of[pr.min(pc)];
            let b = sym.sn_begin[s];
            let k = sym.k(s);
            let m = sym.m(s);
            let pos = if pr >= pc {
                // lower (incl. diagonal): L panel, column pc
                let lr = local_row(s, pr).ok_or("entry outside symbolic pattern")?;
                sym.val_ptr[s] + (pc - b) * m + lr
            } else {
                // upper: row pr in s's columns
                let lc = local_row(s, pc).ok_or("entry outside symbolic pattern")?;
                if lc < k {
                    sym.val_ptr[s] + lc * m + (pr - b) // diagonal block, upper part
                } else {
                    len_l + ut_ptr[s] + (pr - b) * (m - k) + (lc - k)
                }
            };
            amap[fill[s]] = Target { t, pos };
            fill[s] += 1;
        }
    }
    Ok(LuSymbolic {
        ut_ptr,
        amap_ptr,
        amap,
    })
}

struct Ctx<'a> {
    sym: &'a Symbolic,
    lsym: &'a LuSymbolic,
    avals: &'a [f64],
    lp: SyncPtr,
    up: SyncPtr,
    len_l: usize,
    eps: f64,
    nthreads: usize,
    par_flops: f64,
    nperturbed: &'a AtomicUsize,
}

#[inline]
unsafe fn mat_mut<'a>(base: SyncPtr, off: usize, m: usize, k: usize) -> MatMut<'a, f64> {
    MatMut::from_raw_parts_mut(base.get().add(off), m, k, 1, m as isize)
}
#[inline]
unsafe fn mat_ref<'a>(base: SyncPtr, off: usize, m: usize, k: usize) -> MatRef<'a, f64> {
    MatRef::from_raw_parts(base.get().add(off), m, k, 1, m as isize)
}

impl<'a> Ctx<'a> {
    fn subtree_seq(&self, r: usize) {
        for s in self.sym.first_desc[r]..=r {
            self.supernode(s, false);
        }
    }

    /// apply all descendant updates to target columns [c0, c1) of supernode s
    fn apply_updates(
        &self,
        s: usize,
        c0: usize,
        c1: usize,
        map: &[u32],
        cl: &mut Vec<f64>,
        cu: &mut Vec<f64>,
    ) {
        let sym = self.sym;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let ls =
            unsafe { std::slice::from_raw_parts_mut(self.lp.get().add(sym.val_ptr[s]), m * k) };
        let us = unsafe {
            std::slice::from_raw_parts_mut(self.up.get().add(self.lsym.ut_ptr[s]), (m - k) * k)
        };
        for u in sym.updates(s) {
            let d = u.d as usize;
            let kd = sym.k(d);
            let md = sym.m(d);
            let pd = sym.pattern(d);
            let (r0, r1) = (u.r0 as usize, u.r1 as usize);
            let p1 = &pd[r0..r1];
            let lo = r0 + p1.partition_point(|&r| r < b + c0);
            let hi = r0 + p1.partition_point(|&r| r < b + c1);
            if lo >= hi {
                continue;
            }
            let ncol = hi - lo;
            let ld = unsafe { mat_ref(self.lp, sym.val_ptr[d], md, kd) };
            let utd = unsafe { mat_ref(self.up, self.lsym.ut_ptr[d], md - kd, kd) };
            // C_L = L_d[rows r0..] * U_d[:, cols lo..hi]      ((|pd|-r0) x ncol)
            let nl = pd.len() - r0;
            cl.clear();
            cl.resize(nl * ncol, 0.0);
            {
                let c = MatMut::from_column_major_slice_mut(cl, nl, ncol);
                matmul(
                    c,
                    Accum::Replace,
                    ld.subrows(kd + r0, nl),
                    utd.subrows(lo, ncol).transpose(),
                    1.0,
                    Par::Seq,
                );
            }
            for cc in 0..ncol {
                let tcol = pd[lo + cc] - b;
                let dst = &mut ls[tcol * m..(tcol + 1) * m];
                let src = &cl[cc * nl..(cc + 1) * nl];
                for rr in 0..nl {
                    dst[map[pd[r0 + rr]] as usize] -= src[rr];
                }
            }
            // C_U = Ut_d[rows r1..] * L_d[rows lo..hi]^T     ((|pd|-r1) x ncol)
            let nu = pd.len() - r1;
            if nu > 0 {
                cu.clear();
                cu.resize(nu * ncol, 0.0);
                {
                    let c = MatMut::from_column_major_slice_mut(cu, nu, ncol);
                    matmul(
                        c,
                        Accum::Replace,
                        utd.subrows(r1, nu),
                        ld.subrows(kd + lo, ncol).transpose(),
                        1.0,
                        Par::Seq,
                    );
                }
                for cc in 0..ncol {
                    let tcol = pd[lo + cc] - b;
                    let dst = &mut us[tcol * (m - k)..(tcol + 1) * (m - k)];
                    let src = &cu[cc * nu..(cc + 1) * nu];
                    for rr in 0..nu {
                        dst[map[pd[r1 + rr]] as usize - k] -= src[rr];
                    }
                }
            }
        }
    }

    fn supernode(&self, s: usize, top: bool) {
        let sym = self.sym;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        // 1. zero + assemble
        unsafe {
            std::slice::from_raw_parts_mut(self.lp.get().add(sym.val_ptr[s]), m * k).fill(0.0);
            std::slice::from_raw_parts_mut(self.up.get().add(self.lsym.ut_ptr[s]), (m - k) * k)
                .fill(0.0);
        }
        for tg in &self.lsym.amap[self.lsym.amap_ptr[s]..self.lsym.amap_ptr[s + 1]] {
            unsafe {
                if tg.pos < self.len_l {
                    *self.lp.get().add(tg.pos) += self.avals[tg.t];
                } else {
                    *self.up.get().add(tg.pos - self.len_l) += self.avals[tg.t];
                }
            }
        }
        // 2. descendant updates
        let upds = sym.updates(s);
        if !upds.is_empty() {
            let mut map = vec![0u32; 0];
            let mut uflops = 0.0;
            for u in upds {
                let d = u.d as usize;
                uflops += 2.0
                    * ((u.r1 - u.r0) as f64)
                    * ((sym.pattern(d).len() - u.r0 as usize) as f64)
                    * sym.k(d) as f64;
            }
            thread_local! { static MAP: std::cell::RefCell<Vec<u32>> = const { std::cell::RefCell::new(Vec::new()) }; }
            let par = top && self.nthreads > 1 && k >= 64 && uflops >= 4.0 * self.par_flops;
            if par {
                map.resize(sym.n, 0);
                for c in 0..k {
                    map[b + c] = c as u32;
                }
                for (r, &row) in pat.iter().enumerate() {
                    map[row] = (k + r) as u32;
                }
                let bw = 32usize.max(k / (8 * self.nthreads)).min(256);
                let nb = k.div_ceil(bw);
                use rayon::prelude::*;
                let map = &map;
                (0..nb).into_par_iter().for_each(|blk| {
                    let (mut cl, mut cu) = (Vec::new(), Vec::new());
                    self.apply_updates(s, blk * bw, ((blk + 1) * bw).min(k), map, &mut cl, &mut cu);
                });
            } else {
                MAP.with(|mp| {
                    let mut mp = mp.borrow_mut();
                    if mp.len() < sym.n {
                        mp.resize(sym.n, 0);
                    }
                    for c in 0..k {
                        mp[b + c] = c as u32;
                    }
                    for (r, &row) in pat.iter().enumerate() {
                        mp[row] = (k + r) as u32;
                    }
                    let (mut cl, mut cu) = (Vec::new(), Vec::new());
                    self.apply_updates(s, 0, k, &mp, &mut cl, &mut cu);
                });
            }
        }
        // 3. dense LU (no pivoting) of the k x k diagonal block, then off-diagonal panels
        let mut flops_own = 0.0;
        for c in 0..k {
            flops_own += ((m - c - 1) as f64).powi(2);
        }
        let dpar = if top && self.nthreads > 1 && flops_own >= self.par_flops {
            Par::rayon(self.nthreads)
        } else {
            Par::Seq
        };
        let mut lpan = unsafe { mat_mut(self.lp, sym.val_ptr[s], m, k) };
        let np = lu_nopiv(lpan.rb_mut().subrows_mut(0, k), self.eps, dpar);
        if np > 0 {
            self.nperturbed.fetch_add(np, AO::Relaxed);
        }
        if m > k {
            let (a11, a21) = lpan.rb_mut().split_at_row_mut(k);
            let a11 = a11.rb();
            // L21 = A21 U11^{-1}   <=>  U11^T L21^T = A21^T
            solve_lower_triangular_in_place(a11.transpose(), a21.transpose_mut(), dpar);
            // U12 = L11^{-1} A12, stored transposed: Ut12^T
            let ut = unsafe { mat_mut(self.up, self.lsym.ut_ptr[s], m - k, k) };
            solve_unit_lower_triangular_in_place(a11, ut.transpose_mut(), dpar);
        }
    }
}

/// recursive LU without pivoting, in place (L unit lower, U upper); pivots with
/// |p| <= eps are replaced by +-eps. Returns the number of perturbed pivots.
fn lu_nopiv(a: MatMut<'_, f64>, eps: f64, par: Par) -> usize {
    let n = a.nrows();
    if n <= 32 {
        return lu_small(a, eps);
    }
    let n1 = n / 2;
    let (left, right) = a.split_at_col_mut(n1);
    let (mut a11, mut a21) = left.split_at_row_mut(n1);
    let (mut a12, a22) = right.split_at_row_mut(n1);
    let mut np = lu_nopiv(a11.rb_mut(), eps, par);
    // A12 <- L11^{-1} A12 ; A21 <- A21 U11^{-1}
    solve_unit_lower_triangular_in_place(a11.rb(), a12.rb_mut(), par);
    solve_lower_triangular_in_place(a11.rb().transpose(), a21.rb_mut().transpose_mut(), par);
    let mut a22 = a22;
    matmul(a22.rb_mut(), Accum::Add, a21.rb(), a12.rb(), -1.0, par);
    np += lu_nopiv(a22, eps, par);
    np
}

fn lu_small(mut a: MatMut<'_, f64>, eps: f64) -> usize {
    let n = a.nrows();
    let mut np = 0;
    for j in 0..n {
        let mut p = a[(j, j)];
        if p.abs() <= eps {
            p = if p < 0.0 { -eps } else { eps };
            a[(j, j)] = p;
            np += 1;
        }
        let inv = 1.0 / p;
        for i in j + 1..n {
            a[(i, j)] *= inv;
        }
        for jj in j + 1..n {
            let u = a[(j, jj)];
            if u != 0.0 {
                for i in j + 1..n {
                    let l = a[(i, j)];
                    a[(i, jj)] -= l * u;
                }
            }
        }
    }
    np
}

/// Run the factorization tasks as a dependency graph: every subtree (a root in
/// `subtrees`, processed sequentially by `sub`) and every top supernode (`top`, processed
/// by `topf`) starts as soon as everything it depends on is done. A top supernode
/// depends on its children, each of which is a top supernode or a subtree root. Each
/// task's arithmetic does not depend on when it runs, so the result is deterministic.
/// Falls back to "all subtrees, then the top supernodes in order" if the schedule does
/// not have that shape.
pub(crate) fn run_tree<FS, FT>(sym: &Symbolic, subtrees: &[usize], top: &[usize], sub: FS, topf: FT)
where
    FS: Fn(usize) + Sync,
    FT: Fn(usize) + Sync,
{
    use rayon::prelude::*;
    use std::sync::atomic::AtomicU32;
    const NONE: u32 = u32::MAX;
    // slot of each top supernode
    let mut slot = vec![NONE; sym.ns];
    for (q, &s) in top.iter().enumerate() {
        slot[s] = q as u32;
    }
    let mut is_sub = vec![false; sym.ns];
    for &r in subtrees {
        is_sub[r] = true;
    }
    // pending children of each top supernode; every child must be top or a subtree root
    let mut pending = vec![0u32; top.len()];
    let mut ok = true;
    for (q, &s) in top.iter().enumerate() {
        for &c in &sym.sn_children[sym.sn_children_ptr[s]..sym.sn_children_ptr[s + 1]] {
            if slot[c] == NONE && !is_sub[c] {
                ok = false;
            }
            pending[q] += 1;
        }
    }
    for &r in subtrees {
        let p = sym.sn_parent[r];
        if p != usize::MAX && slot[p] == NONE {
            ok = false;
        }
    }
    if !ok {
        subtrees.par_iter().with_max_len(1).for_each(|&r| sub(r));
        for &s in top {
            topf(s);
        }
        return;
    }
    let pending: Vec<AtomicU32> = pending.into_iter().map(AtomicU32::new).collect();
    // finishing a node may make its parent ready
    fn done<'s, FT: Fn(usize) + Sync>(
        sc: &rayon::Scope<'s>,
        sym: &'s Symbolic,
        slot: &'s [u32],
        pending: &'s [AtomicU32],
        topf: &'s FT,
        node: usize,
    ) {
        let p = sym.sn_parent[node];
        if p == usize::MAX {
            return;
        }
        let q = slot[p] as usize;
        if pending[q].fetch_sub(1, AO::AcqRel) == 1 {
            sc.spawn(move |sc| {
                topf(p);
                done(sc, sym, slot, pending, topf, p);
            });
        }
    }
    let (slot, pending, topf, sub) = (&slot, &pending, &topf, &sub);
    rayon::scope(|sc| {
        // top supernodes without children (rare) are ready at once
        for (q, &s) in top.iter().enumerate() {
            if pending[q].load(AO::Acquire) == 0 {
                sc.spawn(move |sc| {
                    topf(s);
                    done(sc, sym, slot, pending, topf, s);
                });
            }
        }
        for &r in subtrees {
            sc.spawn(move |sc| {
                sub(r);
                done(sc, sym, slot, pending, topf, r);
            });
        }
    });
}

pub fn factor(
    sym: &Symbolic,
    lsym: &LuSymbolic,
    avals: &[f64],
    sched: &(Vec<usize>, Vec<usize>),
    nthreads: usize,
    pivot_eps: f64,
    par_flops: f64,
    reuse: Option<LuNumeric>,
) -> LuNumeric {
    let len_l = sym.len_val();
    let len_u = lsym.ut_ptr[sym.ns];
    let (mut l, mut ut) = match reuse {
        Some(o) if o.l.len() == len_l && o.ut.len() == len_u => (o.l, o.ut),
        _ => (vec![0.0; len_l], vec![0.0; len_u]),
    };
    let amax = avals.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
    let eps = pivot_eps * if amax > 0.0 { amax } else { 1.0 };
    let np = AtomicUsize::new(0);
    let ctx = Ctx {
        sym,
        lsym,
        avals,
        lp: SyncPtr(l.as_mut_ptr()),
        up: SyncPtr(ut.as_mut_ptr()),
        len_l,
        eps,
        nthreads,
        par_flops,
        nperturbed: &np,
    };
    let prof = std::env::var_os("FRS_PROFILE").is_some();
    let t0 = std::time::Instant::now();
    run_tree(
        sym,
        &sched.0,
        &sched.1,
        |r| ctx.subtree_seq(r),
        |s| ctx.supernode(s, true),
    );
    if prof {
        eprintln!(
            "lu factor: {} subtrees + {} top supernodes (task graph) {:.3}s",
            sched.0.len(),
            sched.1.len(),
            t0.elapsed().as_secs_f64()
        );
    }
    LuNumeric {
        l,
        ut,
        n_perturbed: np.load(AO::Relaxed),
    }
}

/// x <- U^{-1} L^{-1} x in permuted ordering, parallel over the elimination tree.
///
/// forward  (L y = b): left-looking; supernode s gathers its descendants' contributions
///          into its own entries only, in fixed-size chunks summed in a fixed order, so
///          independent subtrees run in parallel and the result does not depend on the
///          thread count;
/// backward (U x = y): top of the tree first (parallel over the columns of large
///          supernodes), then the independent subtrees in parallel; each supernode only
///          writes its own entries and reads ancestors, which are final.
/// The same algorithm runs for any thread count: the solution is bitwise identical.
pub fn solve_permuted_par(
    sym: &Symbolic,
    lsym: &LuSymbolic,
    num: &LuNumeric,
    x: &mut [f64],
    subtrees: &[usize],
    top: &[usize],
) {
    use rayon::prelude::*;
    let xp = SyncPtr(x.as_mut_ptr());
    let lp = SyncPtr(num.l.as_ptr() as *mut f64);
    let up = SyncPtr(num.ut.as_ptr() as *mut f64);

    let fwd = |s: usize, par: bool| {
        let xq = xp;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let upds = sym.updates(s);
        let gather = |u: &crate::symbolic::Update, acc: &mut [f64]| {
            let d = u.d as usize;
            let kd = sym.k(d);
            let md = sym.m(d);
            let bd = sym.sn_begin[d];
            let pd = sym.pattern(d);
            let (r0, r1) = (u.r0 as usize, u.r1 as usize);
            let ld = unsafe { mat_ref(lp, sym.val_ptr[d], md, kd) };
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
        // SAFETY: x[b..b+k] is owned by supernode s during this call
        let xs = unsafe { std::slice::from_raw_parts_mut(xq.get().add(b), k) };
        let mut acc = vec![0.0f64; k];
        for pc in &partial {
            for (a, v) in acc.iter_mut().zip(pc) {
                *a += v;
            }
        }
        for c in 0..k {
            xs[c] += acc[c];
        }
        let ls = unsafe { mat_ref(lp, sym.val_ptr[s], sym.m(s), k) };
        let xm = MatMut::from_column_major_slice_mut(xs, k, 1);
        solve_unit_lower_triangular_in_place(ls.subrows(0, k), xm, Par::Seq);
    };

    let bwd = |s: usize, par: bool| {
        let xq = xp;
        let b = sym.sn_begin[s];
        let k = sym.k(s);
        let m = sym.m(s);
        let pat = sym.pattern(s);
        let ls = unsafe { mat_ref(lp, sym.val_ptr[s], m, k) };
        let mut v = vec![0.0f64; k];
        if m > k {
            let ut = unsafe { mat_ref(up, lsym.ut_ptr[s], m - k, k) };
            let dot = |c: usize| -> f64 {
                let col = ut.col(c);
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
        solve_upper_triangular_in_place(ls.subrows(0, k), xm, Par::Seq);
    };

    subtrees.par_iter().with_max_len(1).for_each(|&r| {
        for s in sym.first_desc[r]..=r {
            fwd(s, false);
        }
    });
    for &s in top {
        fwd(s, true);
    }
    for &s in top.iter().rev() {
        bwd(s, true);
    }
    subtrees.par_iter().with_max_len(1).for_each(|&r| {
        for s in (sym.first_desc[r]..=r).rev() {
            bwd(s, false);
        }
    });
}

/// y = A x for CSR
/// y = A x for a CSR matrix, rows in parallel. Every row is summed in the same order as
/// in [`csr_spmv`], so the result is bitwise identical to it for any thread count.
pub fn csr_spmv_par(n: usize, rp: &[usize], ci: &[usize], v: &[f64], x: &[f64], y: &mut [f64]) {
    use rayon::prelude::*;
    const ROWS: usize = 2048;
    y[..n].par_chunks_mut(ROWS).enumerate().for_each(|(c, yc)| {
        let r0 = c * ROWS;
        for (k, yr) in yc.iter_mut().enumerate() {
            let r = r0 + k;
            let mut acc = 0.0;
            for t in rp[r]..rp[r + 1] {
                acc += v[t] * x[ci[t]];
            }
            *yr = acc;
        }
    });
}

pub fn csr_spmv(n: usize, rp: &[usize], ci: &[usize], v: &[f64], x: &[f64], y: &mut [f64]) {
    for r in 0..n {
        let mut acc = 0.0;
        for t in rp[r]..rp[r + 1] {
            acc += v[t] * x[ci[t]];
        }
        y[r] = acc;
    }
}
