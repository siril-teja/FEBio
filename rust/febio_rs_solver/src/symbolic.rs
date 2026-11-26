//! Symbolic analysis for the supernodal LDL^T factorization.
//!
//! Pipeline (done once per sparsity pattern):
//!  1. build the full adjacency graph of the symmetric matrix,
//!  2. compress it into supervariables (columns with identical structure, e.g. the
//!     x/y/z dofs of one FE node),
//!  3. fill-reducing ordering of the compressed graph (METIS nested dissection, or AMD),
//!  4. elimination tree + postorder, column counts, relaxed supernodes (faer symbolic),
//!  5. precompute everything the numeric phase needs so that it is pure arithmetic:
//!     A -> panel scatter maps, supernodal tree, and the per-supernode list of
//!     descendant supernodes that update it (left-looking).

use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::sparse::linalg::cholesky::{simplicial, supernodal};
use faer::sparse::linalg::SymbolicSupernodalParams;
use faer::sparse::SymbolicSparseColMatRef;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ordering {
    Metis,
    Amd,
    Natural,
}

/// One descendant update: supernode `d` has rows `pd[r0..r1]` inside the columns of the
/// target supernode; rows `pd[r0..]` all lie in the target's index set.
#[derive(Clone, Copy, Debug)]
pub struct Update {
    pub d: u32,
    pub r0: u32,
    pub r1: u32,
}

pub struct Symbolic {
    pub n: usize,
    /// new -> old
    pub perm: Vec<usize>,
    /// old -> new
    pub iperm: Vec<usize>,
    pub ns: usize,
    /// supernode column ranges [sn_begin[s], sn_begin[s+1])
    pub sn_begin: Vec<usize>,
    /// off-diagonal row pattern of supernode s: rows[rows_ptr[s]..rows_ptr[s+1]] (permuted indices, sorted)
    pub rows_ptr: Vec<usize>,
    pub rows: Vec<usize>,
    /// offset of supernode s's panel (m_s x k_s, column-major) in the value array
    pub val_ptr: Vec<usize>,
    pub sn_parent: Vec<usize>, // usize::MAX = root
    pub sn_children_ptr: Vec<usize>,
    pub sn_children: Vec<usize>,
    pub roots: Vec<usize>,
    /// A entry scatter: for supernode s, entries amap[amap_ptr[s]..amap_ptr[s+1]] = (index into A values, position in panel)
    pub amap_ptr: Vec<usize>,
    pub amap: Vec<(usize, usize)>,
    pub upd_ptr: Vec<usize>,
    pub upd: Vec<Update>,
    /// accumulated flop estimate per subtree (used for scheduling decisions)
    pub subtree_flops: Vec<f64>,
    /// first (smallest) supernode index in the subtree of s (postorder => contiguous range)
    pub first_desc: Vec<usize>,
    pub flops: f64,
    pub nnz_l: usize,
    pub n_supervars: usize,
}

impl Symbolic {
    #[inline]
    pub fn k(&self, s: usize) -> usize {
        self.sn_begin[s + 1] - self.sn_begin[s]
    }
    #[inline]
    pub fn m(&self, s: usize) -> usize {
        self.k(s) + self.rows_ptr[s + 1] - self.rows_ptr[s]
    }
    #[inline]
    pub fn pattern(&self, s: usize) -> &[usize] {
        &self.rows[self.rows_ptr[s]..self.rows_ptr[s + 1]]
    }
    pub fn children(&self, s: usize) -> &[usize] {
        &self.sn_children[self.sn_children_ptr[s]..self.sn_children_ptr[s + 1]]
    }
    pub fn updates(&self, s: usize) -> &[Update] {
        &self.upd[self.upd_ptr[s]..self.upd_ptr[s + 1]]
    }
    pub fn len_val(&self) -> usize {
        self.val_ptr[self.ns]
    }
}

/// Lower-triangular CSC input (0-based), may contain the diagonal.
pub struct LowerCsc<'a> {
    pub n: usize,
    pub colptr: &'a [usize],
    pub rowidx: &'a [usize],
}

/// full symmetric adjacency without diagonal (CSR == CSC), sorted
fn full_adjacency(a: &LowerCsc) -> (Vec<usize>, Vec<usize>) {
    let n = a.n;
    let mut deg = vec![0usize; n + 1];
    for j in 0..n {
        for &i in &a.rowidx[a.colptr[j]..a.colptr[j + 1]] {
            if i != j {
                deg[i] += 1;
                deg[j] += 1;
            }
        }
    }
    let mut xadj = vec![0usize; n + 1];
    for i in 0..n {
        xadj[i + 1] = xadj[i] + deg[i];
    }
    let mut pos = xadj.clone();
    let mut adj = vec![0usize; xadj[n]];
    for j in 0..n {
        for &i in &a.rowidx[a.colptr[j]..a.colptr[j + 1]] {
            if i != j {
                adj[pos[i]] = j;
                pos[i] += 1;
                adj[pos[j]] = i;
                pos[j] += 1;
            }
        }
    }
    for i in 0..n {
        adj[xadj[i]..xadj[i + 1]].sort_unstable();
    }
    (xadj, adj)
}

/// consecutive columns j, j+1 with pattern(j) U {j} == pattern(j+1) U {j+1} form one supervariable
fn supervariables(n: usize, xadj: &[usize], adj: &[usize]) -> Vec<usize> {
    let mut start = vec![0usize];
    let closed_eq = |a: usize, b: usize| -> bool {
        let pa = &adj[xadj[a]..xadj[a + 1]];
        let pb = &adj[xadj[b]..xadj[b + 1]];
        if pa.len() != pb.len() {
            return false;
        }
        // compare merged closed neighborhoods
        let mut ia = 0;
        let mut ib = 0;
        let mut sa = false;
        let mut sb = false;
        loop {
            let va = if !sa && (ia >= pa.len() || a < pa[ia]) {
                Some(a)
            } else if ia < pa.len() {
                Some(pa[ia])
            } else {
                None
            };
            let vb = if !sb && (ib >= pb.len() || b < pb[ib]) {
                Some(b)
            } else if ib < pb.len() {
                Some(pb[ib])
            } else {
                None
            };
            match (va, vb) {
                (None, None) => return true,
                (Some(x), Some(y)) if x == y => {
                    if !sa && x == a {
                        sa = true
                    } else {
                        ia += 1
                    }
                    if !sb && y == b {
                        sb = true
                    } else {
                        ib += 1
                    }
                }
                _ => return false,
            }
        }
    };
    for j in 1..n {
        if !closed_eq(j - 1, j) {
            start.push(j);
        }
    }
    start.push(n);
    start
}

fn order_metis(nv: usize, xadj: &[usize], adj: &[usize], vwgt: &[i32]) -> Vec<usize> {
    use metis_sys as m;
    if nv == 0 {
        return vec![];
    }
    if adj.is_empty() || nv < 8 {
        return (0..nv).collect();
    }
    let mut nvtxs: m::idx_t = nv as m::idx_t;
    let mut x: Vec<m::idx_t> = xadj.iter().map(|&v| v as m::idx_t).collect();
    let mut a: Vec<m::idx_t> = adj.iter().map(|&v| v as m::idx_t).collect();
    let mut w: Vec<m::idx_t> = vwgt.iter().map(|&v| v as m::idx_t).collect();
    let mut perm = vec![0 as m::idx_t; nv];
    let mut iperm = vec![0 as m::idx_t; nv];
    let mut opts = [0 as m::idx_t; m::METIS_NOPTIONS as usize];
    unsafe {
        m::METIS_SetDefaultOptions(opts.as_mut_ptr());
    }
    opts[m::moptions_et_METIS_OPTION_NUMBERING as usize] = 0;
    let ret = unsafe {
        m::METIS_NodeND(
            &mut nvtxs,
            x.as_mut_ptr(),
            a.as_mut_ptr(),
            w.as_mut_ptr(),
            opts.as_mut_ptr(),
            perm.as_mut_ptr(),
            iperm.as_mut_ptr(),
        )
    };
    if ret != m::rstatus_et_METIS_OK {
        return (0..nv).collect();
    }
    // METIS: row i of A' is row perm[i] of A  (perm: new -> old)
    perm.iter().map(|&p| p as usize).collect()
}

fn order_amd(nv: usize, xadj: &[usize], adj: &[usize]) -> Vec<usize> {
    use faer::sparse::linalg::amd;
    if nv == 0 {
        return vec![];
    }
    let mut fwd = vec![0usize; nv];
    let mut inv = vec![0usize; nv];
    let a = SymbolicSparseColMatRef::new_checked(nv, nv, xadj, None, adj);
    let req = amd::order_scratch::<usize>(nv, adj.len());
    let mut mem = MemBuffer::new(req);
    match amd::order(
        &mut fwd,
        &mut inv,
        a,
        Default::default(),
        MemStack::new(&mut mem),
    ) {
        Ok(_) => fwd, // perm_fwd: new -> old? (faer: fwd[i] = old index at new position i)
        Err(_) => (0..nv).collect(),
    }
}

/// upper-triangular CSC (incl diagonal) of P A P^T given iperm (old->new)
fn permuted_upper(a: &LowerCsc, iperm: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let n = a.n;
    let mut cnt = vec![0usize; n + 1];
    for j in 0..n {
        for &i in &a.rowidx[a.colptr[j]..a.colptr[j + 1]] {
            let (pi, pj) = (iperm[i], iperm[j]);
            cnt[pi.max(pj)] += 1;
        }
    }
    let mut cp = vec![0usize; n + 1];
    for j in 0..n {
        cp[j + 1] = cp[j] + cnt[j];
    }
    let mut pos = cp.clone();
    let mut ri = vec![0usize; cp[n]];
    for j in 0..n {
        for &i in &a.rowidx[a.colptr[j]..a.colptr[j + 1]] {
            let (pi, pj) = (iperm[i], iperm[j]);
            let c = pi.max(pj);
            ri[pos[c]] = pi.min(pj);
            pos[c] += 1;
        }
    }
    for j in 0..n {
        ri[cp[j]..cp[j + 1]].sort_unstable();
    }
    (cp, ri)
}

fn etree_postorder(parent: &[isize]) -> Vec<usize> {
    let n = parent.len();
    let mut head = vec![usize::MAX; n];
    let mut next = vec![usize::MAX; n];
    // build child lists in reverse so that children are visited in increasing order
    for j in (0..n).rev() {
        if parent[j] >= 0 {
            let p = parent[j] as usize;
            next[j] = head[p];
            head[p] = j;
        }
    }
    let mut post = Vec::with_capacity(n);
    let mut stack = Vec::new();
    for r in 0..n {
        if parent[r] >= 0 {
            continue;
        }
        stack.push(r);
        while let Some(&top) = stack.last() {
            let c = head[top];
            if c == usize::MAX {
                post.push(top);
                stack.pop();
            } else {
                head[top] = next[c];
                stack.push(c);
            }
        }
    }
    post
}

pub fn analyze(a: &LowerCsc, ordering: Ordering) -> Result<Symbolic, String> {
    analyze_with_hint(a, ordering, None)
}

/// `hint`: a fill-reducing permutation (new -> old) to reuse instead of computing a new
/// ordering, e.g. the one of a previous, slightly different sparsity pattern (contact).
pub fn analyze_with_hint(
    a: &LowerCsc,
    ordering: Ordering,
    hint: Option<&[usize]>,
) -> Result<Symbolic, String> {
    let n = a.n;
    // ---- 1-3: graph, supervariables, ordering
    let prof = std::env::var_os("FRS_PROFILE").is_some();
    let mut tlast = std::time::Instant::now();
    let mut tick = |lab: &str| {
        if prof {
            eprintln!("  analyze {lab}: {:.3}s", tlast.elapsed().as_secs_f64());
        }
        tlast = std::time::Instant::now();
    };
    let (perm0, nsv) = match hint {
        Some(h) if h.len() == n => (h.to_vec(), 0),
        _ => {
            let (xadj, adj) = full_adjacency(a);
            let sv = supervariables(n, &xadj, &adj);
            let nsv = sv.len() - 1;
            let mut sv_of = vec![0usize; n];
            for s in 0..nsv {
                for j in sv[s]..sv[s + 1] {
                    sv_of[j] = s;
                }
            }
            tick("graph+supervars");
            let perm0: Vec<usize> = match ordering {
                Ordering::Natural => (0..n).collect(),
                _ => {
                    // compressed graph
                    let mut cx = vec![0usize; nsv + 1];
                    let mut ca = Vec::with_capacity(adj.len() / 4 + 1);
                    let mut mark = vec![usize::MAX; nsv];
                    for s in 0..nsv {
                        let j = sv[s];
                        mark[s] = s;
                        for &i in &adj[xadj[j]..xadj[j + 1]] {
                            let t = sv_of[i];
                            if mark[t] != s {
                                mark[t] = s;
                                ca.push(t);
                            }
                        }
                        cx[s + 1] = ca.len();
                    }
                    let vw: Vec<i32> = (0..nsv).map(|s| (sv[s + 1] - sv[s]) as i32).collect();
                    let sperm = match ordering {
                        Ordering::Metis => order_metis(nsv, &cx, &ca, &vw),
                        _ => order_amd(nsv, &cx, &ca),
                    };
                    let mut p = Vec::with_capacity(n);
                    for &s in &sperm {
                        p.extend(sv[s]..sv[s + 1]);
                    }
                    p
                }
            };
            drop(adj);
            drop(xadj);
            (perm0, nsv)
        }
    };
    let mut iperm0 = vec![0usize; n];
    for (new, &old) in perm0.iter().enumerate() {
        iperm0[old] = new;
    }

    // ---- 4: etree postorder, then faer supernodal symbolic
    tick("ordering");
    let (cp, ri) = permuted_upper(a, &iperm0);
    let mut etree = vec![0isize; n];
    let mut cc = vec![0usize; n];
    {
        let req = simplicial::prefactorize_symbolic_cholesky_scratch::<usize>(n, ri.len());
        let mut mem = MemBuffer::new(req);
        let au = SymbolicSparseColMatRef::new_checked(n, n, &cp, None, &ri);
        simplicial::prefactorize_symbolic_cholesky::<usize>(
            &mut etree,
            &mut cc,
            au,
            MemStack::new(&mut mem),
        );
    }
    tick("etree1");
    let post = etree_postorder(&etree);
    let perm: Vec<usize> = post.iter().map(|&p| perm0[p]).collect();
    let mut iperm = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        iperm[old] = new;
    }
    // relabel etree and column counts (invariant under a postorder relabeling)
    {
        let mut inv_post = vec![0usize; n];
        for (newi, &p) in post.iter().enumerate() {
            inv_post[p] = newi;
        }
        let mut et2 = vec![-1isize; n];
        let mut cc2 = vec![0usize; n];
        for j in 0..n {
            let nj = inv_post[j];
            et2[nj] = if etree[j] >= 0 {
                inv_post[etree[j] as usize] as isize
            } else {
                -1
            };
            cc2[nj] = cc[j];
        }
        etree = et2;
        cc = cc2;
    }
    let (cp, ri) = permuted_upper(a, &iperm);
    tick("etree2");
    let sym = {
        let req = StackReq::any_of(&[
            supernodal::factorize_supernodal_symbolic_cholesky_scratch::<usize>(n),
        ]);
        let mut mem = MemBuffer::new(req);
        let au = SymbolicSparseColMatRef::new_checked(n, n, &cp, None, &ri);
        let et = unsafe { simplicial::EliminationTreeRef::from_inner(&etree) };
        supernodal::factorize_supernodal_symbolic_cholesky::<usize>(
            au,
            et,
            &cc,
            MemStack::new(&mut mem),
            SymbolicSupernodalParams::default(),
        )
        .map_err(|e| format!("supernodal symbolic failed: {e:?}"))?
    };
    drop(cp);
    drop(ri);

    tick("faer supernodal");
    // ---- 5: supernodal tree and numeric-phase maps
    let ns = sym.n_supernodes();
    let mut sn_begin = Vec::with_capacity(ns + 1);
    sn_begin.extend_from_slice(sym.supernode_begin());
    sn_begin.push(n);
    let mut rows_ptr = vec![0usize; ns + 1];
    let mut rows = Vec::new();
    for s in 0..ns {
        let mut p: Vec<usize> = sym.supernode(s).pattern().to_vec();
        p.sort_unstable();
        rows.extend_from_slice(&p);
        rows_ptr[s + 1] = rows.len();
    }
    let mut sn_of = vec![0usize; n];
    for s in 0..ns {
        for j in sn_begin[s]..sn_begin[s + 1] {
            sn_of[j] = s;
        }
    }
    let mut val_ptr = vec![0usize; ns + 1];
    let mut flops = 0.0f64;
    let mut sflops = vec![0.0f64; ns];
    let mut nnz_l = 0usize;
    for s in 0..ns {
        let k = sn_begin[s + 1] - sn_begin[s];
        let m = k + rows_ptr[s + 1] - rows_ptr[s];
        val_ptr[s + 1] = val_ptr[s] + m * k;
        nnz_l += m * k - k * (k + 1) / 2;
        // LDL^T flops for this supernode's columns: sum over cols of (rows below)^2
        let mut f = 0.0;
        for c in 0..k {
            let below = (m - c - 1) as f64;
            f += below * below + below;
        }
        sflops[s] = f;
        flops += f;
    }
    let mut sn_parent = vec![usize::MAX; ns];
    for s in 0..ns {
        let p = &rows[rows_ptr[s]..rows_ptr[s + 1]];
        if let Some(&r) = p.first() {
            let t = sn_of[r];
            if t <= s {
                return Err("supernodal tree not postordered".into());
            }
            sn_parent[s] = t;
        }
    }
    let mut ccount = vec![0usize; ns + 1];
    let mut roots = Vec::new();
    for s in 0..ns {
        if sn_parent[s] == usize::MAX {
            roots.push(s);
        } else {
            ccount[sn_parent[s]] += 1;
        }
    }
    let mut sn_children_ptr = vec![0usize; ns + 1];
    for s in 0..ns {
        sn_children_ptr[s + 1] = sn_children_ptr[s] + ccount[s];
    }
    let mut cpos = sn_children_ptr.clone();
    let mut sn_children = vec![0usize; sn_children_ptr[ns]];
    for s in 0..ns {
        let p = sn_parent[s];
        if p != usize::MAX {
            sn_children[cpos[p]] = s;
            cpos[p] += 1;
        }
    }
    // subtree flops (postorder: children before parents)
    let mut subtree_flops = sflops.clone();
    for s in 0..ns {
        if sn_parent[s] != usize::MAX {
            let v = subtree_flops[s];
            subtree_flops[sn_parent[s]] += v;
        }
    }

    let mut first_desc: Vec<usize> = (0..ns).collect();
    for s in 0..ns {
        let p = sn_parent[s];
        if p != usize::MAX && first_desc[s] < first_desc[p] {
            first_desc[p] = first_desc[s];
        }
    }

    tick("tree");
    // A -> panel scatter
    let mut acount = vec![0usize; ns + 1];
    for j in 0..n {
        for &i in &a.rowidx[a.colptr[j]..a.colptr[j + 1]] {
            let lo = iperm[i].min(iperm[j]);
            acount[sn_of[lo]] += 1;
        }
    }
    let mut amap_ptr = vec![0usize; ns + 1];
    for s in 0..ns {
        amap_ptr[s + 1] = amap_ptr[s] + acount[s];
    }
    let mut apos = amap_ptr.clone();
    let mut amap = vec![(0usize, 0usize); amap_ptr[ns]];
    // per-supernode local row lookup: rows inside the supernode are direct, others by binary search
    for j in 0..n {
        for t in a.colptr[j]..a.colptr[j + 1] {
            let i = a.rowidx[t];
            let (pi, pj) = (iperm[i], iperm[j]);
            let (lo, hi) = (pi.min(pj), pi.max(pj));
            let s = sn_of[lo];
            let b = sn_begin[s];
            let e = sn_begin[s + 1];
            let k = e - b;
            let m = k + rows_ptr[s + 1] - rows_ptr[s];
            let lr = if hi < e {
                hi - b
            } else {
                let p = &rows[rows_ptr[s]..rows_ptr[s + 1]];
                match p.binary_search(&hi) {
                    Ok(x) => k + x,
                    Err(_) => return Err("A entry outside symbolic pattern".into()),
                }
            };
            let lc = lo - b;
            amap[apos[s]] = (t, lc * m + lr);
            apos[s] += 1;
        }
    }
    tick("amap");
    // updates (left-looking): for each d, group its pattern rows by target supernode
    let mut ucount = vec![0usize; ns + 1];
    for d in 0..ns {
        let p = &rows[rows_ptr[d]..rows_ptr[d + 1]];
        let mut r = 0;
        while r < p.len() {
            let s = sn_of[p[r]];
            let e = sn_begin[s + 1];
            while r < p.len() && p[r] < e {
                r += 1;
            }
            ucount[s] += 1;
        }
    }
    let mut upd_ptr = vec![0usize; ns + 1];
    for s in 0..ns {
        upd_ptr[s + 1] = upd_ptr[s] + ucount[s];
    }
    let mut upos = upd_ptr.clone();
    let mut upd = vec![Update { d: 0, r0: 0, r1: 0 }; upd_ptr[ns]];
    for d in 0..ns {
        let p = &rows[rows_ptr[d]..rows_ptr[d + 1]];
        let mut r = 0;
        while r < p.len() {
            let s = sn_of[p[r]];
            let e = sn_begin[s + 1];
            let r0 = r;
            while r < p.len() && p[r] < e {
                r += 1;
            }
            upd[upos[s]] = Update {
                d: d as u32,
                r0: r0 as u32,
                r1: r as u32,
            };
            upos[s] += 1;
        }
    }

    tick("updates");
    Ok(Symbolic {
        n,
        perm,
        iperm,
        ns,
        sn_begin,
        rows_ptr,
        rows,
        val_ptr,
        sn_parent,
        sn_children_ptr,
        sn_children,
        roots,
        amap_ptr,
        amap,
        upd_ptr,
        upd,
        subtree_flops,
        first_desc,
        flops,
        nnz_l,
        n_supervars: nsv,
    })
}

impl Symbolic {
    /// (critical-path flops, total flops) of the supernodal tree, using own-flops per supernode
    pub fn critical_path(&self) -> (f64, f64) {
        let mut own = vec![0.0f64; self.ns];
        for s in 0..self.ns {
            let k = self.k(s);
            let m = self.m(s);
            for c in 0..k {
                let below = (m - c - 1) as f64;
                own[s] += below * below + below;
            }
        }
        let mut cp = own.clone();
        for s in 0..self.ns {
            let p = self.sn_parent[s];
            if p != usize::MAX {
                let v = cp[s];
                // parent's cp = own + max child cp: computed incrementally
                if v + own[p] > cp[p] {
                    cp[p] = v + own[p];
                }
            }
        }
        let crit = self.roots.iter().map(|&r| cp[r]).fold(0.0, f64::max);
        (crit, own.iter().sum())
    }
}
