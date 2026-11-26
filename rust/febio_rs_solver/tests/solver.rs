//! Correctness tests on FE-like 3D stencil matrices (3 dofs per node, 27-point coupling).
use febio_rs_solver::{lu::csr_spmv, sym_spmv_lower, Solver};

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % 10_000) as f64 / 10_000.0
    }
}

/// full CSR pattern + values of a structurally symmetric matrix on a g^3 grid.
/// `asym` scales the unsymmetric perturbation, `diag` the diagonal.
fn stencil(g: usize, asym: f64, diag: f64, seed: u64) -> (usize, Vec<i32>, Vec<i32>, Vec<f64>) {
    let mut rng = Rng(seed);
    let n = 3 * g * g * g;
    let id = |i: usize, j: usize, k: usize| (k * g + j) * g + i;
    let mut rows: Vec<Vec<(usize, f64)>> = vec![vec![]; n];
    for k in 0..g {
        for j in 0..g {
            for i in 0..g {
                let me = id(i, j, k);
                for dk in -1i64..=1 {
                    for dj in -1i64..=1 {
                        for di in -1i64..=1 {
                            let (ii, jj, kk) = (i as i64 + di, j as i64 + dj, k as i64 + dk);
                            if ii < 0
                                || jj < 0
                                || kk < 0
                                || ii >= g as i64
                                || jj >= g as i64
                                || kk >= g as i64
                            {
                                continue;
                            }
                            let o = id(ii as usize, jj as usize, kk as usize);
                            for a in 0..3 {
                                for b in 0..3 {
                                    let (r, c) = (3 * me + a, 3 * o + b);
                                    if c < r {
                                        continue;
                                    }
                                    if r == c {
                                        rows[r].push((c, diag));
                                    } else {
                                        let v = -0.4 - 0.3 * rng.next();
                                        let d = asym * (rng.next() - 0.5);
                                        rows[r].push((c, v + d));
                                        rows[c].push((r, v - d));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    let mut rp = vec![0i32];
    let mut ci = vec![];
    let mut va = vec![];
    for r in rows.iter_mut() {
        r.sort_by_key(|e| e.0);
        for &(c, v) in r.iter() {
            ci.push(c as i32);
            va.push(v);
        }
        rp.push(ci.len() as i32);
    }
    (n, rp, ci, va)
}

fn lower_of(n: usize, rp: &[i32], ci: &[i32], v: &[f64]) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    // CSR upper == CSC lower for a symmetric matrix
    let mut cp = vec![0i32];
    let (mut ri, mut vv) = (vec![], vec![]);
    for r in 0..n {
        for t in rp[r] as usize..rp[r + 1] as usize {
            if ci[t] as usize >= r {
                ri.push(ci[t]);
                vv.push(v[t]);
            }
        }
        cp.push(ri.len() as i32);
    }
    (cp, ri, vv)
}

fn rel_err(x: &[f64], y: &[f64]) -> f64 {
    let e: f64 = x.iter().zip(y).map(|(a, b)| (a - b) * (a - b)).sum();
    let n: f64 = y.iter().map(|a| a * a).sum();
    (e / n).sqrt()
}

fn xtrue(n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5)
        .collect()
}

fn check_sym(g: usize, diag: f64, threads: usize, tol: f64) {
    let (n, rp, ci, v) = stencil(g, 0.0, diag, 7);
    let (cp, ri, lv) = lower_of(n, &rp, &ci, &v);
    let mut s = Solver::new(threads);
    s.analyze_symmetric(n, &cp, &ri, 0).unwrap();
    let xt = xtrue(n);
    let cpu: Vec<usize> = cp.iter().map(|&x| x as usize).collect();
    let riu: Vec<usize> = ri.iter().map(|&x| x as usize).collect();
    let mut b = vec![0.0; n];
    sym_spmv_lower(n, &cpu, &riu, &lv, &xt, &mut b);
    for _ in 0..2 {
        s.factor(&lv).unwrap();
        let mut x = vec![0.0; n];
        s.solve(&mut x, &b).unwrap();
        let e = rel_err(&x, &xt);
        assert!(
            e < tol,
            "symmetric g={g} diag={diag} T={threads}: rel error {e:e}"
        );
    }
}

fn check_unsym(g: usize, asym: f64, diag: f64, threads: usize, tol: f64) {
    let (n, rp, ci, v) = stencil(g, asym, diag, 11);
    let mut s = Solver::new(threads);
    s.analyze_unsymmetric(n, &rp, &ci, 0).unwrap();
    let xt = xtrue(n);
    let rpu: Vec<usize> = rp.iter().map(|&x| x as usize).collect();
    let ciu: Vec<usize> = ci.iter().map(|&x| x as usize).collect();
    let mut b = vec![0.0; n];
    csr_spmv(n, &rpu, &ciu, &v, &xt, &mut b);
    for _ in 0..2 {
        s.factor(&v).unwrap();
        let mut x = vec![0.0; n];
        s.solve(&mut x, &b).unwrap();
        let e = rel_err(&x, &xt);
        assert!(
            e < tol,
            "unsymmetric g={g} asym={asym} diag={diag} T={threads}: rel error {e:e}"
        );
    }
}

#[test]
fn symmetric_spd() {
    for t in [1, 4] {
        check_sym(9, 30.0, t, 1e-11);
    }
}

#[test]
fn symmetric_indefinite() {
    for t in [1, 4] {
        check_sym(9, -1.5, t, 1e-7);
    }
}

#[test]
fn unsymmetric_diag_dominant() {
    for t in [1, 4] {
        check_unsym(9, 0.6, 30.0, t, 1e-11);
    }
}

#[test]
fn unsymmetric_strongly_nonsymmetric() {
    for t in [1, 4] {
        check_unsym(9, 4.0, 8.0, t, 1e-9);
    }
}

#[test]
fn unsymmetric_indefinite() {
    for t in [1, 4] {
        check_unsym(9, 1.0, -1.5, t, 1e-7);
    }
}

#[test]
fn larger_tree_parallel() {
    check_sym(16, 30.0, 4, 1e-11);
    check_unsym(16, 1.0, 30.0, 4, 1e-11);
}
