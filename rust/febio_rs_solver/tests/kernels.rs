//! The solver relies on faer's multithreaded dense kernels being bitwise identical to their
//! sequential versions (parallelism over output blocks only, no split reductions).
//! These tests guard that property, e.g. across faer upgrades.
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::linalg::cholesky::ldlt::factor::{
    cholesky_in_place, cholesky_in_place_scratch, LdltParams, LdltRegularization,
};
use faer::linalg::matmul::matmul;
use faer::linalg::triangular_solve::solve_unit_lower_triangular_in_place;
use faer::{Accum, Mat, Par, Spec};

fn same(a: &Mat<f64>, b: &Mat<f64>) -> bool {
    (0..a.ncols()).all(|j| (0..a.nrows()).all(|i| a[(i, j)].to_bits() == b[(i, j)].to_bits()))
}

#[test]
fn gemm_parallel_matches_sequential() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for &(m, n, k) in &[
        (1200usize, 300usize, 200usize),
        (300, 300, 1500),
        (2000, 48, 1500),
    ] {
        let a = Mat::<f64>::from_fn(m, k, |i, j| ((i * 7 + j * 13) % 101) as f64 / 101.0 - 0.5);
        let b = Mat::<f64>::from_fn(k, n, |i, j| ((i * 3 + j * 17) % 97) as f64 / 97.0 - 0.5);
        let mut c1 = Mat::<f64>::zeros(m, n);
        let mut c2 = Mat::<f64>::zeros(m, n);
        matmul(
            c1.as_mut(),
            Accum::Replace,
            a.as_ref(),
            b.as_ref(),
            1.0,
            Par::Seq,
        );
        pool.install(|| {
            matmul(
                c2.as_mut(),
                Accum::Replace,
                a.as_ref(),
                b.as_ref(),
                1.0,
                Par::rayon(4),
            )
        });
        assert!(
            same(&c1, &c2),
            "gemm {m}x{n}x{k} differs between Seq and rayon"
        );
    }
}

#[test]
fn ldlt_and_trsm_parallel_match_sequential() {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let n = 1500;
    let a = Mat::<f64>::from_fn(n, n, |i, j| {
        if i == j {
            n as f64
        } else {
            (((i * j) % 37) as f64) / 37.0 - 0.5 + (((i + j) % 11) as f64) / 50.0
        }
    });
    let a = Mat::<f64>::from_fn(n, n, |i, j| 0.5 * (a[(i, j)] + a[(j, i)]));
    let fac = |par: Par| {
        let mut l = a.clone();
        let params: Spec<LdltParams, f64> = Default::default();
        let mut mem = MemBuffer::new(cholesky_in_place_scratch::<f64>(n, par, params));
        let reg = LdltRegularization {
            dynamic_regularization_signs: None,
            dynamic_regularization_delta: 1e-12,
            dynamic_regularization_epsilon: 1e-12,
        };
        cholesky_in_place(l.as_mut(), reg, par, MemStack::new(&mut mem), params).unwrap();
        l
    };
    let l1 = fac(Par::Seq);
    let l2 = pool.install(|| fac(Par::rayon(4)));
    assert!(same(&l1, &l2), "LDLT differs between Seq and rayon");
    let b = Mat::<f64>::from_fn(n, 500, |i, j| ((i * 5 + j) % 23) as f64 - 11.0);
    let (mut x1, mut x2) = (b.clone(), b.clone());
    solve_unit_lower_triangular_in_place(l1.as_ref(), x1.as_mut(), Par::Seq);
    pool.install(|| solve_unit_lower_triangular_in_place(l1.as_ref(), x2.as_mut(), Par::rayon(4)));
    assert!(
        same(&x1, &x2),
        "triangular solve differs between Seq and rayon"
    );
}
