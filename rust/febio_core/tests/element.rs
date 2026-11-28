//! Element tests (ROADMAP 4.3): shape-function identities, integration rules, volumes and
//! the patch test (an affine deformation gives the exact F at every integration point).
use febio_core::element::{det3, shape, shape_deriv, ElementTraits};
use febio_core::ElemType;

fn kinds() -> Vec<ElemType> {
    vec![ElemType::Hex8, ElemType::Tet4, ElemType::Tet10]
}

/// parametric coordinates of the nodes, in FEBio order
fn node_coords(k: &ElemType) -> Vec<[f64; 3]> {
    match k {
        ElemType::Hex8 => vec![
            [-1., -1., -1.],
            [1., -1., -1.],
            [1., 1., -1.],
            [-1., 1., -1.],
            [-1., -1., 1.],
            [1., -1., 1.],
            [1., 1., 1.],
            [-1., 1., 1.],
        ],
        ElemType::Tet4 => vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        ElemType::Tet10 => {
            let v = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
            let mid = |a: usize, b: usize| {
                [
                    (v[a][0] + v[b][0]) / 2.,
                    (v[a][1] + v[b][1]) / 2.,
                    (v[a][2] + v[b][2]) / 2.,
                ]
            };
            vec![
                v[0],
                v[1],
                v[2],
                v[3],
                mid(0, 1),
                mid(1, 2),
                mid(2, 0),
                mid(0, 3),
                mid(1, 3),
                mid(2, 3),
            ]
        }
        _ => unreachable!(),
    }
}

#[test]
fn kronecker_delta_and_partition_of_unity() {
    for k in kinds() {
        let xi = node_coords(&k);
        for (a, p) in xi.iter().enumerate() {
            let h = shape(&k, p[0], p[1], p[2]);
            for (b, v) in h.iter().enumerate() {
                let e = if a == b { 1.0 } else { 0.0 };
                assert!((v - e).abs() < 1e-14, "{k:?} H{b} at node {a} = {v}");
            }
        }
        for p in [[0.1, 0.2, 0.3], [-0.4, 0.15, 0.05]] {
            let s: f64 = shape(&k, p[0], p[1], p[2]).iter().sum();
            assert!((s - 1.0).abs() < 1e-14);
            let d = shape_deriv(&k, p[0], p[1], p[2]);
            for c in 0..3 {
                assert!(d.iter().map(|x| x[c]).sum::<f64>().abs() < 1e-13);
            }
        }
    }
}

#[test]
fn derivatives_match_finite_differences() {
    let h = 1e-6;
    for k in kinds() {
        let p = [0.13, 0.21, 0.17];
        let d = shape_deriv(&k, p[0], p[1], p[2]);
        for c in 0..3 {
            let mut pp = p;
            let mut pm = p;
            pp[c] += h;
            pm[c] -= h;
            let hp = shape(&k, pp[0], pp[1], pp[2]);
            let hm = shape(&k, pm[0], pm[1], pm[2]);
            for a in 0..hp.len() {
                let fd = (hp[a] - hm[a]) / (2.0 * h);
                assert!((fd - d[a][c]).abs() < 1e-8, "{k:?} dH{a}/dxi{c}");
            }
        }
    }
}

#[test]
fn reference_volumes() {
    for (k, v) in [
        (ElemType::Hex8, 8.0),
        (ElemType::Tet4, 1.0 / 6.0),
        (ElemType::Tet10, 1.0 / 6.0),
    ] {
        let tr = ElementTraits::new(&k);
        let x0 = node_coords(&k);
        assert!((tr.volume(&x0) - v).abs() < 1e-12, "{k:?}");
        assert!(
            (tr.weights.iter().sum::<f64>() - v).abs() < 1e-8,
            "{k:?} weights"
        );
    }
}

/// a distorted reference configuration (curved edges for tet10)
fn distorted(k: &ElemType) -> Vec<[f64; 3]> {
    node_coords(k)
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let w = 0.07 * ((i * 37 % 11) as f64 - 5.0) / 5.0;
            [
                2.0 * p[0] + 0.3 * p[1] + w,
                1.5 * p[1] - 0.2 * p[2] - 0.5 * w,
                1.2 * p[2] + 0.1 * p[0] + 0.3 * w,
            ]
        })
        .collect()
}

#[test]
fn patch_test_affine_deformation_is_exact() {
    let a = [[1.1, 0.2, -0.05], [0.03, 0.9, 0.1], [-0.07, 0.04, 1.2]];
    let c = [0.5, -0.3, 2.0];
    for k in kinds() {
        let tr = ElementTraits::new(&k);
        let x0 = distorted(&k);
        let x: Vec<[f64; 3]> = x0
            .iter()
            .map(|p| {
                let mut q = c;
                for i in 0..3 {
                    for j in 0..3 {
                        q[i] += a[i][j] * p[j];
                    }
                }
                q
            })
            .collect();
        for g in 0..tr.gauss_points() {
            let f = tr.deformation_gradient(&x0, &x, g);
            for i in 0..3 {
                for j in 0..3 {
                    assert!((f[i][j] - a[i][j]).abs() < 1e-12, "{k:?} g{g} F[{i}][{j}]");
                }
            }
            assert!((det3(&f) - det3(&a)).abs() < 1e-12);
        }
        // the current volume is det(A) times the reference volume
        let v0 = tr.volume(&x0);
        assert!(
            (tr.volume(&x) - det3(&a) * v0).abs() < 1e-11 * v0.abs().max(1.0),
            "{k:?} volume"
        );
    }
}
