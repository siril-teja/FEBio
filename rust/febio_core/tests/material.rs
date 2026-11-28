//! Material tests (ROADMAP 4.4): stress-free reference state, small-strain limit, and the
//! tangent against a finite-difference Truesdell rate of the Kirchhoff stress.
use febio_core::element::{det3, Mat3};
use febio_core::material::{ElasticMaterial, NeoHookean, OrthoElastic};

fn ortho() -> OrthoElastic {
    // Finley et al. cortical bone
    OrthoElastic::new(
        8000.0, 8000.0, 12000.0, 2857.0, 3677.0, 4444.0, 0.4, 0.3, 0.35,
    )
    .unwrap()
}

fn mats() -> Vec<(&'static str, Box<dyn ElasticMaterial>)> {
    vec![
        ("neo-Hookean", Box::new(NeoHookean { e: 3.4, v: 0.45 })),
        ("orthotropic", Box::new(ortho())),
    ]
}

const I3: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn mul(a: &Mat3, b: &Mat3) -> Mat3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

#[test]
fn stress_free_reference() {
    for (n, m) in mats() {
        let s = m.stress(&I3);
        assert!(s.iter().flatten().all(|v| v.abs() < 1e-12), "{n}");
    }
}

#[test]
fn small_strain_limit() {
    let eps = 1e-7;
    let h: Mat3 = [[0.3, 0.1, -0.2], [0.1, -0.5, 0.25], [-0.2, 0.25, 0.4]]; // symmetric
    let mut f = I3;
    for i in 0..3 {
        for j in 0..3 {
            f[i][j] += eps * h[i][j];
        }
    }
    // neo-Hookean: sigma = lambda tr(e) I + 2 mu e
    let (e, v) = (3.4, 0.45);
    let lam = v * e / ((1.0 + v) * (1.0 - 2.0 * v));
    let mu = 0.5 * e / (1.0 + v);
    let s = NeoHookean { e, v }.stress(&f);
    let tr = h[0][0] + h[1][1] + h[2][2];
    for i in 0..3 {
        for j in 0..3 {
            let lin = eps * (lam * tr * I3[i][j] + 2.0 * mu * h[i][j]);
            assert!(
                (s[i][j] - lin).abs() < 1e-6 * eps * mu,
                "neo-Hookean {i}{j}"
            );
        }
    }
    // orthotropic in the global axes: e = C^-1 sigma with the engineering constants
    let s = ortho().stress(&f);
    let (e1, e2, e3, g12, g23, g31, v12, v23, v31) = (
        8000.0, 8000.0, 12000.0, 2857.0, 3677.0, 4444.0, 0.4, 0.3, 0.35,
    );
    let ex = s[0][0] / e1 - v12 / e1 * s[1][1] - v31 / e3 * s[2][2];
    let ey = -v12 / e1 * s[0][0] + s[1][1] / e2 - v23 / e2 * s[2][2];
    let ez = -v31 / e3 * s[0][0] - v23 / e2 * s[1][1] + s[2][2] / e3;
    let tol = 1e-6 * eps;
    assert!(
        (ex - eps * h[0][0]).abs() < tol
            && (ey - eps * h[1][1]).abs() < tol
            && (ez - eps * h[2][2]).abs() < tol
    );
    assert!((s[0][1] - 2.0 * g12 * eps * h[0][1]).abs() < 1e-6 * eps * g12);
    assert!((s[1][2] - 2.0 * g23 * eps * h[1][2]).abs() < 1e-6 * eps * g23);
    assert!((s[0][2] - 2.0 * g31 * eps * h[0][2]).abs() < 1e-6 * eps * g31);
}

/// J c : d = d(tau)/dt - l tau - tau l^T with l = d (symmetric), tau = J sigma.
#[test]
fn tangent_matches_truesdell_rate() {
    let f0: Mat3 = [[1.08, 0.12, -0.05], [0.03, 0.93, 0.07], [-0.06, 0.04, 1.11]];
    let d: Mat3 = [[0.2, 0.05, -0.1], [0.05, -0.3, 0.15], [-0.1, 0.15, 0.25]];
    let h = 1e-6;
    for (n, m) in mats() {
        let tau = |eps: f64| {
            let mut g = I3;
            for i in 0..3 {
                for j in 0..3 {
                    g[i][j] += eps * d[i][j];
                }
            }
            let f = mul(&g, &f0);
            let s = m.stress(&f);
            let jj = det3(&f);
            s.map(|r| r.map(|x| x * jj))
        };
        let (tp, tm, t0) = (tau(h), tau(-h), tau(0.0));
        let c = m.tangent(&f0);
        let j0 = det3(&f0);
        let mut worst: f64 = 0.0;
        let scale = t0
            .iter()
            .flatten()
            .fold(0.0f64, |a, b| a.max(b.abs()))
            .max(1.0);
        for i in 0..3 {
            for k in 0..3 {
                let rate = (tp[i][k] - tm[i][k]) / (2.0 * h);
                let mut lt = 0.0;
                for r in 0..3 {
                    lt += d[i][r] * t0[r][k] + t0[i][r] * d[k][r];
                }
                let mut cd = 0.0;
                for p in 0..3 {
                    for q in 0..3 {
                        cd += c[i][k][p][q] * d[p][q];
                    }
                }
                worst = worst.max((rate - lt - j0 * cd).abs() / scale);
            }
        }
        assert!(worst < 1e-6, "{n}: tangent mismatch {worst:e}");
    }
}
