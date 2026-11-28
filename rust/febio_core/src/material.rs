//! Elastic materials (ROADMAP 4.4): neo-Hookean and orthotropic elastic, with FEBio's
//! formulations (FEBioMech/FENeoHookean.cpp, FEOrthoElastic.cpp). Each material returns
//! the Cauchy stress and the spatial elasticity tensor for a deformation gradient F.
use crate::element::{det3, inv3, Mat3};
use crate::model::Material;

/// Fourth-order tensor, full storage: c[i][j][k][l].
pub type Tens4 = [[[[f64; 3]; 3]; 3]; 3];

pub trait ElasticMaterial: Send + Sync {
    /// Cauchy stress at deformation gradient F
    fn stress(&self, f: &Mat3) -> Mat3;
    /// spatial elasticity tensor (FEBio's `Tangent`) at deformation gradient F
    fn tangent(&self, f: &Mat3) -> Tens4;
}

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

fn transpose(a: &Mat3) -> Mat3 {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}

const I3: Mat3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// left Cauchy-Green tensor b = F F^T
pub fn left_cauchy_green(f: &Mat3) -> Mat3 {
    mul(f, &transpose(f))
}

/// right Cauchy-Green tensor C = F^T F
pub fn right_cauchy_green(f: &Mat3) -> Mat3 {
    mul(&transpose(f), f)
}

/// a ⊗ b + b ⊗ a  (FEBio dyad1s(a, b))
fn dyad1s(a: &Mat3, b: &Mat3, s: f64, c: &mut Tens4) {
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    c[i][j][k][l] += s * (a[i][j] * b[k][l] + b[i][j] * a[k][l]);
                }
            }
        }
    }
}

/// (a_ik b_jl + a_il b_jk + b_ik a_jl + b_il a_jk) / 2  (FEBio dyad4s(a, b))
fn dyad4s(a: &Mat3, b: &Mat3, s: f64, c: &mut Tens4) {
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    c[i][j][k][l] += s
                        * 0.5
                        * (a[i][k] * b[j][l]
                            + a[i][l] * b[j][k]
                            + b[i][k] * a[j][l]
                            + b[i][l] * a[j][k]);
                }
            }
        }
    }
}

/// Compressible neo-Hookean (FEBio "neo-Hookean"):
/// sigma = mu/J (b - I) + lambda ln(J)/J I.
#[derive(Clone, Debug)]
pub struct NeoHookean {
    pub e: f64,
    pub v: f64,
}

impl NeoHookean {
    fn lame(&self) -> (f64, f64) {
        let lam = self.v * self.e / ((1.0 + self.v) * (1.0 - 2.0 * self.v));
        let mu = 0.5 * self.e / (1.0 + self.v);
        (lam, mu)
    }
}

impl ElasticMaterial for NeoHookean {
    fn stress(&self, f: &Mat3) -> Mat3 {
        let j = det3(f);
        let (lam, mu) = self.lame();
        let b = left_cauchy_green(f);
        let mut s = [[0.0; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                s[i][k] = (b[i][k] - I3[i][k]) * (mu / j) + I3[i][k] * (lam * j.ln() / j);
            }
        }
        s
    }

    fn tangent(&self, f: &Mat3) -> Tens4 {
        let j = det3(f);
        let (lam, mu) = self.lame();
        let lam1 = lam / j;
        let mu1 = (mu - lam * j.ln()) / j;
        let mut c = [[[[0.0; 3]; 3]; 3]; 3];
        // dyad1s(I) = I ⊗ I (single argument: no factor 2); dyad4s(I) = (I_ik I_jl + I_il I_jk)/2
        for i in 0..3 {
            for j2 in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        c[i][j2][k][l] = lam1 * I3[i][j2] * I3[k][l]
                            + 2.0 * mu1 * 0.5 * (I3[i][k] * I3[j2][l] + I3[i][l] * I3[j2][k]);
                    }
                }
            }
        }
        c
    }
}

/// Orthotropic elastic (FEBio "orthotropic elastic", a St. Venant-Kirchhoff-type
/// hyperelastic law in the material axes). `axes` holds the material axes a0_i as columns
/// (identity: global axes, FEBio's default when no mat_axis is given).
#[derive(Clone, Debug)]
pub struct OrthoElastic {
    pub lam: [[f64; 3]; 3],
    pub mu: [f64; 3],
    pub axes: Mat3,
}

impl OrthoElastic {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        e1: f64,
        e2: f64,
        e3: f64,
        g12: f64,
        g23: f64,
        g31: f64,
        v12: f64,
        v23: f64,
        v31: f64,
    ) -> Result<OrthoElastic, String> {
        // FEOrthoElastic::EvaluateLameCoefficients
        let mu = [g12 + g31 - g23, g12 - g31 + g23, -g12 + g31 + g23];
        let c: Mat3 = [
            [1.0 / e1, -v12 / e1, -v31 / e3],
            [-v12 / e1, 1.0 / e2, -v23 / e2],
            [-v31 / e3, -v23 / e2, 1.0 / e3],
        ];
        // positive definite (leading principal minors, as FEBio's eigenvalue check)
        let d = det3(&c);
        if !(c[0][0] > 0.0 && c[0][0] * c[1][1] - c[0][1] * c[1][0] > 0.0 && d > 0.0) {
            return Err("orthotropic elastic: compliance matrix is not positive definite".into());
        }
        let k = inv3(&c, d);
        let mut lam = k;
        for i in 0..3 {
            lam[i][i] = k[i][i] - 2.0 * mu[i];
        }
        Ok(OrthoElastic { lam, mu, axes: I3 })
    }

    /// K_i, A_i (current texture tensors) for deformation gradient F
    fn texture(&self, f: &Mat3) -> ([f64; 3], [Mat3; 3]) {
        let c = right_cauchy_green(f);
        let mut k = [0.0; 3];
        let mut a = [[[0.0; 3]; 3]; 3];
        for i in 0..3 {
            let a0 = [self.axes[0][i], self.axes[1][i], self.axes[2][i]];
            let mut ca0 = [0.0; 3];
            for r in 0..3 {
                for s in 0..3 {
                    ca0[r] += c[r][s] * a0[s];
                }
            }
            k[i] = a0[0] * ca0[0] + a0[1] * ca0[1] + a0[2] * ca0[2];
            let mut ai = [0.0; 3];
            for r in 0..3 {
                for s in 0..3 {
                    ai[r] += f[r][s] * a0[s];
                }
                ai[r] /= k[i].sqrt();
            }
            for r in 0..3 {
                for s in 0..3 {
                    a[i][r][s] = ai[r] * ai[s];
                }
            }
        }
        (k, a)
    }
}

impl ElasticMaterial for OrthoElastic {
    fn stress(&self, f: &Mat3) -> Mat3 {
        let j = det3(f);
        let b = left_cauchy_green(f);
        let (k, a) = self.texture(f);
        let mut bmi = b;
        for (i, row) in bmi.iter_mut().enumerate() {
            row[i] -= 1.0;
        }
        let mut s = [[0.0; 3]; 3];
        for i in 0..3 {
            // (A_i (b - I)).sym() * 2 mu_i K_i
            let ab = mul(&a[i], &bmi);
            for r in 0..3 {
                for c2 in 0..3 {
                    s[r][c2] += 0.5 * (ab[r][c2] + ab[c2][r]) * 2.0 * self.mu[i] * k[i];
                }
            }
            for jj in 0..3 {
                for r in 0..3 {
                    for c2 in 0..3 {
                        s[r][c2] += self.lam[i][jj]
                            * ((k[i] - 1.0) * k[jj] * a[jj][r][c2]
                                + (k[jj] - 1.0) * k[i] * a[i][r][c2])
                            / 2.0;
                    }
                }
            }
        }
        for row in s.iter_mut() {
            for v in row.iter_mut() {
                *v /= 2.0 * j;
            }
        }
        s
    }

    fn tangent(&self, f: &Mat3) -> Tens4 {
        let j = det3(f);
        let b = left_cauchy_green(f);
        let (k, a) = self.texture(f);
        let mut c = [[[[0.0; 3]; 3]; 3]; 3];
        for i in 0..3 {
            dyad4s(&a[i], &b, self.mu[i] * k[i], &mut c);
            for jj in 0..3 {
                dyad1s(&a[i], &a[jj], self.lam[i][jj] * k[i] * k[jj] / 2.0, &mut c);
            }
        }
        for x in c.iter_mut().flatten().flatten().flatten() {
            *x /= j;
        }
        c
    }
}

/// Build a material from a `.feb` material definition.
pub fn from_model(m: &Material) -> Result<Box<dyn ElasticMaterial>, String> {
    let p = |n: &str| {
        m.param(n)
            .ok_or_else(|| format!("material {:?} ({}): missing parameter {n}", m.name, m.kind))
    };
    match m.kind.as_str() {
        "neo-Hookean" => Ok(Box::new(NeoHookean {
            e: p("E")?,
            v: p("v")?,
        })),
        "orthotropic elastic" => {
            if m.tree.child("mat_axis").is_some() {
                return Err("orthotropic elastic: mat_axis is not supported yet".into());
            }
            Ok(Box::new(OrthoElastic::new(
                p("E1")?,
                p("E2")?,
                p("E3")?,
                p("G12")?,
                p("G23")?,
                p("G31")?,
                p("v12")?,
                p("v23")?,
                p("v31")?,
            )?))
        }
        k => Err(format!("material type '{k}' is not supported yet")),
    }
}
