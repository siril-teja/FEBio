//! Solid elements (ROADMAP 4.3): shape functions, integration rules and kinematics for
//! hex8, tet4 and tet10, with FEBio's node ordering and default integration rules
//! (hex8: 2x2x2 Gauss, FE_HEX8G8; tet4: 1 point, FE_TET4G1; tet10: FEBio's 8-point
//! rule, FE_TET10G8). Matches FECore/FESolidElementShape.cpp and FEElementTraits.cpp.
use crate::model::ElemType;

/// A 3x3 matrix, row-major: m[i][j].
pub type Mat3 = [[f64; 3]; 3];

/// Shape functions and their parametric derivatives at the integration points of one
/// element type.
#[derive(Clone, Debug)]
pub struct ElementTraits {
    pub kind: ElemType,
    pub nodes: usize,
    /// integration points (r, s, t) and weights
    pub points: Vec<[f64; 3]>,
    pub weights: Vec<f64>,
    /// H[g][a]
    pub h: Vec<Vec<f64>>,
    /// dH/d(r,s,t): dh[g][a] = [Hr, Hs, Ht]
    pub dh: Vec<Vec<[f64; 3]>>,
}

/// Shape function values at (r, s, t).
pub fn shape(kind: &ElemType, r: f64, s: f64, t: f64) -> Vec<f64> {
    match kind {
        ElemType::Hex8 => vec![
            0.125 * (1.0 - r) * (1.0 - s) * (1.0 - t),
            0.125 * (1.0 + r) * (1.0 - s) * (1.0 - t),
            0.125 * (1.0 + r) * (1.0 + s) * (1.0 - t),
            0.125 * (1.0 - r) * (1.0 + s) * (1.0 - t),
            0.125 * (1.0 - r) * (1.0 - s) * (1.0 + t),
            0.125 * (1.0 + r) * (1.0 - s) * (1.0 + t),
            0.125 * (1.0 + r) * (1.0 + s) * (1.0 + t),
            0.125 * (1.0 - r) * (1.0 + s) * (1.0 + t),
        ],
        ElemType::Tet4 => vec![1.0 - r - s - t, r, s, t],
        ElemType::Tet10 => {
            let (r1, r2, r3, r4) = (1.0 - r - s - t, r, s, t);
            vec![
                r1 * (2.0 * r1 - 1.0),
                r2 * (2.0 * r2 - 1.0),
                r3 * (2.0 * r3 - 1.0),
                r4 * (2.0 * r4 - 1.0),
                4.0 * r1 * r2,
                4.0 * r2 * r3,
                4.0 * r3 * r1,
                4.0 * r1 * r4,
                4.0 * r2 * r4,
                4.0 * r3 * r4,
            ]
        }
        k => panic!("shape functions not implemented for {k:?}"),
    }
}

/// Parametric derivatives [Hr, Hs, Ht] of the shape functions at (r, s, t).
pub fn shape_deriv(kind: &ElemType, r: f64, s: f64, t: f64) -> Vec<[f64; 3]> {
    match kind {
        ElemType::Hex8 => {
            let c = [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ];
            c.iter()
                .map(|[a, b, g]| {
                    [
                        0.125 * a * (1.0 + b * s) * (1.0 + g * t),
                        0.125 * b * (1.0 + a * r) * (1.0 + g * t),
                        0.125 * g * (1.0 + a * r) * (1.0 + b * s),
                    ]
                })
                .collect()
        }
        ElemType::Tet4 => vec![
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        ElemType::Tet10 => vec![
            [
                -3.0 + 4.0 * r + 4.0 * (s + t),
                -3.0 + 4.0 * s + 4.0 * (r + t),
                -3.0 + 4.0 * t + 4.0 * (r + s),
            ],
            [4.0 * r - 1.0, 0.0, 0.0],
            [0.0, 4.0 * s - 1.0, 0.0],
            [0.0, 0.0, 4.0 * t - 1.0],
            [4.0 - 8.0 * r - 4.0 * (s + t), -4.0 * r, -4.0 * r],
            [4.0 * s, 4.0 * r, 0.0],
            [-4.0 * s, 4.0 - 8.0 * s - 4.0 * (r + t), -4.0 * s],
            [-4.0 * t, -4.0 * t, 4.0 - 8.0 * t - 4.0 * (r + s)],
            [4.0 * t, 0.0, 4.0 * r],
            [0.0, 4.0 * t, 4.0 * s],
        ],
        k => panic!("shape functions not implemented for {k:?}"),
    }
}

/// FEBio's default integration rule for an element type: points (r, s, t) and weights.
pub fn default_rule(kind: &ElemType) -> (Vec<[f64; 3]>, Vec<f64>) {
    match kind {
        ElemType::Hex8 => {
            let a = 1.0 / 3f64.sqrt();
            let p = vec![
                [-a, -a, -a],
                [a, -a, -a],
                [a, a, -a],
                [-a, a, -a],
                [-a, -a, a],
                [a, -a, a],
                [a, a, a],
                [-a, a, a],
            ];
            (p, vec![1.0; 8])
        }
        ElemType::Tet4 => (vec![[0.25, 0.25, 0.25]], vec![1.0 / 6.0]),
        ElemType::Tet10 => {
            // FE_TET10G8 (FEElementTraits.cpp), constants as given there
            let (a, b) = (0.0158359099, 0.3280546970);
            let (c, d) = (0.6791431780, 0.1069522740);
            let (w1, w2) = (0.138527967 / 6.0, 0.111472033 / 6.0);
            let p = vec![
                [a, b, b],
                [b, a, b],
                [b, b, a],
                [b, b, b],
                [c, d, d],
                [d, c, d],
                [d, d, c],
                [d, d, d],
            ];
            (p, vec![w1, w1, w1, w1, w2, w2, w2, w2])
        }
        k => panic!("integration rule not implemented for {k:?}"),
    }
}

impl ElementTraits {
    pub fn new(kind: &ElemType) -> ElementTraits {
        let (points, weights) = default_rule(kind);
        let h = points
            .iter()
            .map(|p| shape(kind, p[0], p[1], p[2]))
            .collect();
        let dh = points
            .iter()
            .map(|p| shape_deriv(kind, p[0], p[1], p[2]))
            .collect();
        ElementTraits {
            kind: kind.clone(),
            nodes: kind.nodes().unwrap(),
            points,
            weights,
            h,
            dh,
        }
    }

    pub fn gauss_points(&self) -> usize {
        self.points.len()
    }

    /// Jacobian dX/d(r,s,t) at integration point g for nodal positions x:
    /// J[i][j] = sum_a x_a[i] * dH_a/dxi_j
    pub fn jacobian(&self, x: &[[f64; 3]], g: usize) -> Mat3 {
        let mut j = [[0.0; 3]; 3];
        for (a, d) in self.dh[g].iter().enumerate() {
            for i in 0..3 {
                for k in 0..3 {
                    j[i][k] += x[a][i] * d[k];
                }
            }
        }
        j
    }

    /// Spatial shape-function gradients dH_a/dX at integration point g, and det(dX/dxi).
    pub fn gradients(&self, x0: &[[f64; 3]], g: usize) -> (Vec<[f64; 3]>, f64) {
        let j = self.jacobian(x0, g);
        let det = det3(&j);
        let ji = inv3(&j, det);
        // dH/dX_i = sum_k dH/dxi_k * dxi_k/dX_i = sum_k dH/dxi_k * Ji[k][i]
        let grads = self.dh[g]
            .iter()
            .map(|d| {
                let mut gr = [0.0; 3];
                for (i, gi) in gr.iter_mut().enumerate() {
                    *gi = d[0] * ji[0][i] + d[1] * ji[1][i] + d[2] * ji[2][i];
                }
                gr
            })
            .collect();
        (grads, det)
    }

    /// Deformation gradient F = dx/dX at integration point g, from reference positions
    /// x0 and current positions x.
    pub fn deformation_gradient(&self, x0: &[[f64; 3]], x: &[[f64; 3]], g: usize) -> Mat3 {
        let (grads, _) = self.gradients(x0, g);
        let mut f = [[0.0; 3]; 3];
        for (a, ga) in grads.iter().enumerate() {
            for i in 0..3 {
                for k in 0..3 {
                    f[i][k] += x[a][i] * ga[k];
                }
            }
        }
        f
    }

    /// Reference volume: sum_g det(J0_g) w_g.
    pub fn volume(&self, x0: &[[f64; 3]]) -> f64 {
        (0..self.gauss_points())
            .map(|g| det3(&self.jacobian(x0, g)) * self.weights[g])
            .sum()
    }
}

pub fn det3(m: &Mat3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// inverse of m, given its determinant
pub fn inv3(m: &Mat3, det: f64) -> Mat3 {
    let d = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ]
}
