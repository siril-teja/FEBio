//! Recompute each element's mean J and, where the material is supported, its mean Cauchy
//! stress from a .feb mesh and FEBio's nodal displacements, and compare with FEBio's values
//! (see rust/bench/core/jcheck.py; all inputs have 12 significant digits).
//! usage: jcheck FILE.feb DUMP
use febio_core::element::{det3, ElementTraits, Mat3};
use febio_core::material::from_model;
use std::collections::HashMap;
use std::io::Read;

struct Reader {
    buf: Vec<u8>,
    o: usize,
}

impl Reader {
    fn i64(&mut self) -> i64 {
        let v = i64::from_le_bytes(self.buf[self.o..self.o + 8].try_into().unwrap());
        self.o += 8;
        v
    }
    fn f64(&mut self) -> f64 {
        let v = f64::from_le_bytes(self.buf[self.o..self.o + 8].try_into().unwrap());
        self.o += 8;
        v
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = febio_core::read_file(&a[1]).unwrap();
    let mut buf = vec![];
    std::fs::File::open(&a[2])
        .unwrap()
        .read_to_end(&mut buf)
        .unwrap();
    let mut r = Reader { buf, o: 0 };
    // displacements by node id
    let nn = r.i64() as usize;
    let mut u: HashMap<i64, [f64; 3]> = HashMap::new();
    for _ in 0..nn {
        let id = r.i64();
        u.insert(id, [r.f64(), r.f64(), r.f64()]);
    }
    let x: Vec<[f64; 3]> = m
        .nodes
        .iter()
        .zip(&m.node_ids)
        .map(|(p, id)| {
            let d = u[id];
            [p[0] + d[0], p[1] + d[1], p[2] + d[2]]
        })
        .collect();
    // FEBio's J and stress (xx, yy, zz, xy, yz, xz) by element id
    let ne = r.i64() as usize;
    let mut ev: HashMap<i64, [f64; 7]> = HashMap::new();
    for _ in 0..ne {
        let id = r.i64();
        let mut v = [0.0; 7];
        for x in v.iter_mut() {
            *x = r.f64();
        }
        ev.insert(id, v);
    }
    let (mut jworst, mut count, mut jmin, mut jmax) = (0.0f64, 0usize, f64::MAX, f64::MIN);
    let (mut stress_elems, mut stress_scale, mut stress_diffs) = (0usize, 0.0f64, Vec::new());
    for b in m.blocks.iter().filter(|b| b.kind.is_solid()) {
        let tr = ElementTraits::new(&b.kind);
        // material of the block: 2.5 id or 4.0 name
        let mat = b
            .material
            .as_ref()
            .and_then(|k| {
                m.materials.iter().find(|x| {
                    x.id.map(|i| i.to_string()).as_deref() == Some(k.as_str())
                        || x.name.as_deref() == Some(k.as_str())
                })
            })
            .and_then(|x| from_model(x).ok());
        if mat.is_some() {
            stress_elems += b.len();
        }
        let ng = tr.gauss_points();
        for e in 0..b.len() {
            let theirs = ev[&b.ids[e]];
            let conn = b.element(e);
            let x0: Vec<[f64; 3]> = conn.iter().map(|&n| m.nodes[n]).collect();
            let xt: Vec<[f64; 3]> = conn.iter().map(|&n| x[n]).collect();
            let fs: Vec<Mat3> = (0..ng)
                .map(|g| tr.deformation_gradient(&x0, &xt, g))
                .collect();
            let j: f64 = fs.iter().map(det3).sum::<f64>() / ng as f64;
            jworst = jworst.max((j - theirs[0]).abs() / theirs[0].abs().max(1e-12));
            jmin = jmin.min(j);
            jmax = jmax.max(j);
            count += 1;
            if let Some(mat) = mat.as_ref() {
                let mut s = [[0.0; 3]; 3];
                for f in &fs {
                    let sg = mat.stress(f);
                    for (row, srow) in s.iter_mut().zip(sg.iter()) {
                        for (v, w) in row.iter_mut().zip(srow.iter()) {
                            *v += w / ng as f64;
                        }
                    }
                }
                let ours = [s[0][0], s[1][1], s[2][2], s[0][1], s[1][2], s[0][2]];
                let scale = theirs[1..].iter().fold(0.0f64, |a, v| a.max(v.abs()));
                let diff = ours
                    .iter()
                    .zip(&theirs[1..])
                    .fold(0.0f64, |a, (p, q)| a.max((p - q).abs()));
                stress_scale = stress_scale.max(scale);
                stress_diffs.push((diff, scale));
            }
        }
    }
    // stress: (a) largest difference relative to the peak stress of the model, and (b) the
    // largest relative difference over elements stressed above 1% of the peak (in nearly
    // unstressed elements, strains ~1e-9, the 12-digit displacements set the noise level)
    let speak = stress_diffs
        .iter()
        .fold(0.0f64, |a, (d, _)| a.max(d / stress_scale));
    let sworst = stress_diffs
        .iter()
        .filter(|(_, s)| *s > 1e-2 * stress_scale)
        .fold(0.0f64, |a, (d, s)| a.max(d / s));
    println!(
        "{count} solid elements, J in [{jmin:.6}, {jmax:.6}], J max rel diff {jworst:.3e}; \
         {stress_elems} elements with supported materials, stress max rel diff {sworst:.3e} (max diff / peak stress {speak:.3e})"
    );
}
