//! Recompute each element's mean J from a .feb mesh and nodal displacements, and compare
//! with FEBio's plotted "relative volume" (see rust/bench/core/jcheck.py).
//! usage: jcheck FILE.feb DUMP
use febio_core::element::{det3, ElementTraits};
use std::io::Read;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = febio_core::read_file(&a[1]).unwrap();
    let mut buf = vec![];
    std::fs::File::open(&a[2])
        .unwrap()
        .read_to_end(&mut buf)
        .unwrap();
    let mut o = 0usize;
    let i64_ = |o: &mut usize| {
        let v = i64::from_le_bytes(buf[*o..*o + 8].try_into().unwrap());
        *o += 8;
        v
    };
    let nn = i64_(&mut o) as usize;
    assert_eq!(nn, m.node_count(), "node count");
    let f64s = |o: &mut usize, n: usize| -> Vec<f64> {
        let v = (0..n)
            .map(|k| f64::from_le_bytes(buf[*o + 8 * k..*o + 8 * k + 8].try_into().unwrap()))
            .collect();
        *o += 8 * n;
        v
    };
    let d = f64s(&mut o, 3 * nn);
    let x: Vec<[f64; 3]> = (0..nn)
        .map(|i| {
            [
                m.nodes[i][0] + d[3 * i],
                m.nodes[i][1] + d[3 * i + 1],
                m.nodes[i][2] + d[3 * i + 2],
            ]
        })
        .collect();
    let nd = i64_(&mut o) as usize;
    // plot-file domains (element ids and mean J), matched to element blocks by element ids
    let mut doms: std::collections::HashMap<Vec<i64>, Vec<f64>> = Default::default();
    for _ in 0..nd {
        let ne = i64_(&mut o) as usize;
        let ids: Vec<i64> = (0..ne).map(|_| i64_(&mut o)).collect();
        let nj = i64_(&mut o) as usize;
        let j = f64s(&mut o, nj);
        if nj == ne {
            doms.insert(ids, j);
        }
    }
    let (mut worst, mut count, mut jmin, mut jmax) = (0.0f64, 0usize, f64::MAX, f64::MIN);
    let mut skipped = 0usize;
    for b in &m.blocks {
        if !b.kind.is_solid() {
            skipped += b.len();
            continue;
        }
        let jref = doms
            .get(&b.ids)
            .unwrap_or_else(|| panic!("no plot domain with the elements of block {:?}", b.name));
        let tr = ElementTraits::new(&b.kind);
        for (e, &jr) in jref.iter().enumerate() {
            let conn = b.element(e);
            let x0: Vec<[f64; 3]> = conn.iter().map(|&n| m.nodes[n]).collect();
            let xt: Vec<[f64; 3]> = conn.iter().map(|&n| x[n]).collect();
            let ng = tr.gauss_points();
            let j: f64 = (0..ng)
                .map(|g| det3(&tr.deformation_gradient(&x0, &xt, g)))
                .sum::<f64>()
                / ng as f64;
            worst = worst.max((j - jr).abs() / jr.abs().max(1e-12));
            jmin = jmin.min(j);
            jmax = jmax.max(j);
            count += 1;
        }
    }
    println!("{count} solid elements ({skipped} shell elements skipped), J in [{jmin:.6}, {jmax:.6}], max rel diff {worst:.3e}");
}
