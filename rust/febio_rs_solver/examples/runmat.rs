//! Load a dumped symmetric matrix (lower CSC) and factor/solve it with the Rust solver.
use febio_rs_solver::{sym_spmv_lower, Solver};
use std::io::Read;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let nt: usize = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(0);
    let mut buf = vec![];
    std::fs::File::open(&a[1])
        .unwrap()
        .read_to_end(&mut buf)
        .unwrap();
    let rd8 = |o: usize| i64::from_le_bytes(buf[o..o + 8].try_into().unwrap());
    let n = rd8(0) as usize;
    let nnz = rd8(8) as usize;
    let mut o = 16;
    let cp: Vec<i32> = (0..=n)
        .map(|i| i32::from_le_bytes(buf[o + 4 * i..o + 4 * i + 4].try_into().unwrap()))
        .collect();
    o += 4 * (n + 1);
    let ri: Vec<i32> = (0..nnz)
        .map(|i| i32::from_le_bytes(buf[o + 4 * i..o + 4 * i + 4].try_into().unwrap()))
        .collect();
    o += 4 * nnz;
    let v: Vec<f64> = (0..nnz)
        .map(|i| f64::from_le_bytes(buf[o + 8 * i..o + 8 * i + 8].try_into().unwrap()))
        .collect();
    println!("n={n} nnz={nnz}");
    let mut s = Solver::new(nt);
    let t = Instant::now();
    s.analyze_symmetric(n, &cp, &ri, 0).unwrap();
    println!(
        "analyze {:.3}s nnz(L)={} sn={} flops={:.3e}",
        t.elapsed().as_secs_f64(),
        s.stats.nnz_l,
        s.stats.n_supernodes,
        s.stats.flops
    );
    for r in 0..2 {
        let t = Instant::now();
        s.factor(&v).unwrap();
        println!(
            "factor[{r}] {:.3}s perturbed={}",
            t.elapsed().as_secs_f64(),
            s.stats.n_perturbed
        );
    }
    let cpu: Vec<usize> = cp.iter().map(|&x| x as usize).collect();
    let riu: Vec<usize> = ri.iter().map(|&x| x as usize).collect();
    let xt: Vec<f64> = (0..n)
        .map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5)
        .collect();
    let mut b = vec![0.0; n];
    sym_spmv_lower(n, &cpu, &riu, &v, &xt, &mut b);
    let mut x = vec![0.0; n];
    let t = Instant::now();
    s.solve(&mut x, &b).unwrap();
    let ts = t.elapsed().as_secs_f64();
    let en = x
        .iter()
        .zip(&xt)
        .map(|(a, c)| (a - c) * (a - c))
        .sum::<f64>()
        .sqrt()
        / xt.iter().map(|a| a * a).sum::<f64>().sqrt();
    println!("solve {ts:.4}s rel error {en:.3e}");
}
