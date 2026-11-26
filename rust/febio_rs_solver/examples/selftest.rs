//! Self test: 3-dof-per-node 27-point stencil on an n^3 grid (elasticity-like sparsity),
//! symmetric indefinite shift optional. Reports residual and timings.
use febio_rs_solver::{sym_spmv_lower, Solver};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let g: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(20);
    let nt: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(0);
    let shift: f64 = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(0.0);
    let ord: f64 = args.get(4).map(|s| s.parse().unwrap()).unwrap_or(0.0);
    let nn = g * g * g;
    let n = 3 * nn;
    let id = |i: usize, j: usize, k: usize| (k * g + j) * g + i;
    // lower CSC
    let mut cp = vec![0i32];
    let mut ri: Vec<i32> = vec![];
    let mut v: Vec<f64> = vec![];
    let mut seed = 12345u64;
    let mut rnd = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % 1000) as f64 / 1000.0
    };
    for k in 0..g {
        for j in 0..g {
            for i in 0..g {
                let me = id(i, j, k);
                for a in 0..3 {
                    let col = 3 * me + a;
                    let mut ent: Vec<(usize, f64)> = vec![];
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
                                for b in 0..3 {
                                    let row = 3 * o + b;
                                    if row < col {
                                        continue;
                                    }
                                    let val = if row == col {
                                        30.0 - shift
                                    } else {
                                        -0.4 - 0.3 * rnd()
                                    };
                                    ent.push((row, val));
                                }
                            }
                        }
                    }
                    ent.sort_by_key(|e| e.0);
                    for (r, x) in ent {
                        ri.push(r as i32);
                        v.push(x);
                    }
                    cp.push(ri.len() as i32);
                }
            }
        }
    }
    println!("n = {n}, nnz(lower) = {}", ri.len());
    if let Ok(path) = std::env::var("FRS_DUMP") {
        use std::io::Write;
        let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
        f.write_all(&(n as i64).to_le_bytes()).unwrap();
        f.write_all(&(ri.len() as i64).to_le_bytes()).unwrap();
        for x in &cp {
            f.write_all(&x.to_le_bytes()).unwrap();
        }
        for x in &ri {
            f.write_all(&x.to_le_bytes()).unwrap();
        }
        for x in &v {
            f.write_all(&x.to_le_bytes()).unwrap();
        }
        return;
    }
    let mut s = Solver::new(nt);
    s.ordering = match ord as i32 {
        1 => febio_rs_solver::Ordering::Amd,
        _ => febio_rs_solver::Ordering::Metis,
    };
    let t = Instant::now();
    s.analyze_symmetric(n, &cp, &ri, 0).unwrap();
    {
        let sy = s.symbolic().unwrap();
        let (c, t) = sy.critical_path();
        println!(
            "critical path {:.3e} of {:.3e} flops ({:.1}%)",
            c,
            t,
            100.0 * c / t
        );
    }
    println!(
        "analyze {:.3}s  nnz(L)={}  supernodes={} supervars={} flops={:.3e}",
        t.elapsed().as_secs_f64(),
        s.stats.nnz_l,
        s.stats.n_supernodes,
        s.stats.n_supervars,
        s.stats.flops
    );
    for rep in 0..3 {
        let t = Instant::now();
        s.factor(&v).unwrap();
        let tf = t.elapsed().as_secs_f64();
        println!(
            "factor[{rep}] {:.3}s  ({:.2} GFlop/s) perturbed={}",
            tf,
            s.stats.flops / tf * 1e-9,
            s.stats.n_perturbed
        );
    }
    let xt: Vec<f64> = (0..n)
        .map(|i| ((i * 7919) % 101) as f64 / 101.0 - 0.5)
        .collect();
    let cpu: Vec<usize> = cp.iter().map(|&x| x as usize).collect();
    let riu: Vec<usize> = ri.iter().map(|&x| x as usize).collect();
    let mut b = vec![0.0; n];
    sym_spmv_lower(n, &cpu, &riu, &v, &xt, &mut b);
    let mut x = vec![0.0; n];
    let t = Instant::now();
    s.solve(&mut x, &b).unwrap();
    let ts = t.elapsed().as_secs_f64();
    let mut r = vec![0.0; n];
    sym_spmv_lower(n, &cpu, &riu, &v, &x, &mut r);
    let rn = r
        .iter()
        .zip(&b)
        .map(|(a, c)| (a - c) * (a - c))
        .sum::<f64>()
        .sqrt();
    let bn = b.iter().map(|a| a * a).sum::<f64>().sqrt();
    let en = x
        .iter()
        .zip(&xt)
        .map(|(a, c)| (a - c) * (a - c))
        .sum::<f64>()
        .sqrt()
        / xt.iter().map(|a| a * a).sum::<f64>().sqrt();
    println!(
        "solve {:.4}s  rel residual {:.3e}  rel error {:.3e}",
        ts,
        rn / bn,
        en
    );
}
