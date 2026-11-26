# Rust solver benchmarks

## Setup

- **Machine:** 4-core Intel Xeon @ 2.1 GHz (AVX-512), 15 GB RAM, Linux VM (cloud container, not dedicated hardware; run-to-run noise ≈ ±5%).
- **FEBio:** this fork @ `552138d` (4.13), Release, GCC, OpenMP.
- **Baselines:** MKL PARDISO from MKL 2024.2 (the default solver in official FEBio builds) and Skyline (the default when FEBio is built without MKL).
- **Rust solver:** `febio_rs_solver` 0.1 (rustc 1.97, `-C target-cpu=native`, faer 0.24.4, METIS 5 vendored).
- **Threads:** `OMP_NUM_THREADS` = 1 or 4, used by both PARDISO and the Rust solver.
- **Models:** generated with `rust/bench/gen_model.py`. Each is a block [0,1]×[0,1]×[0,2] with the bottom face fixed and the top face given 25% compression plus 0.3 shear in 5 steps (3 steps for n40). That is large deformation.
  - `hex8_nN`: N×N×2N hex8 elements, neo-Hookean (E = 1, ν = 0.3).
  - `hex20_nN`: quadratic hex20 elements. `_mr` = nearly incompressible Mooney–Rivlin (κ/c₁ = 100).
  - `_bfgs`: FEBio defaults (BFGS, max_ups 10, so few factorizations and many back-solves).
  - `_newton`: full Newton (`max_ups 0`, so one factorization per iteration).
- **Driver:** `rust/bench/run_bench.py`. Raw data is in `rust/bench/results/`.

## Accuracy

Every run converged with **exactly the same number of equilibrium iterations, reformations and RHS evaluations** as PARDISO. The final nodal displacement fields agree to within **2·10⁻¹²** of max|u|. That is the floor set by FEBio's 12-significant-digit log output, so the agreement is at least that good. On the dumped n30 stiffness matrix the solve's relative error against a known solution is 1·10⁻¹⁴ for Rust and 3·10⁻¹⁵ for PARDISO. PARDISO always runs one iterative-refinement step under FEBio's settings. The Rust solver refines only when pivots were perturbed (`refine` parameter to force it).

## End-to-end FEBio: time in linear solver (s)

| model | dofs | T | skyline | PARDISO | Rust | speed-up vs PARDISO | iters (P/R) | max rel. Δu |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| hex8_n10_bfgs | 7 139 | 1 | 1.18 | 0.59 | 0.37 | **1.57×** | 18/18 | 2e-14 |
| hex8_n10_bfgs | 7 139 | 4 | – | 0.33 | 0.24 | **1.41×** | 18/18 | 2e-14 |
| hex8_n10_newton | 7 139 | 1 | 3.03 | 1.48 | 0.75 | **1.98×** | 15/15 | 2e-12 |
| hex8_n10_newton | 7 139 | 4 | – | 0.77 | 0.57 | **1.33×** | 15/15 | 2e-12 |
| hex8_n20_bfgs | 52 479 | 1 | 170.4 | 13.94 | 10.72 | **1.30×** | 18/18 | 2e-12 |
| hex8_n20_bfgs | 52 479 | 4 | – | 5.08 | 4.25 | **1.20×** | 18/18 | 2e-12 |
| hex8_n20_newton | 52 479 | 1 | – | 40.16 | 24.35 | **1.65×** | 15/15 | 2e-13 |
| hex8_n20_newton | 52 479 | 4 | – | 13.72 | 9.82 | **1.40×** | 15/15 | 2e-12 |
| hex8_n30_bfgs | 172 019 | 1 | – | 118.2 | 80.5 | **1.47×** | 18/18 | 2e-12 |
| hex8_n30_bfgs | 172 019 | 4 | – | 36.09 | 32.17 | **1.12×** | 18/18 | 2e-12 |
| hex8_n30_newton | 172 019 | 1 | – | 345.5 | 205.8 | **1.68×** | 15/15 | 2e-12 |
| hex8_n30_newton | 172 019 | 4 | – | 100.8 | 83.0 | **1.21×** | 15/15 | 2e-12 |
| hex20_n10_newton | 27 379 | 1 | – | 21.10 | 11.99 | **1.76×** | 15/15 | 2e-13 |
| hex20_n10_newton | 27 379 | 4 | – | 9.28 | 5.87 | **1.58×** | 15/15 | 2e-13 |
| hex20_n14_mr_newton | 72 435 | 1 | – | 136.6 | 82.0 | **1.67×** | 18/18 | 2e-12 |
| hex20_n14_mr_newton | 72 435 | 4 | – | 44.54 | 33.31 | **1.34×** | 18/18 | 2e-12 |
| hex8_n40_newton | 401 759 | 4 | – | 302.4 | 273.8 | **1.10×** | 10/10 | 2e-12 |

Skyline was run only up to n20, where it was already 12× slower than PARDISO. The two cheapest models were run once more on the final build and reproduced these timings within noise.

## Where the speed-up comes from

On the factorization kernel alone (the n30 FEBio stiffness matrix: 172k dofs, 6.7M nonzeros in the lower triangle; `runmat` vs `pardiso_bench`, both with METIS ordering):

| threads | Rust factor | PARDISO factor | Rust solve | PARDISO solve (incl. 1 refinement step) |
|---:|---:|---:|---:|---:|
| 1 | **12.8 s** | 18.7–19.1 s | 0.60 s | 0.56 s |
| 2 | **8.9 s** | 10.0–10.4 s | 0.41 s | 0.33 s |
| 4 | 5.25 s | 5.2–5.3 s | 0.28 s | 0.19 s |

| | Rust | PARDISO |
|---|---:|---:|
| nnz(L) | 179.8 M | 187.9 M |
| symbolic analysis | 2.0 s | 1.3–1.6 s |

So the end-to-end gains have two sources:

1. **Faster numerical factorization at low thread counts.** 1.47× at 1 thread and 1.15× at 2. At 4 threads both reach about 180 GFLOP/s, roughly two-thirds of this CPU's FP64 peak (~270 GFLOP/s), so they are at parity.
2. **Symbolic analysis is done once per sparsity pattern.** FEBio's `PardisoSolver::Factor()` reruns PARDISO's reordering and symbolic phase (phase 11) at every stiffness reformation. That costs about 1.4 s per reformation at n30, about 21 s of the 101 s in `hex8_n30_newton`. This accounts for most of the 4-thread advantage. The same saving is available in C++ by skipping phase 11 in `PardisoSolver` when the matrix profile hasn't changed.

The non-solver part of FEBio (assembly, residuals, material updates, I/O) is 5–25% of the runtime on the larger models. Porting the Newton driver or element assembly to Rust would therefore gain little: the linear solver is where the time is.

## Limitations

- Only 4 cores were available, and parallel scaling beyond 4 threads is untested. The top of the elimination tree, about 40–50% of the flops on these meshes, relies on parallel dense kernels.
- At 4 threads the triangular solve is about 1.4× slower than PARDISO's. That matters only for BFGS runs with many back-solves.
- Symmetric matrices use LDLᵀ with static pivot perturbation (1e-8·max|A|) and no Bunch–Kaufman pivoting. That is fine for displacement-based solid mechanics, and auto-refinement kicks in when pivots are perturbed. It has not been validated on saddle-point systems such as Lagrange-multiplier contact or mixed u–p formulations.
- The unsymmetric path (faer sparse LU) works and matches PARDISO, but it is slower than PARDISO (biphasic and fluid models were not benchmarked).
- The benchmark models are synthetic blocks. Contact, multiphysics and FEBio's official test-suite models were not run.
