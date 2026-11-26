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

---

# Spine benchmark: L4–L5 functional spinal unit (surrogate)

The published Finley et al. (2018) L1–L5 FEBio models ([SimTK febiolumbar](https://simtk.org/projects/febiolumbar)) can only be downloaded with a SimTK login, so they could not be fetched from the benchmark machine. Instead `rust/bench/spine/gen_fsu.py` builds an **idealized L4–L5 FSU** with the ingredients that dominate solver cost in published FEBio spine models:

- **Vertebrae:** rigid bodies. L5 is fixed; L4 gets a ±7.5 N·m pure moment in 10 load steps.
- **Disc:** elliptic 50 × 34 × 10 mm, bonded to the vertebrae through rigid node sets.
  - nucleus: Mooney–Rivlin, c₁ = 0.12, c₂ = 0.03 MPa;
  - annulus: uncoupled mixture of Mooney–Rivlin (c₁ = 0.18, c₂ = 0.045 MPa) with two ±30° fiber families (`fiber-exp-pow-uncoupled`, ξ = 8 MPa, β = 2.5).
- **Facet joints:** left and right cartilage pads (neo-Hookean, E = 10 MPa) with frictionless sliding-elastic contact. Facets are oriented 45° to the sagittal plane with a 30° sagittal inclination.
- **Ligaments:** ALL, PLL, LF, ISL, SSL, ITL and CL as tension-only linear springs. Stiffnesses are representative (order of Pintar et al. 1992) and were not calibrated.
- **Biphasic variant (`--biphasic`):** the disc phases become biphasic with constant isotropic permeability, analysed transiently. This produces an **unsymmetric** system, which is what biphasic and multiphasic disc models do in FEBio.

Material and geometric parameters are literature-typical, **not calibrated**. The only adjustment was the facet sagittal inclination (0° → 30°): with vertical facets, extension came out at 11.2°.

## Range of motion vs. literature (4 threads, n24 mesh)

The literature values are from Niosi et al. 2006 (Eur Spine J 15:913): intact specimens, ±7.5 N·m, no preload, n = 10. They are for **L3–L4**, the closest level reported with this protocol in an openly readable table.

| motion | surrogate (PARDISO) | surrogate (Rust) | in vitro, mean ± SD |
|---|---:|---:|---:|
| flexion | 2.350° | 2.350° | 3.7 ± 1.5° |
| extension | 5.319° | 5.319° | 3.3 ± 1.5° |
| lateral bending (one side) | 5.567° | 5.567° | 3.8 ± 1.4° |
| axial rotation (one side) | 3.350° | 3.350° | 2.1 ± 0.9° |

All four motions fall within about 1.3 SD of the in-vitro means (flexion within 1 SD). Flexion ROM changes by only 0.5% from 15k to 148k dofs (2.350° → 2.360° → 2.362°). This is a sanity check of the surrogate, not a validation. The solver comparison does not depend on it.

## PARDISO vs. Rust on the spine surrogate

Times are in seconds. "lin." is FEBio's "time in linear solver". ROM is the L4 rotation at 7.5 N·m.

| model | dofs | T | PARDISO lin. | Rust lin. | speed-up | PARDISO total | Rust total | ROM P / R (deg) | iters P/R | max rel. Δu |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| fsu_n24_flex | 15,075 | 1 | 3.67 | 1.77 | **2.08×** | 9.0 | 7.1 | -2.350392 / -2.350392 | 36/36 | 5e-13 |
| fsu_n24_flex | 15,075 | 4 | 1.97 | 1.09 | **1.82×** | 4.4 | 3.4 | -2.350392 / -2.350392 | 36/36 | 5e-13 |
| fsu_n24_ext | 15,075 | 1 | 4.75 | 2.80 | **1.69×** | 13.0 | 10.1 | 5.319273 / 5.319273 | 70/70 | 2e-12 |
| fsu_n24_ext | 15,075 | 4 | 2.48 | 1.64 | **1.52×** | 5.8 | 4.6 | 5.319273 / 5.319273 | 70/70 | 2e-12 |
| fsu_n24_lb | 15,075 | 1 | 4.18 | 2.41 | **1.73×** | 11.0 | 9.3 | 5.567042 / 5.567042 | 57/57 | 2e-13 |
| fsu_n24_lb | 15,075 | 4 | 2.41 | 1.59 | **1.52×** | 5.4 | 4.6 | 5.567042 / 5.567042 | 57/57 | 2e-13 |
| fsu_n24_ar | 15,075 | 1 | 4.15 | 2.32 | **1.79×** | 11.4 | 9.5 | 3.349529 / 3.349529 | 56/56 | 4e-12 |
| fsu_n24_ar | 15,075 | 4 | 2.41 | 1.46 | **1.66×** | 5.6 | 4.5 | 3.349529 / 3.349529 | 56/56 | 4e-13 |
| fsu_n40_flex | 57,423 | 1 | 26.40 | 18.29 | **1.44×** | 49.2 | 41.2 | -2.359780 / -2.359780 | 43/43 | 5e-13 |
| fsu_n40_flex | 57,423 | 4 | 11.15 | 8.16 | **1.37×** | 21.5 | 17.8 | -2.359780 / -2.359780 | 43/43 | 5e-13 |
| fsu_n40_ar | 57,423 | 1 | 31.95 | 20.35 | **1.57×** | 62.3 | 47.1 | 3.368879 / 3.368879 | 54/54 | 4e-13 |
| fsu_n40_ar | 57,423 | 4 | 12.11 | 11.25 | **1.08×** | 24.9 | 22.6 | 3.368879 / 3.368879 | 54/54 | 4e-12 |
| fsu_n56_flex | 148,155 | 4 | 45.72 | 34.19 | **1.34×** | 74.5 | 61.0 | -2.362030 / -2.362030 | 43/43 | 5e-13 |
| fsu_bp_n24_flex (biphasic) | 20,700 | 1 | 27.38 | 6.35 | **4.31×** | 37.5 | 15.7 | -2.145760 / -2.145760 | 35/35 | 5e-13 |
| fsu_bp_n24_flex (biphasic) | 20,700 | 4 | 20.38 | 4.19 | **4.86×** | 24.8 | 8.3 | -2.145760 / -2.145760 | 35/35 | 5e-13 |
| fsu_bp_n40_flex (biphasic) | 79,276 | 4 | 268.5 | 31.8 | **8.44×** | 286.5 | 48.3 | -2.216770 / -2.216770 | 35/35 | 5e-13 |

**Numerical equivalence:** every run has identical equilibrium-iteration and reformation counts, and identical ROM to all 9 printed digits. Nodal displacement fields agree to within 4·10⁻¹² of max|u|, the precision floor of FEBio's logfile.

### What drives the spine speed-ups

1. **Unsymmetric systems (biphasic disc).** These go to the new supernodal LU (static pivoting with 10⁻¹³ perturbation, PARDISO's own default, plus residual check). PARDISO's unsymmetric mode (mtype 11) does weighted matching, scaling and pivoting, and FEBio reruns its ordering and symbolic analysis at every reformation. No pivot was perturbed, no partial-pivoting fallback was triggered, and the residuals were about 10⁻¹⁴ without refinement.
2. **Contact models rebuild the matrix profile every time step.** FEBio calls `Destroy()` + `PreProcess()` each time. The Rust solver keeps its symbolic analysis and compares patterns: an identical pattern is reused (9 of 10 rebuilds in `fsu_n40_flex`), and a pattern changed by less than 5% reuses the ordering. Before this change, symmetric spine runs at 4 threads were about 15% *slower* than PARDISO (12.7 s vs 11.1 s at n40).
3. **Faster symmetric factorization** at 1–2 threads, as in the block benchmarks above.

The smallest symmetric speed-up is `fsu_n40_ar` at 4 threads (1.08×). Axial rotation engages facet contact earlier, so the profile changes more often and some reformations need a full analysis.

---

# Determinism

**Sources of run-to-run variation found** (verified by repeated runs and byte comparison of the output):

1. **FEBio assembly.** The element loops (`#pragma omp parallel for`) add element matrices and vectors into the global system with `#pragma omp atomic`. The summation order follows thread timing.
2. **MKL PARDISO.** Multithreaded factorization without conditional numerical reproducibility (CNR).
3. **Rust solver.**
   - The parallel forward-solve gather used a `rayon` fold/reduce, whose split points depend on timing.
   - METIS keeps its random-number generator in global state, so two concurrent analyses in one process interleave draws.
   - The 1-thread solve used a different (right-looking) algorithm from the parallel one.
   - The symmetric update kernel split work differently for different thread counts.

**Fixes.**
- **Rust solver:**
  - the gather sums fixed-size chunks in a fixed order;
  - METIS calls are serialized;
  - one solve algorithm is used for every thread count;
  - the symmetric update uses a single gemm aligned to the update's first row.
  - The tests verify that faer's parallel gemm, LDLᵀ and triangular solve are bitwise identical to their sequential versions. The solver is now bitwise reproducible and **thread-count independent** (tests: `deterministic_*`, `thread_count_independent_results`, `tests/kernels.rs`).
- **FEBio:** deterministic mode, on by default (`FECore/fecore_threads.*`). FEBio's OpenMP loops run on one thread. Linear solvers keep the user's thread count. MKL runs with `mkl_cbwr_set(MKL_CBWR_AUTO)`. `-nondeterministic` / `FEBIO_DETERMINISTIC=0` restores the upstream behavior.

**Verification.** Spine FSU with facet contact, 4 threads:
- the Rust solver gives 1 distinct output over 3 runs, and byte-identical output for 1, 3 and 4 threads;
- PARDISO gives 1 distinct output over 3 runs;
- with `-nondeterministic`, 2 runs give 2 distinct outputs.

The CI regression gate repeats every model for both solvers, and runs the Rust solver on 1 thread, requiring byte-identical output.

**Cost of deterministic mode** (4 threads, mean of 2 runs, total wall time):

| model | solver | non-deterministic | deterministic | change |
|---|---|---:|---:|---:|
| spine FSU n40, flexion | Rust | 14.8 s | 25.9 s | +76% |
| spine FSU n40, flexion | PARDISO | 15.5 s | 29.0 s | +87% |
| biphasic FSU n24 | Rust | 6.7 s | 10.6 s | +58% |
| hex8 n20, Newton | Rust | 11.3 s | 13.2 s | +17% |
| hex8 n20, Newton | PARDISO | 12.4 s | 16.5 s | +33% |

The extra time is FEBio's element loops running serially. PARDISO's linear-solve time also rises (e.g. 8.2 → 10.8 s) because of MKL CNR. Recovering this without giving up determinism is OPEN_ITEMS P12.
