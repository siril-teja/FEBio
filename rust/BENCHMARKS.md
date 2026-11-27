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

# Published spine models: Finley et al. (2018)

These are the subject-specific FEBio lumbar models of Finley et al. (2018, CMBBE 21(6):444–452), stored in `models/finley2018/`.
- **Load case:** each model carries a pressure preload and then a 7.5 N·m pure moment over 100 load steps.
- **Format conversion:** the FEBio 2.5 files are converted for FEBio 4 by `rust/bench/finley/convert_feb25.py`. The converter changes no values, apart from the documented solver-setting mapping for 59F L1–L5 below, and produces identical input for every configuration; see `models/finley2018/README.md`.
- **Driver:** `rust/bench/finley/run_finley.py`; raw data in `rust/bench/results/finley_l45.json` (L4–L5) and `finley_l15.json` (L1–L5).
- **Machine:** 4 threads, the same 4-core VM as above. Run-to-run noise on this VM is about ±10%.

| config | executable | linear solver | mode |
|---|---|---|---|
| native | upstream FEBio 4.13 (`552138d`) | MKL PARDISO | upstream (non-deterministic) |
| Rust | FEBio-Rust (0.4.9 + #17, i.e. 0.4.10) | Rust | deterministic (default) |
| Rust + PARDISO | FEBio-Rust | MKL PARDISO (reproducible mode) | deterministic (default) |

## L4–L5 functional spinal units

| model | equations | native | Rust | Rust + PARDISO | Rust vs native | iterations / reformations | rotation of the loaded vertebra (°) |
|---|---:|---:|---:|---:|---:|---|---|
| 59F L4–L5 hybrid | 49,971 | 126.3 s (lin. 87.6) | **120.3 s** (lin. 64.1) | 152.7 s | **1.05×** | 378 / 100, all three | 4.0433, all three |
| 49F L4–L5 hex | 110,508 | 288.0 s (lin. 223.3) | **214.0 s** (lin. 134.6) | 342.4 s | **1.35×** | 368 / 100, all three | 4.0414, all three |
| 49F L4–L5 hybrid | 112,764 | 67.8 s (lin. 50.2) | 73.9 s (lin. 46.3) | 83.1 s | 0.92× | 97 / 15, all three | 5.0045, all three |

The 49F hybrid Rust row was measured with 0.4.10, which removed unnecessary LU refinement (#17); before 0.4.10 it took 101.4 s.

- **Accuracy:**
  - All three configurations take the same number of equilibrium iterations and reformations.
  - The final rigid-body rotations agree to the 12 digits FEBio prints.
  - The authors' FEBio 2.x logs report 378 iterations for 59F and 365 for 49F hex.
- **Linear solver:** the Rust solver is 8–40% faster than native PARDISO. FEBio-Rust's own PARDISO runs are slower than native, because MKL's reproducible (CNR) mode costs about 30%.
- **Remaining gap:**
  - **Stiffness and residual evaluation:** about 2× slower than native. FEBio-Rust adds element contributions in a fixed order, and native uses atomics (OPEN_ITEMS P13, ROADMAP 3.11).
  - **Serial loops:** contact and spring loops still run serially in deterministic mode (ROADMAP 3.4/3.5).
  - These costs dominate on BFGS-heavy runs with few reformations (49F hybrid).

## L1–L5 lumbar spines

| model | equations | native | Rust | Rust vs native | iterations / reformations | rotation of the loaded vertebra (°) |
|---|---:|---:|---:|---:|---|---|
| 49F L1–L5 hex | 280,467 | 1042 s (lin. 802) | **971 s** (lin. 623) | **1.07×** | 958 / 101, identical | 18.2622, identical |
| 59F L1–L5 hybrid | 170,478 | **703 s** (lin. 520) | 947 s (lin. 668) | 0.74× | 867 / 101, identical | 19.5712, identical |
| 59F L1–L5 hybrid, 1.1.4 (rerun) | 170,478 | 981 s (lin. 709) | **937 s** (lin. 636) | **1.05×** | 867 / 101, identical | 19.5712, identical |
| 59F L1–L5 hybrid, 1.1.8 (rerun) | 170,478 | 975 s (lin. 702) | **884 s** (lin. 642) | **1.10×** | 867 / 101, identical | 19.5712, identical |

Raw data is in `rust/bench/results/finley_l15.json`; the 1.1.4 and 1.1.8 reruns are in `finley_l15_59_v114.json` and `finley_l15_59_v118.json`. The rerun ran both builds back to back on a VM that was about 40% slower than for the earlier rows, so compare within a row only. Its IDP and facet forces are identical too (see below). The final rigid-body quaternions agree to every printed digit.

- **49F hex (symmetric matrix):** the Rust LDLᵀ is 22% faster than PARDISO. Stiffness, residual and model update are 1.5–2× slower than native, which keeps the end-to-end gain at 7%.
- **59F hybrid:** the model sets `symmetric_stiffness 0`, so it solves a non-symmetric system with LU. In this table (1.0.0) PARDISO's LU is faster: 101 factorizations of a matrix with 2.1·10⁸ nonzeros in L. Iterative refinement is not the cause: 4 refinement steps in 80 solves over the first 12 steps. ROADMAP 3.12 (1.1.2–1.1.3) has since made the Rust LU faster than PARDISO; see below.
- **59F hybrid after the LU work (1.1.3) and colored assembly (1.1.4):** timings vary by up to ±15% between runs on this VM, so builds were compared back to back.
  - **Full run, native then 1.1.3:** native 969 s (lin. 700), Rust 1046 s (lin. 644). The Rust LU now beats PARDISO, but stiffness, residual and model update cost 285 s against native's 166 s. The model also ran slower on this VM than for the table above (native 969 s vs 703 s).
  - **Full run, native then 1.1.4:** native 981 s (lin. 709, stiffness 41, residual 55, update 71), Rust 937 s (lin. 636, stiffness 41, residual 71, update 75): **1.05×**. Range of motion, IDP and facet forces identical.
  - **Full run, native then 1.1.8:** native 975 s, Rust 884 s (**1.10×**). Reformation 62 → 20 s (static profile reuse); residual 54 vs 65 s, the remaining gap. IDP and facet forces identical.
  - **Linear solver limits (1.1.8):** factorization is about 60% of Rust linear-solver time and runs near dense-kernel peak; the triangular solves (about 22%) stream the LU factors at memory bandwidth. Mixed precision was checked and rejected: f32 factors with f64 refinement need 7–8 solves to reach the accuracy of one f64 solve, and this BFGS run does about 8.6 solves per factorization.
  - **First 20 steps, alternated twice (s):**

    | build | total | linear solver | stiffness | residual | model update |
    |---|---:|---:|---:|---:|---:|
    | native (PARDISO) | 196 / 178 | 146 / 131 | 8.7 / 7.9 | 8.5 / 7.7 | 10.5 / 10.3 |
    | 1.1.3 (ordered assembly) | 194 / 182 | 122 / 113 | 17.1 / 16.5 | 17.4 / 16.9 | 13.7 / 13.0 |
    | **1.1.4 (colored assembly)** | **172 / 166** | 117 / 113 | **8.4 / 8.4** | **11.0 / 11.0** | 11.8 / 11.4 |

    128 iterations and identical rotation in every run. Colored assembly halves stiffness evaluation (now equal to native) and cuts the residual by 36%; end to end, 1.1.4 is 6–12% faster than native.
- **59F solver setting (assumption):** the model sets `qnmethod 2`. FEBio 4 reads that as JFNK, which needs an iterative solver; FEBio 2.5/2.6, which produced the authors' log, had only BFGS (0) and Broyden (1). The converter maps it to BFGS:

  | setting | iterations | reformations |
  |---|---:|---:|
  | authors' log (679 iterations) | 679 | 100 |
  | BFGS | 867 | 101 |
  | Broyden | 1121 | 112 |

  BFGS matches the authors' reformation count and was FEBio 2.x's default. Both settings converge to essentially the same rotation: 19.571° (BFGS) vs 19.568° (Broyden).

## Spine quantities of interest

The benchmark driver and the CI gate both compare three quantities between configurations, computed by `rust/bench/qoi.py` from the log and plot files:

| quantity | definition | tolerance |
|---|---|---|
| range of motion | rotation of each rigid body (vertebra) at the end of the run | 1e-9 relative |
| intradiscal pressure (IDP) | mean −(σxx+σyy+σzz)/3 over each nucleus's elements; positive = compressive | 1e-6 relative to the group maximum (float32 plot data) |
| facet contact force | magnitude of FEBio's total contact force on each contact surface | 1e-6 relative to the group maximum |

Native FEBio 4.13 against FEBio-Rust (Rust solver), 4 threads. The values are native's; FEBio-Rust's are **identical in every case**:

| model | range of motion | IDP per disc (MPa) | facet contact force (N) |
|---|---:|---|---|
| 59F L4–L5 hybrid | 4.0433° | L4–L5 −0.754 | no facet in contact |
| 49F L4–L5 hex | 4.0414° | L4–L5 −0.688 | L4–L5 right 0.447 |
| 49F L4–L5 hybrid | 5.0045° | L4–L5 0.108 | no facet in contact |
| 49F L1–L5 hex | 18.2622° | L1–L2 0.196, L2–L3 0.214, L3–L4 0.166, L4–L5 0.106 | L1–L2 R 1.627, L2–L3 R 0.282, L3–L4 R 0.895 |
| 59F L1–L5 hybrid | 19.5712° | L1–L2 −0.685, L2–L3 −0.654, L3–L4 −0.717, L4–L5 −0.682 | L1–L2 L 0.057, L2–L3 L 0.227, L3–L4 L 0.038 |

The negative IDP values mean the nucleus ends up in mean tension under this load. The 49F L4–L5 hybrid FEBio-Rust run is the pre-0.4.10 one; refinement changes since then do not change its rotation (identical to 12 digits).

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
