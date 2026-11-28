# FEBio-Rust Roadmap

## Objective

Build an open-source nonlinear FE solver for spine biomechanics that is as fast as possible on bare-metal multicore machines. It must not give up accuracy or reproducibility to get there.

## Targets

| | Target | Measured on |
|---|---|---|
| Speed | ≥ 5× faster end to end than upstream FEBio 4.13 + PARDISO | Finley et al. 2018 L4–L5 and L1–L5 models, up to 32 cores + one NVIDIA L40 GPU |
| Accuracy | Same iteration counts; displacements and range of motion within 1e-9, intradiscal pressure and facet contact forces within 1e-6 (relative) of upstream FEBio | same models |
| Determinism | Byte-identical output for every run and every thread count | CI gate |
| Build | One `cargo build`, no vendor libraries | Linux, macOS, Windows |
| Use | Python API for scripting, parameter studies and calibration | — |

Completing every item below should reach these targets.

## How we work

- **One item = one pull request** into `main`, gated by CI: the regression check, the determinism check and the version and changelog policy.
- **Every speed item comes with benchmark numbers** and leaves the output unchanged.
- **Until the Rust core is complete, C++ FEBio is the reference:** every Rust-core item must match it.

Status: ✅ done · 🔄 in progress · ⬜ planned

## Phase 1 — Fast linear solver and reproducible runs (done)

| # | Item | Status |
|---|---|---|
| 1.1 | Parallel sparse direct solver in Rust (LDLᵀ / LU) | ✅ |
| 1.2 | CI gate: PARDISO vs Rust regression | ✅ |
| 1.3 | Bitwise-reproducible runs by default | ✅ #4 |
| 1.4 | PARDISO symbolic reuse | ✅ #5 |
| 1.5 | Versioning and PR policy | ✅ #7 |
| 1.6 | Fiber-direction memoization | ✅ #6 |
| 1.7 | Deterministic parallel element loops (solid, biphasic) | ✅ #8 |
| 1.8 | OpenMP idle threads sleep during the Rust solve | ✅ #9 |

## Phase 2 — Real models and baseline

| # | Item | Status |
|---|---|---|
| 2.1 | Finley 2018 lumbar models in `models/` + benchmark driver | ✅ |
| 2.2 | Trim to the solid mechanics we need (drop fluid, FSI, RVE, optimization) | ✅ #11 |
| 2.3 | CI smoke test on a published FSU model | ✅ |
| 2.4 | Scaling study, 4–32 cores and NUMA | ⬜ |
| 2.5 | Python driver for the `febio4` executable (run, parse log and plot files) | ✅ |
| 2.6 | Rust LDLᵀ pivot threshold for stiff/soft mixed systems | ✅ |
| 2.7 | Rust LU: refine only while it helps | ✅ |
| 2.8 | Spine quantities of interest (range of motion, intradiscal pressure, facet contact force) in benchmarks and the CI gate | ✅ |

## Phase 3 — Remove remaining bottlenecks in the C++ host

| # | Item | Status |
|---|---|---|
| 3.1 | Cheaper restart snapshot per time step | ⬜ |
| 3.2 | Reuse the stiffness sparsity pattern across time steps | ✅ #28 |
| 3.3 | Type-tag material point lookup (no `dynamic_cast`) | ⬜ |
| 3.4 | Deterministic parallel contact | ✅ |
| 3.5 | Deterministic parallel springs, rigid bodies and surface loads | ✅ |
| 3.6 | Faster contact search (spatial index) | ⬜ |
| 3.7 | Parallel triangular solve | ✅ #22 |
| 3.8 | Top-of-tree parallel factorization | ✅ #23 |
| 3.9 | Portable release builds (runtime CPU dispatch) | ✅ |
| 3.10 | Pipelined ordered assembly (compute next chunk while one thread assembles) | ✅ |
| 3.11 | Graph-colored parallel assembly (deterministic, thread-count independent) | ✅ |
| 3.12 | Faster Rust LU on large non-symmetric systems (≥ PARDISO) | ✅ |
| 3.13 | GPU mixed-precision LU: FP32 factorization and solves on the GPU (NVIDIA L40, 48 GB), FP64 refinement on the CPU, automatic FP64 fallback; opt-in, QoIs within tolerance | ⬜ |
| 3.14 | GPU FP64 direct solver (NVIDIA cuDSS) as an optional linear solver; compare with 3.13 and the 32-core CPU solver | ⬜ |
| 3.15 | Overlap CPU assembly with GPU factorization and solves | ⬜ |

## Phase 4 — Rust FE core

| # | Item | Status |
|---|---|---|
| 4.1 | Core: mesh, degrees of freedom, assembly, Newton/BFGS, time stepping | ⬜ |
| 4.2 | `.feb` reader (2.5 / 4.0) | ✅ |
| 4.3 | Elements: hex8, tet4, tet10 | ⬜ |
| 4.4 | Materials: neo-Hookean, orthotropic elastic | ⬜ |
| 4.5 | Materials: Holmes–Mow + fiber-exp-pow mixtures, local axes | ⬜ |
| 4.6 | Rigid bodies, nonlinear springs, pressure loads | ⬜ |
| 4.7 | Sliding-elastic contact | ⬜ |
| 4.8 | `sliding2` contact | ⬜ |
| 4.9 | Log and `.xplt` output (FEBio Studio compatible) | ⬜ |
| 4.10 | Python bindings (PyO3, NumPy arrays) | ⬜ |
| 4.11 | Gate: Finley models match C++ FEBio | ⬜ |

## Phase 5 — Beyond the C++ baseline

| # | Item | Status |
|---|---|---|
| 5.1 | SIMD and batched element kernels | ⬜ |
| 5.2 | NUMA-aware memory and thread placement | ⬜ |
| 5.3 | Parallel load cases and parameter sweeps from Python | ✅ |
| 5.4 | Rust core becomes the default; C++ kept as reference only | ⬜ |
| 5.5 | Target check: ≥ 5× on up to 32 cores + L40, accuracy and determinism held | ⬜ |
