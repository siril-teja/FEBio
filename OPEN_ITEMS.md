# Open items: Rust solver and FEBio performance work

Status as of 2026-11-26 (`main`). Background and results: `rust/README.md`, `rust/BENCHMARKS.md`.

## Validation / models

| # | Item | Why it matters | Next step |
|---|---|---|---|
| V1 | **Run the published FEBio lumbar spine models** (Finley et al. 2018, CMBBE 21(6):444–452; L1–L5, no sacrum; `49F.zip` / `59F.zip`) | The spine benchmark so far uses an idealized L4–L5 surrogate (`rust/bench/spine/gen_fsu.py`). The published models need a SimTK login, which the benchmark machine doesn't have. | Download from https://simtk.org/projects/febiolumbar, commit to e.g. `models/spine/` (CC BY-NC 3.0, keep attribution), then run `rust/bench/spine/run_spine.py`-style comparisons against the paper's ROM. |
| V2 | **Literature comparison at the right level** | The ROM corridor used (Niosi et al. 2006) is for L3–L4 and the surrogate is uncalibrated. Its ROM lies within about 1.3 SD (flexion within 1 SD). | Use L4–L5 in-vitro data at ±7.5 N·m (e.g. Heuer et al. 2007, from tables) once the published model is in. |
| V3 | **Saddle-point / indefinite systems** | The symmetric LDLᵀ uses static pivot perturbation without Bunch–Kaufman pivoting. It is not yet validated on Lagrange-multiplier contact or mixed u–p elements. | Add models with `laugon` = Lagrange multipliers and u–p (`three-field`) elements to `ci_regression.py`. Add intra-supernode 2×2 pivoting if needed. |
| V4 | **Solid-mechanics coverage** | The fork is solid-mechanics only (1.0.0). The regression set covers hex8/hex20 solids, fiber mixtures, rigid bodies, facet contact (frictionless and frictional) and springs. Shells, beams, prestrain and follower pressure are untested. | Add one small model each to `ci_regression.py`. |
| V5 | **FEBio's official test suite** | `ci/Linux/test.sh` expects `TestSuite/` (separate repo, not in this fork). | Obtain the test suite and run it with `linear_solver type="rust"` as default. |

## Performance

| # | Item | Notes |
|---|---|---|
| P1 | **Scaling beyond 4 cores** | All benchmarks ran on a 4-core VM. The top of the elimination tree relies on parallel dense kernels, so check 16–64 cores (and NUMA). |
| P2 | **Triangular-solve speed at 4 threads** | About 1.4× slower than PARDISO on the n30 block matrix (0.28 s vs 0.19 s). This matters for BFGS runs with many back-solves. The LU solve is sequential. |
| P3 | **PARDISO wrapper re-analyzes on every Factor()** | `NumCore/PardisoSolver.cpp` runs phase 11 (reordering + symbolic) at each reformation. Running it only on the first `Factor()` after `PreProcess()` measured −30% end to end on Newton runs without contact. It is a C++-only change. |
| P4 | **Non-solver time** | Profiled end-to-end with `perf`. Candidates:<br>• fiber-direction recomputation: ~9% of samples in fiber-reinforced spine models;<br>• OpenMP spin-wait competing with the Rust solver threads;<br>• the items P6–P9. |
| P11 | ~~Run-to-run determinism~~ **Done** | Deterministic mode is on by default: FEBio's OpenMP loops run on one thread; the Rust solver is bitwise reproducible and independent of thread count; PARDISO runs in MKL CNR mode. Enforced by the CI determinism gate. |
| P12 | **Deterministic mode cost** (partly done) | Solid and three-field domains now run material updates on all threads. Their element matrices and vectors are computed in parallel and assembled in element order, so results are bitwise identical to serial: −12% to −32% end to end. Still serial in deterministic mode: other domains (shells, beams, rigid, discrete), contact, and surface loads. Convert them the same way (`FECoreThreads::ComputeThenAssembleInOrder`). |
| P5 | **Unsymmetric fallback path** | The partial-pivoting fallback (faer LU) is slower than PARDISO. It never triggered in the tests, but it would on badly scaled systems. |
| P6 | **Static matrix profile rebuilt every time step** | `FENewtonSolver::ReformStiffness` calls `CreateStiffness(m_niter == 0)`, so the static (element) profile is rebuilt at the first iteration of every time step. That is 2–7% of spine runtime in `SparseMatrixProfile::ColumnProfile::insertRow`. It could be kept across time steps within an analysis step, but erosion (FEAMR `FEErosionAdaptor`) and mesh adaptors change connectivity between steps. A safe version needs a topology-change counter that those code paths increment. |
| P7 | **Restart snapshot each time step** | With `max_retries > 0`, `FEAnalysis::Solve` serializes the whole model into a `DumpMemStream` before every time step (`fem.Serialize(dmp)`). That is about 8–10% of spine runtime, dominated by per-material-point virtual `Serialize` calls. Options: a bulk (POD) snapshot of material-point arrays, or snapshotting only state that the step can modify. |
| P8 | **`dynamic_cast` in material-point access** | `FEMaterialPoint::ExtractData<T>()` walks the data chain with `dynamic_cast`. It is still a few % in `ElementInternalForce` and mixtures. A per-class type tag (enum or static id) would make it a pointer compare. |
| P9 | **Line searches** | In the spine runs, line searches (`FELineSearch::UpdateModel`) re-evaluate all material points and are a large share of the model-update time. This is algorithmic: changing `lstol` changes the iteration path, so it is a modeling choice, not a free speed-up. |
| P10 | **Rust vs. PARDISO at 4 threads (non-contact models)** | Once P3 is fixed, PARDISO is expected to be about 10–15% faster than the Rust solver at 4 threads on non-contact block models. Rust still wins at 1–2 threads, on contact models and on unsymmetric systems. To close the gap: top-of-tree parallel updates (multifrontal-style aggregation of small descendants) and the parallel solve (P2). |

## Engineering / CI

| # | Item | Notes |
|---|---|---|
| E1 | **Upstream FEBio workflows stay pending in this fork** | `*-pull-request.yml` call `febiosoftware/febio-workflows` with upstream secrets and runners. Decide whether to disable them in the fork. `Rust solver CI` is the active gate. |
| E2 | **MKL CMake integration** | `FindDependencies.cmake`/`CMakeLists.txt` only look for `lib*.a`. CI uses `ci/rust/setup-mkl.sh` (symlinks `.so` under `.a` names, MKL 2024.2 pinned because 2025 dropped `mkl_dcsrtrsv`). A proper `.so` search would be cleaner. |
| E3 | **Static linking of the Rust crate** | faer's x86 gemm crashes when its static archive is embedded in `libnumcore.so`, so the crate is a `cdylib`. Revisit with newer faer and file upstream. |
| E4 | **Portability** | The crate builds with `-C target-cpu=native`. Release binaries need a portable target (or runtime dispatch). Windows/macOS builds of the Rust path are untested. |
| E5 | **Default solver** | `rust` is opt-in (`<linear_solver type="rust"/>`). Making it the default needs V3–V5 first. |
