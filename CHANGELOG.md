# Changelog

All notable changes to FEBio-Rust. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Every pull request into `main` bumps `VERSION` according to its [conventional-commit](https://www.conventionalcommits.org/en/v1.0.0/) type:

| type | bump |
|---|---|
| breaking change (`!` or `BREAKING CHANGE:`) | major |
| `feat` | minor |
| everything else (`fix`, `perf`, `refactor`, `docs`, `test`, `build`, `ci`, `chore`, `revert`) | patch |

The FEBio version this fork is based on (`FEBioLib/version.h`, currently 4.13.0) is upstream's and is not changed here.

## [1.4.0] - 2026-11-27

### Added
- **Parallel parameter sweeps and load cases from Python** (ROADMAP 5.3).
  - `febio_rs.grid({path: [values]})` builds the cartesian product of `.feb` edits.
  - `febio_rs.sweep(feb, cases, workdir, threads, jobs)` runs one edited copy per case, `jobs` at a time with `threads` each (e.g. 8 × 4 on 32 cores), and returns the results in case order.
  - Every run is deterministic, so results do not depend on `jobs`. A test checks that parallel and sequential sweeps give identical plot data.

## [1.3.1] - 2026-11-27

### Build
- **Portable builds by default** (ROADMAP 3.9).
  - The Rust solver was built with `target-cpu=native` (`rust/febio_rs_solver/.cargo/config.toml`), so a binary could fail with an illegal instruction on a different CPU. It now builds for generic x86-64. faer's dense kernels select AVX2/AVX-512 at run time, so the same binary runs at full speed anywhere.
  - Measured on the Finley 59F L1–L5 matrix at 4 threads, the generic build gives a bitwise-identical solution (same hash). Factorization is 3.37–3.41 s vs 3.25–3.35 s native; solve is 0.158–0.162 s vs 0.164–0.168 s.
  - `-DRUST_TARGET_CPU=native` (or e.g. `x86-64-v3`) restores a host-tuned build.
  - CI no longer sets `target-cpu`.

## [1.3.0] - 2026-11-27

### Added
- **Scaling-study driver** `rust/bench/finley/scaling.py` (ROADMAP 2.4). For each model and thread count (default 4, 8, 16, 32) it runs native FEBio (PARDISO) and FEBio-Rust. It records total, linear-solver and phase times, iterations and QoI differences, and writes JSON plus a Markdown table. `--steps N` runs only the first N steps.

### Documentation
- OPEN_ITEMS P7 (restart snapshot) re-measured: 0.15 s per step, 1.7% of the full Finley 59F L1–L5 run.

## [1.2.1] - 2026-11-27

### Documentation
- **Roadmap: GPU items 3.13–3.15** for the target server (up to 32 cores and one NVIDIA L40).
  - 3.13: mixed-precision LU, with FP32 factorization and solves on the GPU and FP64 refinement on the CPU.
  - 3.14: optional cuDSS FP64 solver.
  - 3.15: overlap CPU assembly with GPU work.
- The speed target is now measured on up to 32 cores plus the L40.
- Items 3.2, 3.7 and 3.8 are marked done (#28, #22, #23).

## [1.2.0] - 2026-11-27

### Added
- **Python driver `febio_rs`** (ROADMAP 2.5), in `python/febio_rs`. Standard library only.
  - `run()` runs `febio4` with a chosen linear solver, thread count and determinism mode, and returns a `Result`.
  - `Result` gives the run summary (iterations, reformations, timings), the data blocks logged to the log file, the plot file (`read_xplt`) and the spine QoIs (range of motion, IDP, facet forces).
  - `set_values()` / `get_value()` edit `.feb` inputs for parameter studies.
  - CLI: `python -m febio_rs run model.feb --solver rust --qoi`.
  - The plot-file reader and QoI code moved from `rust/bench` into the package; the benchmark scripts import them through thin shims.
  - CI runs the package tests, including a check that the plot output is identical for 1 and 2 threads.

### Documentation
- The scaling study (ROADMAP 2.4, OPEN_ITEMS P1) now targets 4–32 cores instead of 4–64.

## [1.1.10] - 2026-11-27

### CI
- **Spine smoke test on a published model** (ROADMAP 2.3). `rust/bench/finley/ci_smoke.py` runs the first 12 time steps of the Finley et al. (2018) 59F L4–L5 functional spinal unit (49,971 equations, from `models/finley2018`) with PARDISO and with the Rust solver.
  - It checks that both terminate normally, with equal iterations and reformations.
  - Range of motion must agree within 1e-9, and IDP and facet forces within 1e-6.
  - `--repeat` adds a reproducibility check.
  - In CI it runs on 2 threads alongside the regression suite (about 1 min) and writes its table to the job summary.
## [1.1.9] - 2026-11-27

### Documentation
- **Benchmarks: full Finley 59F L1–L5 rerun with 1.1.8.** FEBio-Rust 884 s vs native 975 s: 1.10× (1.05× with 1.1.4). Range of motion, IDP and facet forces identical.
- **Linear-solver limits recorded:** factorization runs near dense-kernel peak, the solves at memory bandwidth. Mixed precision was rejected: refinement needs 7–8 solves per right-hand side.

## [1.1.8] - 2026-11-27

### Performance
- **The static matrix profile is reused across time steps** (OPEN_ITEMS P6).
  - Before, every reformation at the first iteration of a time step rebuilt the element ("static") part of the sparse-matrix profile. On the Finley models that is every reformation.
  - `FEGlobalMatrix::Create` now keeps it while `FEModel::StaticProfileKey()` is unchanged. The key is a hash of element connectivity and activity, nodal equation numbers, rigid-body attachment, node flags and, in `FEMechModel`, rigid-body equations.
  - Anything that changes the profile (erosion, mesh adaptation, a new equation numbering) changes the key and forces a rebuild. Models with linear constraints always rebuild.
  - Finley 59F L1–L5, first 6 steps, 4 threads: reformation time 4.14 s → 1.85 s. The matrix is unchanged (same nonzeros) and the output is identical to 1.1.7, value by value, for this model and all ten regression models.

## [1.1.7] - 2026-11-27

### Performance
- **Deterministic parallel surface loads and per-node solver loops** (ROADMAP 3.5).
  - **Surface-load vectors** (`FESurface::LoadVector`, e.g. pressure) use the colored assembly. The coloring cache moved from `FESolidDomain` to `FEMeshPartition`, so domains and surfaces share it.
  - **Per-node loops** run on all requested threads also in deterministic mode: displacement scatter, nodal positions, the biphasic pressure and velocity update, the nodal reaction loads, the increment update and the start-of-step node update. Each node writes only its own data.
  - **Springs and rigid-body coupling** stay serial. On Finley 59F L1–L5 (292 springs) they cost about 0.01 s per 6 steps, less than the threading overhead.
  - **Results:** the pressure-load summation order changed, so results differ in the last bits. They are still bitwise reproducible for any thread count.
  - **Finley 59F L1–L5, first 6 steps, 4 threads, alternated ×2 against 1.1.4:** residual 3.75/3.85 s → 3.26/3.45 s, model update 3.88/3.74 s → 3.51/3.69 s. The comparison includes the 1.1.6 contact change.

## [1.1.6] - 2026-11-27

### Performance
- **Deterministic parallel contact projection** (ROADMAP 3.4).
  - `sliding2`: the integration-point projection loop now runs in parallel. Its temporaries were shared between iterations; they are now per point.
  - `sliding-elastic`: the projection and contact-pressure loops already had OpenMP pragmas but ran on one thread in deterministic mode. They now use all requested threads.
  - Each integration point writes only its own data, so the output is unchanged. Plot-file data is identical to 1.1.4 for all ten regression models, and for the first 6 steps of Finley 59F L1–L5.
  - Finley 59F L1–L5, first 6 steps, 4 threads: contact update 0.69 s → 0.27 s.
  - Contact forces and stiffness stay serial. They cost 1.7 s of the full 59F run, too little to parallelize.

## [1.1.5] - 2026-11-27

### Documentation
- **Benchmarks: full Finley 59F L1–L5 rerun with 1.1.4.** FEBio-Rust 937 s vs native FEBio 981 s (1.05×; it was 0.74× with 1.0.0). Linear solver 636 s vs 709 s. Iterations, range of motion, IDP and facet forces are identical.

## [1.1.4] - 2026-11-27

### Performance
- **Graph-colored parallel element assembly** (ROADMAP 3.11) for solid, three-field and biphasic domains. It replaces the ordered pipeline (3.10), whose single assembling thread bounded the element loops.
  - Elements are grouped into blocks of up to 16 consecutive elements (fewer for small domains). Blocks are colored first-fit in index order so that blocks of one color share no node or rigid body.
  - Blocks run in parallel. A block waits only for the earlier-colored blocks it shares a node with, so there is no barrier between colors.
  - Every global entry receives its contributions in color order, so output is bitwise identical from run to run and for any thread count. It differs from earlier versions in the last bits, because the summation order changed.
  - Elements without rigid nodes no longer take the rigid-coupling lock.
  - Finley 59F L1–L5, first 20 steps, 4 threads: stiffness 16.5–17.1 s → 8.4 s (native FEBio: 7.9–8.7 s), residual 16.9–17.4 s → 11.0 s, total 182–194 s → 166–172 s. Native: 178–196 s. Iterations and rotation identical.
- **Rust solver `omp_blocktime` default 0 → 5 ms.** OpenMP threads now stay awake between the back-to-back per-domain loops of one assembly and still sleep before the Rust factorization.
- ROADMAP 3.12 is marked done: on the full 59F L1–L5 run, Rust linear-solver time is now below PARDISO's (644 s vs 700 s).

## [1.1.3] - 2026-11-27

### Performance
- **Rust LU factorization runs the elimination tree as a task graph** (ROADMAP 3.12, part 2). Every top-of-tree supernode starts as soon as its children are done, instead of the whole top running in sequence after all subtrees.
  - The factorization is unchanged, so results are bitwise identical.
  - On the Finley 59F L1–L5 matrix: 3.4–3.8 s instead of 4.0–4.3 s per factorization at 4 threads.
- **Tuning checked and rejected:**
  - A finer subtree split was slower.
  - AMD ordering roughly triples the factorization work.
  - Other METIS options increase fill.
  - The LU kernels run close to the dense-kernel peak on one core.

## [1.1.2] - 2026-11-27

### Performance
- **Rust LU triangular solves are now parallel and deterministic** (ROADMAP 3.12, part 1). They use the same elimination-tree method as the LDLᵀ path: chunked left-looking forward gather, then a top-down backward pass. Results are bitwise identical for any thread count.
  - On the Finley 59F L1–L5 matrix (170,478 equations), a solve takes 0.16–0.19 s instead of 0.40 s at 4 threads.
  - On the first 20 steps of that model: linear solver −22%, total −16%. Iterations, range of motion, IDP and facet forces are identical.
- **`FRS_DUMP`** now also dumps unsymmetric matrices, and the new `runlu` example times them.

## [1.1.1] - 2026-11-27

### CI
- **Linux only:** removed the upstream FEBio Linux, macOS and Windows workflows and the release-tag workflow. They needed `febiosoftware` infrastructure and never ran in this fork. Also removed their scripts: `ci/Linux`, `ci/macOS`, `ci/Windows`, the plugin scripts and `ci/febio.xml`.
- **`Rust solver CI`:**
  - runs once per pull request into `main`; there are no more duplicate push runs. A merge to `main` only rebuilds, without tests, to refresh the caches that later PRs restore;
  - skips documentation-only changes;
  - runs the crate and FEBio jobs in parallel;
  - caches compiled Rust dependencies, using a portable `x86-64-v3` target in CI, plus ccache C++ objects and the MKL download;
  - installs only `ccache` from apt.

## [1.1.0] - 2026-11-27

### Added
- **Spine quantities of interest**, compared between configurations:
  - range of motion (rotation of each rigid body);
  - intradiscal pressure (mean hydrostatic pressure per nucleus);
  - facet contact force (FEBio's total contact force per contact surface).
- **Code:**
  - `rust/bench/qoi.py` computes the quantities.
  - `rust/bench/xplt.py` is a minimal reader for FEBio `.xplt` plot files, including compressed states.
- **CI regression gate:** now fails when range of motion differs by more than 1e-9, or IDP or facet force by more than 1e-6, between PARDISO and the Rust solver. The spine surrogate writes `contact force` to its plot file for this.
- **Finley benchmark driver:** records the quantities and compares them against the first configuration.
- **`rust/BENCHMARKS.md`:** the quantities for all five Finley models are identical between native FEBio 4.13 and FEBio-Rust.

## [1.0.0] - 2026-11-27

### Breaking
- FEBio-Rust is now **solid-mechanics focused**. Removed modules: `FEBioFluid` (CFD, FSI, fluid–solutes, thermo- and polar fluids), `FEBioRVE`, `FEBioOpt`, `FEImgLib` and `FEBioTest` (material diagnostics, `-d`). Models that use them no longer load. `FEBioMix` is kept: published FEBio spine models (Finley et al. 2018) use its `sliding2` facet contact.
- Removed the `febio::RunMaterialTest` library API (it depended on `FEBioTest`).
- The `febio4 -d` (diagnose) and `-s` (optimize) options now stop with "Don't know how to do task".
- Removed the LEVMAR, NLOPT and FFTW build options (used only by the removed modules), `Documentation/` and the Xcode project.

### Changed
- Spine surrogate (`gen_fsu.py`): new `--friction mu` (frictional facets give a non-symmetric system; they engage in axial rotation and extension) and `--contact sliding2`. The CI regression set adds `fsu_ar_sliding2` and `fsu_friction_nonsym`.

## [0.5.1] - 2026-11-26

### Fixed
- **`convert_feb25.py`:** `<qnmethod>2</qnmethod>` (59F L1–L5) is mapped to BFGS. FEBio 4 reads 2 as JFNK, which upstream FEBio 4.13 rejects with a direct solver; FEBio 2.5/2.6 had only BFGS and Broyden. This is a documented assumption: it matches the authors' reformation count, and both methods converge to the same rotation.

### Documentation
- **`rust/BENCHMARKS.md`:** Finley L1–L5 results. FEBio-Rust is 1.07× native on 49F hex and 0.74× on 59F hybrid; the latter is a large non-symmetric LU, now ROADMAP 3.12. Results are identical to native.

## [0.5.0] - 2026-11-26

### Added
- **Models:** the published Finley et al. (2018) FEBio lumbar spine models are now in `models/finley2018/`, moved from the repository root, with citation, license and conversion notes.
- **`rust/bench/finley/convert_feb25.py`:** makes the FEBio 2.5 files readable by FEBio 4 with no value changes:
  - fiber angles become fiber generators;
  - material axes move from fiber components to their mixture;
  - the old `sliding-elastic` contact name is updated;
  - `sliding2` models get the biphasic module scope, keeping the symmetric stiffness of FEBio 2.x.
- **`rust/bench/finley/run_finley.py`:** benchmark driver.
- **`rust/BENCHMARKS.md`:** L4–L5 results. FEBio-Rust matches native FEBio 4.13's iterations and rotations exactly and is 0.92–1.35× its speed at 4 threads.

## [0.4.10] - 2026-11-26

### Performance
- **Rust solver, LU path:** iterative refinement now runs only while the relative residual is above 1e-10, and stops when a step does not halve it. A step that makes the residual worse is undone.
  - **Before:** refinement ran whenever the residual was above 1e-12, which is the rounding level for these systems.
  - **Published Finley 49F L4–L5 hybrid model:** 149 refinement solves were added on top of 97 regular ones. They are now gone. Linear solver time drops from 79.8 s to 46.3 s (PARDISO: 50.2 s), with the same iteration count and a final rotation identical to native FEBio to 12 digits.

## [0.4.9] - 2026-11-26

### Performance
- **Pipelined ordered assembly:** while one thread adds a chunk of element matrices/vectors to the global system in element order, the other threads compute the next chunk. Element buffers are reused across calls. Output is byte-identical. Stiffness evaluation is about 10% faster on the Finley 59F L4–L5 model at 4 threads; end to end, the difference is within run-to-run noise.
- **Rigid-body coupling:** elements without rigid nodes skip the rigid-body coupling, which saves two heap allocations per element and assembly.

## [0.4.8] - 2026-11-26

### Documentation
- `README.md`: shorter, restructured, with CI badges and links to the roadmap, benchmarks and changelog.

## [0.4.7] - 2026-11-26

### CI
- PR policy checks the pull request's current title and description, read from the GitHub API after a short wait, instead of the event snapshot. A newer run for the same PR cancels the older one. A description edited right after opening therefore no longer leaves a stale failed check. The rules themselves are unchanged.

## [0.4.6] - 2026-11-26

### Fixed
- Rust solver: the default pivot-perturbation threshold for LDLᵀ is now 1e-13·max|Aᵢⱼ| (was 1e-8), the same as for LU. On the published Finley 2018 L4–L5 model, the old default perturbed genuine pivots: rigid-body and penalty rows are much stiffer than soft tissue. The nonlinear solve then needed 60 instead of 20 iterations over the first 4 steps. Now the Rust solver matches PARDISO's iteration count. Results on the regression set are unchanged.

## [0.4.5] - 2026-11-26

### Documentation
- `ROADMAP.md`: objective, targets and the phased plan; each item is one pull request. Replaces upstream FEBio's roadmap.

## [0.4.4] - 2026-11-26

### CI
- FEBio build job: drop the runner image's third-party apt sources (packages.microsoft.com) before `apt-get update`, so an outage there no longer fails the regression gate.

## [0.4.3] - 2026-11-26

### Performance
- Rust solver: idle Intel OpenMP threads sleep instead of spinning while the Rust thread pool runs (`omp_blocktime`, default 0). Rust linear-solver time −3 to −15%; output byte-identical.

## [0.4.2] - 2026-11-26

### Performance
- Deterministic parallel element loops for the solid, three-field and biphasic domains:
  - material updates run on all threads;
  - element matrices and vectors are computed in parallel and assembled in element order.
  - −12 to −32% end to end; output byte-identical, also to a single-threaded run.

## [0.4.1] - 2026-11-26

### Performance
- Fiber materials memoize the reference fiber direction per material point when all material axes and fiber generators depend only on the reference geometry. Model update −3 to −18%, output byte-identical.

## [0.4.0] - 2026-11-26

First versioned release of FEBio-Rust. It covers everything merged so far.

### Added
- Rust sparse direct linear solver, `<linear_solver type="rust"/>` (#1):
  - supernodal LDLᵀ for symmetric systems and supernodal LU for unsymmetric ones;
  - METIS ordering, tree-parallel and dense-kernel-parallel factorization;
  - symbolic analysis reused across FEBio's matrix rebuilds.
- Deterministic mode, on by default (#4):
  - byte-identical output from run to run;
  - with the Rust solver, also identical for any thread count;
  - `-nondeterministic` / `FEBIO_DETERMINISTIC=0` opts out.
- CI gate: Rust tests, FEBio build with MKL and the Rust solver, PARDISO-vs-Rust regression, determinism checks (#1, #4).
- Project README and `OPEN_ITEMS.md` (#3).
- Versioning: `VERSION`, this changelog, the version in the `febio4` banner, and a CI policy check for PR titles, commits, version bumps and changelog entries.

### Changed
- PARDISO: reordering and symbolic factorization run once per sparsity pattern instead of at every reformation: −23 to −25% on Newton runs without contact (#5).

### Fixed
- CI no longer caches Rust build artifacts across runners. They are compiled for the runner's CPU (`-C target-cpu=native`), and reusing them on another CPU crashed the tests with SIGILL.
