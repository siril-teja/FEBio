# Changelog

All notable changes to FEBio-Rust. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Every pull request into `main` bumps `VERSION` according to its [conventional-commit](https://www.conventionalcommits.org/en/v1.0.0/) type:

| type | bump |
|---|---|
| breaking change (`!` or `BREAKING CHANGE:`) | major |
| `feat` | minor |
| everything else (`fix`, `perf`, `refactor`, `docs`, `test`, `build`, `ci`, `chore`, `revert`) | patch |

The FEBio version this fork is based on (`FEBioLib/version.h`, currently 4.13.0) is upstream's and is not changed here.

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
