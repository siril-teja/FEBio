# Changelog

All notable changes to FEBio-Rust. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/).

Every pull request into `main` bumps `VERSION` according to its [conventional-commit](https://www.conventionalcommits.org/en/v1.0.0/) type:

| type | bump |
|---|---|
| breaking change (`!` or `BREAKING CHANGE:`) | major |
| `feat` | minor |
| everything else (`fix`, `perf`, `refactor`, `docs`, `test`, `build`, `ci`, `chore`, `revert`) | patch |

The FEBio version this fork is based on (`FEBioLib/version.h`, currently 4.13.0) is upstream's and is not changed here.

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
