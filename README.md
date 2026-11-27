<div align="center">

# FEBio-Rust

**Fast, reproducible nonlinear finite elements for spine biomechanics.**

[![Rust solver CI](https://github.com/siril-teja/FEBio-Rust/actions/workflows/rust-solver-ci.yml/badge.svg)](https://github.com/siril-teja/FEBio-Rust/actions/workflows/rust-solver-ci.yml)
[![PR policy](https://github.com/siril-teja/FEBio-Rust/actions/workflows/pr-policy.yml/badge.svg)](https://github.com/siril-teja/FEBio-Rust/actions/workflows/pr-policy.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

[Roadmap](ROADMAP.md) · [Benchmarks](rust/BENCHMARKS.md) · [Changelog](CHANGELOG.md) · [Open items](OPEN_ITEMS.md)

</div>

FEBio-Rust is a solid-mechanics fork of [FEBio](https://github.com/febiosoftware/FEBio). It reads the same input files and gives the same answers, with a parallel linear solver written in Rust. The goal is the fastest open-source solver for spine models on bare-metal multicore machines; see the [roadmap](ROADMAP.md).

## Highlights

| | |
|---|---|
| **Rust sparse solver** | Supernodal LDLᵀ and LU, multithreaded. It reuses its analysis whenever the sparsity pattern is unchanged. |
| **Same answers** | Same iteration counts as MKL PARDISO. Displacements agree to about 1e-12. On spine models, range of motion, intradiscal pressure and facet contact forces are identical. |
| **Reproducible** | Byte-identical output on every run (default). With the Rust solver, also for any thread count. |
| **Gated** | Every change is a pull request; CI checks accuracy, reproducibility and versioning. |

## Quick start

Requires CMake ≥ 3.16, a C++17 compiler and Rust (`cargo`). MKL is optional; it adds the PARDISO solver.

```sh
git clone https://github.com/siril-teja/FEBio-Rust.git && cd FEBio-Rust
cmake . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_RUST_SOLVER=ON
cmake --build build -j
./build/bin/febio4 -i model.feb
```

Select the Rust solver in a model:

```xml
<solver type="solid">
  <linear_solver type="rust"/>
</solver>
```

Or use it for every model, via the config file passed with `-cnf`: `<default_linear_solver type="rust"/>`.

Solver options: `nthreads`, `ordering` (0 = METIS, 1 = AMD), `refine`, `pivot_eps`, `print_level`. For an MKL build, see [BUILD.md](BUILD.md) and `ci/rust/setup-mkl.sh`.

## Reproducibility

| Mode | Switch | Output |
|---|---|---|
| deterministic (default) | — | byte-identical on every run |
| non-deterministic | `-nondeterministic` or `FEBIO_DETERMINISTIC=0` | as upstream FEBio; last-bit differences between runs |

How the deterministic mode works:
- **Element loops:** elements are colored so that elements of one color share no node or rigid body. Each color is computed and assembled in parallel, the colors in a fixed order, so every global entry is summed in the same order for any thread count.
- **Rust solver:** its arithmetic does not depend on thread scheduling.
- **PARDISO:** runs in MKL's reproducible (CNR) mode.

## Scope

**Kept:** FECore, FEBioMech (solids, rigid bodies, contact, springs), FEBioMix (biphasic and the `sliding2` contact used by published spine models), FEBioXML, FEBioPlot, NumCore, FEAMR, FEBioLib, `febio4`.

**Removed (1.0.0):** fluid/FSI, multiscale (RVE), optimization, image analysis and material diagnostics. Use upstream FEBio for those.

## Repository

| Path | Contents |
|---|---|
| `rust/febio_rs_solver/` | Rust solver crate (C ABI, tests) |
| `NumCore/RustSolver.*` | FEBio wrapper for the Rust solver |
| `rust/bench/` | benchmark models, regression and benchmark scripts |
| `models/` | published spine models (Finley et al. 2018) used for benchmarks |
| everything else | FEBio (upstream) |

## Contributing

- **One change per pull request** into `main`, squash-merged.
- **Titles and commits:** [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).
- **Version bump:** bump `VERSION` (`feat` → minor, `!` → major, otherwise patch) and add a `CHANGELOG.md` entry.
- **CI must pass:**
  - the Rust crate checks;
  - the FEBio build;
  - the PARDISO-vs-Rust regression: iterations equal, displacements within 1e-9, and spine quantities of interest within tolerance (range of motion, intradiscal pressure, facet contact force);
  - the spine smoke test: the first 12 steps of the published Finley 59F L4–L5 model, PARDISO vs Rust;
  - the determinism gate;
  - the PR policy.

## Credits and license

FEBio is developed at the University of Utah and Columbia University ([febio.org](https://febio.org)). Please cite: Maas SA, Ellis BJ, Ateshian GA, Weiss JA. *FEBio: finite elements for biomechanics.* J Biomech Eng 134(1), 2012. [doi:10.1115/1.4005694](https://dx.doi.org/10.1115%2F1.4005694).

MIT license ([LICENSE](LICENSE)). The Rust solver builds on [faer](https://codeberg.org/sarah-quinones/faer), [rayon](https://github.com/rayon-rs/rayon) and METIS.
