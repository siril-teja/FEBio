# FEBio-Rust

Version: see [`VERSION`](VERSION) and [`CHANGELOG.md`](CHANGELOG.md). It is also shown in the `febio4` start-up banner.

FEBio-Rust is a fork of [FEBio](https://github.com/febiosoftware/FEBio), the nonlinear finite element solver for biomechanics. It keeps FEBio's input format, materials and physics, and adds a **parallel sparse direct linear solver written in Rust**. The focus is faster spine and soft-tissue models with **bitwise-reproducible results**.

## What is different from FEBio

- **Rust linear solver**: select it in any model with `<linear_solver type="rust"/>`.
  - Symmetric systems use a supernodal LDLᵀ factorization.
  - Unsymmetric systems (biphasic tissue, non-symmetric contact) use a supernodal LU.
  - It reuses its analysis when FEBio rebuilds a matrix with an unchanged sparsity pattern, which is common with contact.
- **Same answers as MKL PARDISO**: identical iteration counts and displacements to ~1e-12 on every benchmark model (block models and an L4–L5 spine segment).
- **Faster where it matters**: 1.3–2× faster linear solves on contact models at 1–4 threads, and 5–8× on biphasic (unsymmetric) disc models. Full numbers in [`rust/BENCHMARKS.md`](rust/BENCHMARKS.md).
- **Bitwise-reproducible runs (default)**: the same input gives byte-identical output on every run.
  - With the Rust solver, the output is also identical for any thread count.
  - Opt out with `-nondeterministic`; see [Determinism](#determinism).
- **CI-gated changes**: every change to `main` goes through a pull request. CI builds FEBio with both solvers and checks that results are unchanged and reproducible.

## Quick start

Requirements: CMake ≥ 3.16, a C++17 compiler, Rust (`cargo`). MKL is optional: it enables the PARDISO solver for comparison.

```sh
git clone https://github.com/siril-teja/FEBio-Rust.git
cd FEBio-Rust
cmake . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_RUST_SOLVER=ON
cmake --build build -j
./build/bin/febio4 -i model.feb
```

Use the Rust solver in a model:

```xml
<Control>
  <solver type="solid">
    <linear_solver type="rust"/>
  </solver>
</Control>
```

Optional solver parameters: `nthreads`, `ordering` (0 = METIS, 1 = AMD), `refine`, `pivot_eps`, `print_level`.

For an MKL build (PARDISO), see `ci/rust/setup-mkl.sh` and [BUILD.md](BUILD.md).

## Determinism

FEBio's element loops normally add contributions to the global matrix and vectors from several threads with atomic operations. The order of floating-point additions then varies from run to run, and so do the last bits of the results. FEBio-Rust removes that variation by default:

- **FEBio's own loops** (assembly, residual, material updates) run on one thread in deterministic mode.
- **The Rust solver** uses all requested threads (`OMP_NUM_THREADS` or `nthreads`). Its arithmetic is fixed regardless of scheduling, so the solution is bitwise identical run to run *and* for any thread count.
- **MKL PARDISO** runs in MKL's conditional-numerical-reproducibility mode (`MKL_CBWR_AUTO`), which makes it reproducible run to run for a fixed thread count.

| mode | how | result |
|---|---|---|
| deterministic (default) | — | byte-identical output on every run |
| non-deterministic | `febio4 -nondeterministic …` or `FEBIO_DETERMINISTIC=0` | FEBio's parallel atomic assembly as upstream; last-bit differences between runs |

Deterministic mode currently costs 17–90% wall time, because FEBio's element loops are serial; see `rust/BENCHMARKS.md`. Making those loops parallel *and* deterministic is the next item in `OPEN_ITEMS.md`.

## Repository layout

| path | contents |
|---|---|
| `rust/febio_rs_solver/` | the Rust solver crate (C ABI, unit and determinism tests) |
| `NumCore/RustSolver.*` | FEBio `LinearSolver` wrapper for the Rust solver |
| `rust/bench/` | benchmark models (block, L4–L5 spine segment), benchmark and regression scripts |
| `rust/BENCHMARKS.md` | methods and results |
| `OPEN_ITEMS.md` | known limitations and next steps |
| everything else | FEBio (upstream) |

## Development workflow

1. Branch from `main`, one change (one speed improvement) per branch.
2. Use [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) for commit messages and the PR title, e.g. `perf(pardiso): …`, `feat(solver): …`, `fix: …`. No AI attribution in titles, descriptions or commits.
3. Bump `VERSION` in the PR (and the version in `rust/febio_rs_solver/Cargo.toml`), and add a `CHANGELOG.md` entry:
   - `feat` → minor;
   - breaking (`!`) → major;
   - anything else → patch.
4. Open a pull request into `main`. It is squash-merged as one conventional commit.
5. CI must pass:
   - `pr-policy.yml`: conventional title and commits, no AI attribution, correct version bump, changelog entry;
   - `rust-solver-ci.yml`, as follows:
   - `cargo fmt`, `clippy` and the Rust tests;
   - an FEBio build with MKL and the Rust solver;
   - the regression gate: each model is run with PARDISO and with the Rust solver, and must produce identical iteration counts and displacements within 1e-9;
   - the determinism gate: repeated runs must produce byte-identical output (both solvers), and the Rust solver's output must be byte-identical on 1 and N threads.

## Status and limitations

See [`OPEN_ITEMS.md`](OPEN_ITEMS.md). In short:
- The published FEBio lumbar spine models (Finley et al. 2018) are not yet in the benchmark set.
- Saddle-point formulations (Lagrange-multiplier contact, mixed u–p) are not yet validated.
- Scaling was only measured up to 4 cores.

## Credits and license

FEBio is developed by the University of Utah and Columbia University (https://febio.org). Please cite FEBio as described there: Maas SA, Ellis BJ, Ateshian GA, Weiss JA. *FEBio: finite elements for biomechanics.* J Biomech Eng 134(1), 2012. [doi:10.1115/1.4005694](https://dx.doi.org/10.1115%2F1.4005694).

FEBio and this fork are released under the MIT license ([LICENSE](LICENSE)). The Rust solver uses [faer](https://codeberg.org/sarah-quinones/faer), [rayon](https://github.com/rayon-rs/rayon) and METIS (via `metis-sys`).
