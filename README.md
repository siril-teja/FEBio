# FEBio-Rust

FEBio-Rust is a fork of [FEBio](https://github.com/febiosoftware/FEBio), the nonlinear finite element solver for biomechanics. It keeps FEBio's input format, materials and physics, and adds a **parallel sparse direct linear solver written in Rust**. The focus is faster spine and soft-tissue models with **bitwise-reproducible results**.

## What is different from FEBio

- **Rust linear solver**: select it in any model with `<linear_solver type="rust"/>`.
  - Symmetric systems use a supernodal LDLᵀ factorization.
  - Unsymmetric systems (biphasic tissue, non-symmetric contact) use a supernodal LU.
  - It reuses its analysis when FEBio rebuilds a matrix with an unchanged sparsity pattern, which is common with contact.
- **Same answers as MKL PARDISO**: identical iteration counts and displacements to ~1e-12 on every benchmark model (block models and an L4–L5 spine segment).
- **Faster where it matters**: 1.3–2× faster linear solves on contact models at 1–4 threads, and 5–8× on biphasic (unsymmetric) disc models. Full numbers in [`rust/BENCHMARKS.md`](rust/BENCHMARKS.md).
- **CI-gated changes**: every change to `main` goes through a pull request. CI builds FEBio with both solvers and checks that results are unchanged.

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
2. Open a pull request into `main`.
3. CI (`.github/workflows/rust-solver-ci.yml`) must pass:
   - `cargo fmt`, `clippy` and the Rust tests;
   - an FEBio build with MKL and the Rust solver;
   - the regression gate: each model is run with PARDISO and with the Rust solver, and must produce identical iteration counts and displacements within 1e-9.

## Status and limitations

See [`OPEN_ITEMS.md`](OPEN_ITEMS.md). In short:
- The published FEBio lumbar spine models (Finley et al. 2018) are not yet in the benchmark set.
- Saddle-point formulations (Lagrange-multiplier contact, mixed u–p) are not yet validated.
- Scaling was only measured up to 4 cores.

## Credits and license

FEBio is developed by the University of Utah and Columbia University (https://febio.org). Please cite FEBio as described there: Maas SA, Ellis BJ, Ateshian GA, Weiss JA. *FEBio: finite elements for biomechanics.* J Biomech Eng 134(1), 2012. [doi:10.1115/1.4005694](https://dx.doi.org/10.1115%2F1.4005694).

FEBio and this fork are released under the MIT license ([LICENSE](LICENSE)). The Rust solver uses [faer](https://codeberg.org/sarah-quinones/faer), [rayon](https://github.com/rayon-rs/rayon) and METIS (via `metis-sys`).
