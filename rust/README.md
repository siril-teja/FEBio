# Rust linear solver for FEBio

`febio_rs_solver` is a parallel sparse direct solver written in Rust and plugged into
FEBio's `NumCore` as a regular `LinearSolver` (`NumCore/RustSolver.{h,cpp}`).
Select it per model with

```xml
<linear_solver type="rust"/>
```

Optional parameters: `nthreads` (default: `OMP_NUM_THREADS`), `ordering` (0 = METIS
nested dissection, 1 = AMD), `refine` (0 = iterative refinement only if pivots were
perturbed, n > 0 = always n steps), `pivot_eps` (default 1e-8), `print_level`.

## What it does

Symmetric matrices (FEBio `CompactSymmMatrix`, lower CSC == PARDISO's upper CSR):

1. **Analysis (once per sparsity pattern):** supervariable detection (identical
   column structure, e.g. the 3 dofs of a node), METIS nested-dissection ordering of
   the compressed graph, elimination-tree postorder, relaxed supernodes (faer's symbolic
   routines), and precomputed scatter maps so the numeric phase is pure arithmetic.
2. **Numeric LDLᵀ (every reformation):** left-looking supernodal factorization.
   Parallelism follows the Geist–Ng scheme: the elimination tree is cut into balanced
   independent subtrees that are factored concurrently (sequential dense kernels), and
   the remaining top of the tree is processed with multithreaded dense kernels (faer
   gemm/trsm/LDLᵀ) and column-block-parallel descendant updates.
   Tiny pivots are perturbed to ±1e-8·max|Aᵢⱼ| (PARDISO's strategy); iterative refinement
   is applied automatically when that happens.
3. **Solve:** supernodal forward/diagonal/backward substitution.

Unsymmetric matrices (`CRSSparseMatrix`): faer's supernodal sparse LU with partial
pivoting (not tuned; included so `symmetric_stiffness = 0` models still run).

Dependencies: [`faer`](https://crates.io/crates/faer) (dense kernels, symbolic Cholesky
helpers, sparse LU), `rayon`, and `metis-sys` (vendored METIS 5 C source, the same
ordering family PARDISO uses internally).

## Building

```sh
cmake . -B build -DUSE_RUST_SOLVER=ON [-DUSE_MKL=ON ...]
cmake --build build -j
```

CMake invokes `cargo build --release` and copies `libfebio_rs_solver.so` into the
build's `lib/` directory. The crate is compiled with `-C target-cpu=native`
(`rust/febio_rs_solver/.cargo/config.toml`); remove that for portable binaries.

The solver is linked as a shared library: when faer's x86 gemm backend is embedded as a
static archive inside another shared object (`libnumcore.so`) it crashes.

## Benchmarks

`rust/bench/gen_model.py` generates parametric hex8/hex20 block models (large-deformation
neo-Hookean / Mooney–Rivlin, compression + shear), `rust/bench/run_bench.py` runs them
with each solver and thread count and compares the final displacement fields against
MKL PARDISO. `rust/bench/pardiso_bench.c` and `examples/selftest.rs` / `examples/runmat.rs`
benchmark the factorization alone on a dumped matrix (`FRS_DUMP=file` makes the solver
write the matrix it factors).

See `rust/BENCHMARKS.md` for results.
