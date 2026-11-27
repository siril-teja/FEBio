# Finley et al. (2018) FEBio lumbar spine models

Subject-specific FEBio models of the human lumbar spine, L1–L5, without the sacrum. They come from two female donors, aged 49 and 59.

> Finley SM, Brodke DS, Spina NT, DeDen CA, Ellis BJ. *FEBio finite element models of the human lumbar spine.* Computer Methods in Biomechanics and Biomedical Engineering 21(6):444–452, 2018. doi:10.1080/10255842.2018.1478967

- **Source:** SimTK project `febiolumbar` (https://simtk.org/projects/febiolumbar), files as downloaded.
- **License:** SimTK lists the models under **CC BY-NC 3.0**: non-commercial use only, with attribution. Check the project page for the current terms before redistributing.

## Files

| archive | inner archive | models |
|---|---|---|
| `Model149y.o.F-latest.zip` | `49F.zip` | `49F_L12/L23/L34/L45_{hex,hybrid}`; `49F_L15_{hex,hybrid,tet10}` |
| `Model259y.o.F-latest.zip` | `59F.zip` | `59F_L12/L23/L34/L45/L15_hybrid` |

- `Lxy` is a functional spinal unit; `L15` is the whole lumbar spine.
- Meshes: `hex` is all-hexahedral, `hybrid` mixes hex and tet4, and `tet10` uses quadratic tetrahedra.
- Most models ship with the authors' FEBio 2.x log, which gives iteration counts and timings.

### Model contents
- **Vertebrae:**
  - cortical shell: orthotropic elastic;
  - trabecular core and posterior elements: neo-Hookean;
  - bony end plates where loads and boundary conditions are applied: rigid bodies;
  - cartilaginous end plates: neo-Hookean.
- **Intervertebral disc:**
  - nucleus: neo-Hookean;
  - annulus: Holmes–Mow ground matrix with ±fiber-exp-pow lamellae.
- **Facets:** cartilage pads with `sliding2` or sliding-elastic contact.
- **Ligaments:** ALL, PLL, LF, ISL, SSL, ITL and capsular ligaments, as nonlinear springs.
- **Loading:**
  - 100 load steps;
  - the pressure loads ramp in over the first time unit;
  - a 7.5 N·m pure moment on the superior rigid body then ramps up smoothly;
  - the inferior body is fixed.

## Running them with FEBio 4

The files are `febio_spec 2.5`. FEBio 4, both upstream and this fork, cannot read them unchanged. `rust/bench/finley/convert_feb25.py` makes these rewrites. The first four change no values; the fifth is a documented assumption:

1. **Fiber angles:** `fiber-exp-pow` `<theta>/<phi>` become `<fiber type="angles">`. FEBio 4 evaluates these with the same formula as FEBio 2.x.
2. **Material axes:** a fiber component's `<mat_axis>` moves to its parent `solid mixture`. This is exact here: each mixture's fibers share one axis, and the matrix is isotropic. The script checks both and refuses otherwise.
3. **Contact name:** `sliding-tension-compression` becomes `sliding-elastic`. These are the same class under its FEBio 3 name.
4. **Module:** files that use `sliding2` declare `<Module type="biphasic"/>`. FEBio 4 registers `sliding2` only in that module. There are no biphasic materials, so the equation count is unchanged.
5. **Quasi-Newton method (assumption, 59F L1–L5 only):** `<qnmethod>2</qnmethod>` becomes BFGS. FEBio 4 reads 2 as JFNK, which needs an iterative linear solver; FEBio 2.5/2.6 had only BFGS and Broyden. BFGS matches the authors' reformation count. Both methods converge to the same rotation.

The benchmark driver is `rust/bench/finley/run_finley.py`. Results are in `rust/BENCHMARKS.md`.
