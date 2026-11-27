# febio_rs: Python driver

Runs the `febio4` executable of this repository and reads its results. Standard library only (Python ≥ 3.8).

```python
import sys; sys.path.insert(0, "python")   # or: export PYTHONPATH=python
import febio_rs as fr

r = fr.run("model.feb", solver="rust", threads=4, check=True)
r.summary()        # ok, equations, iterations, reformations, timings (s)
r.data()           # data blocks written to the log: [{time, data, rows}]
p = r.plot()       # plot file, last state: p.vars, p.domains, p.surfaces, p.states
r.qoi()            # spine QoIs: range of motion (deg), IDP per nucleus, facet force per surface

# parameter study: edit a copy of the input
for c1 in (0.05, 0.1, 0.2):
    f = fr.set_values("model.feb", f"c1_{c1}.feb", {"Material/material[@name='nucleus']/c1": c1})
    print(c1, fr.run(f, solver="rust").qoi()["rom_deg"])
```

Parameter sweeps and load cases run side by side, e.g. 8 runs × 4 threads on 32 cores:

```python
cases = fr.grid({"Material/material[@name='nucleus']/c1": [0.05, 0.1, 0.2],
                 "Control/time_steps": [50, 100]})
results = fr.sweep("model.feb", cases, "runs", threads=4, jobs=8)   # Results in case order
```

- **Finding `febio4`:** `febio=` argument, then `$FEBIO4`, then `febio4` on `PATH`, then `build/bin` or `cmbuild/bin` in the repository.
- **`run` options:**
  - `solver` writes a config file with that default linear solver; without it the model's own setting is used.
  - `threads` sets `OMP_NUM_THREADS` and `MKL_NUM_THREADS`.
  - `deterministic=False` passes `-nondeterministic`.
  - `workdir` runs on a copy of the model in that directory.
- **QoIs:**
  - need `stress` and `contact force` in the model's plot file, and rigid-body rotations (`qx;qy;qz;qw`) in the log;
  - see `febio_rs/qoi.py`.
- **Command line:** `python -m febio_rs run model.feb --solver rust --threads 4 --qoi` prints the summary as JSON.
- **Tests:** `python3 -m unittest discover -s python/tests`. The run tests need `febio4`.
