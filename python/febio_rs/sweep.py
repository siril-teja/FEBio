"""Parameter sweeps and load cases: many febio4 runs side by side (ROADMAP 5.3)."""
import itertools
import os
from concurrent.futures import ThreadPoolExecutor

from .feb import set_values
from .run import run


def grid(params):
    """Cartesian product of parameter values: {path: [v1, v2, ...]} -> [ {path: v}, ... ]
    (last path varies fastest)."""
    keys = list(params)
    return [dict(zip(keys, vals)) for vals in itertools.product(*(params[k] for k in keys))]


def sweep(feb, cases, workdir, threads=1, jobs=None, solver="rust", febio=None, timeout=None):
    """Run one copy of `feb` per case, `jobs` at a time, each on `threads` threads.

    cases    list of {path: value} edits (see febio_rs.set_values), e.g. from grid();
             a case may also be {} (the unmodified model)
    workdir  case i runs in workdir/case_<i>/
    jobs     concurrent runs (default: CPU count // threads, at least 1)
    Returns the Results in case order. Each run is deterministic, so the results do not
    depend on jobs or on scheduling.
    """
    if jobs is None:
        jobs = max(1, (os.cpu_count() or 1) // max(1, threads))
    os.makedirs(workdir, exist_ok=True)
    name = os.path.basename(feb)

    def one(i):
        d = os.path.join(workdir, f"case_{i}")
        os.makedirs(d, exist_ok=True)
        dst = os.path.join(d, name)
        set_values(feb, dst, cases[i])
        return run(dst, febio=febio, solver=solver, threads=threads, timeout=timeout)

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(one, range(len(cases))))
