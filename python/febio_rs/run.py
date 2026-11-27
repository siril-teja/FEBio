"""Run the febio4 executable and collect its results."""
import os
import shutil
import subprocess
import time

from . import log as _log
from . import qoi as _qoi
from .xplt import read_xplt

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def find_febio(febio=None):
    """Path of febio4: the argument, $FEBIO4, febio4 on PATH, or a build in the repository."""
    cands = [febio, os.environ.get("FEBIO4"), shutil.which("febio4"),
             os.path.join(_ROOT, "build", "bin", "febio4"), os.path.join(_ROOT, "cmbuild", "bin", "febio4")]
    for c in cands:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    raise FileNotFoundError("febio4 not found: pass febio=..., set FEBIO4, or build into build/bin")


class Result:
    """Outcome of one run. Files are read lazily."""

    def __init__(self, feb, returncode, wall, stdout, stderr):
        self.feb = feb
        stem = os.path.splitext(feb)[0]
        self.log_path = stem + ".log"
        self.xplt_path = stem + ".xplt"
        self.returncode = returncode
        self.wall = wall
        self.stdout = stdout
        self.stderr = stderr
        self._log = None

    @property
    def log(self):
        """Text of the log file ('' if there is none)."""
        if self._log is None:
            self._log = ""
            if os.path.exists(self.log_path):
                with open(self.log_path, encoding="latin-1") as f:
                    self._log = f.read()
        return self._log

    @property
    def ok(self):
        return _log.NORMAL in self.log

    def summary(self):
        """Termination, iteration counts and timings (see febio_rs.log.summary)."""
        return dict(_log.summary(self.log), wall=self.wall)

    def data(self):
        """Data blocks written to the log (see febio_rs.log.data)."""
        return _log.data(self.log)

    def plot(self, last_only=True):
        """The plot file (febio_rs.xplt.Plot)."""
        return read_xplt(self.xplt_path, last_only=last_only)

    def qoi(self):
        """Spine quantities of interest: range of motion, IDP, facet forces (febio_rs.qoi)."""
        return _qoi.qoi(self.feb, self.xplt_path, self.log)

    def __repr__(self):
        s = self.summary()
        return (f"Result({os.path.basename(self.feb)}: ok={s['ok']}, iterations={s['iterations']}, "
                f"total={s['time_total']} s)")


def run(feb, febio=None, solver=None, threads=None, deterministic=True, config=None,
        workdir=None, timeout=None, check=False):
    """Run febio4 on a .feb file.

    solver        default linear solver written to a config file ("rust", "pardiso", ...);
                  None uses the model's own setting (and no config file).
    threads       OMP_NUM_THREADS / MKL_NUM_THREADS (None: inherit).
    deterministic False passes -nondeterministic (upstream FEBio threading).
    config        path of a febio.xml to pass with -cnf (instead of solver).
    workdir       run in this directory (the .feb is copied there); default: its own.
    check         raise RuntimeError unless the run terminates normally.
    """
    exe = find_febio(febio)
    feb = os.path.abspath(feb)
    if workdir:
        os.makedirs(workdir, exist_ok=True)
        dst = os.path.join(os.path.abspath(workdir), os.path.basename(feb))
        if dst != feb:
            shutil.copyfile(feb, dst)
        feb = dst
    cwd = os.path.dirname(feb)
    args = [exe, "-i", os.path.basename(feb), "-silent"]
    if config:
        args += ["-cnf", os.path.abspath(config)]
    elif solver:
        cfg = os.path.splitext(feb)[0] + "_febio.xml"
        with open(cfg, "w") as f:
            f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_config version="3.0">\n'
                    f'<default_linear_solver type="{solver}"></default_linear_solver>\n</febio_config>\n')
        args += ["-cnf", cfg]
    else:
        args += ["-noconfig"]
    if not deterministic:
        args += ["-nondeterministic"]
    env = dict(os.environ)
    if threads:
        env.update(OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    t0 = time.time()
    p = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True, timeout=timeout)
    r = Result(feb, p.returncode, time.time() - t0, p.stdout, p.stderr)
    if check and not r.ok:
        raise RuntimeError(f"febio4 did not terminate normally on {feb}: {p.stderr[-500:]}")
    return r
