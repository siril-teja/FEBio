#!/usr/bin/env python3
"""Benchmark the published Finley et al. (2018) lumbar spine models.

Runs each model with several (executable, linear solver) configurations and records
wall time, FEBio's timing breakdown, iteration and reformation counts, and the spine
quantities of interest (rust/bench/qoi.py): range of motion (rotation of every rigid
body), intradiscal pressure per nucleus and facet contact force per contact surface.
Configurations after the first are compared with the first on those quantities.

The models are FEBio 2.5 files. They are unpacked from models/finley2018/*.zip and
rewritten for FEBio 4 by convert_feb25.py: purely syntactic changes, identical for every
configuration. The only other change is the log request: rigid-body kinematics are
written for all rigid bodies instead of body 1. That is output only and does not
affect the solution.

usage:
  run_finley.py --workdir DIR --config NAME=FEBIO4[:SOLVER] [--config ...]
                [--models 49F_L45_hex,...] [--threads 4] [--reps 1] [--out results.json]
SOLVER is the default linear solver written to a FEBio config file (default: pardiso).
"""
import argparse
import glob
import json
import math
import os
import re
import subprocess
import sys
import time
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, HERE)
from convert_feb25 import convert  # noqa: E402
sys.path.insert(0, os.path.join(HERE, ".."))
import qoi as Q  # noqa: E402

ZIPS = sorted(glob.glob(os.path.join(ROOT, "models", "finley2018", "*.zip")))


def unpack(workdir):
    """Unpack the nested archives (outer zip -> 49F.zip / 59F.zip -> .feb, .log)."""
    src = os.path.join(workdir, "src")
    os.makedirs(src, exist_ok=True)
    for z in ZIPS:
        with zipfile.ZipFile(z) as outer:
            for inner_name in outer.namelist():
                with outer.open(inner_name) as fh, zipfile.ZipFile(fh) as inner:
                    for name in inner.namelist():
                        base = os.path.basename(name)
                        if base and not os.path.exists(os.path.join(src, base)):
                            with inner.open(name) as fi, open(os.path.join(src, base), "wb") as fo:
                                fo.write(fi.read())
    return src


def prepare(src, workdir, model):
    out = os.path.join(workdir, model + ".feb")
    if not os.path.exists(out):
        text, *_ = convert(open(os.path.join(src, model + ".feb"), encoding="latin-1").read())
        text, n = re.subn(r'<rigid_body_data data="qx;qy;qz;qw">[^<]*</rigid_body_data>',
                          '<rigid_body_data data="qx;qy;qz;qw"></rigid_body_data>', text)
        assert n == 1, model
        open(out, "w", encoding="latin-1").write(text)
    return out


def tsec(log, key):
    m = re.search(re.escape(key) + r"[ .]*:\s*\d+:\d+:\d+ \((\d+(?:\.\d+)?) sec\)", log)
    return float(m.group(1)) if m else None


def ival(log, key):
    m = re.search(re.escape(key) + r"[ .]*:\s*(\d+)", log)
    return int(m.group(1)) if m else None


def last_quats(log):
    """Final 'qx;qy;qz;qw' record: {body id: (qx,qy,qz,qw)}."""
    blocks = re.findall(r"Data = qx;qy;qz;qw\n(.*?)\n\n", log, re.S)
    if not blocks:
        return {}
    q = {}
    for line in blocks[-1].strip().splitlines():
        v = line.split()
        q[int(v[0])] = tuple(float(x) for x in v[1:5])
    return q


def angle_deg(q):
    return 2.0 * math.degrees(math.atan2(math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2), q[3]))


def reference(src, model):
    p = os.path.join(src, model + ".log")
    if not os.path.exists(p):
        return None
    log = open(p, encoding="latin-1").read()
    return dict(iters=ival(log, "Total number of equilibrium iterations"),
                refs=ival(log, "Total number of stiffness reformations"),
                total=tsec(log, "Total elapsed time"), linsolve=tsec(log, "time in linear solver"))


def run(workdir, feb, cname, exe, solver, threads, rep):
    d = os.path.join(workdir, f"{cname}_r{rep}")
    os.makedirs(d, exist_ok=True)
    stem = os.path.splitext(os.path.basename(feb))[0]
    dst = os.path.join(d, stem + ".feb")
    if not os.path.exists(dst):
        os.symlink(feb, dst)
    cfg = os.path.join(d, "febio.xml")
    with open(cfg, "w") as f:
        f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<febio_config version="3.0">\n'
                f'<default_linear_solver type="{solver}"></default_linear_solver>\n</febio_config>\n')
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    t0 = time.time()
    p = subprocess.run([exe, "-i", stem + ".feb", "-cnf", "febio.xml", "-silent"], cwd=d, env=env,
                       capture_output=True, text=True)
    wall = time.time() - t0
    logp = os.path.join(d, stem + ".log")
    log = open(logp, encoding="latin-1").read() if os.path.exists(logp) else ""
    q = last_quats(log)
    return dict(model=stem, config=cname, solver=solver, threads=threads, rep=rep, wall=wall, rc=p.returncode,
                ok="N O R M A L   T E R M I N A T I O N" in log,
                neq=ival(log, "Nr of equations"),
                iters=ival(log, "Total number of equilibrium iterations"),
                refs=ival(log, "Total number of stiffness reformations"),
                rhs=ival(log, "Total number of right hand evaluations"),
                total=tsec(log, "Total elapsed time"), linsolve=tsec(log, "time in linear solver"),
                reform=tsec(log, "reforming stiffness"), stiff=tsec(log, "evaluating stiffness"),
                resid=tsec(log, "evaluating residual"), update=tsec(log, "model update"),
                quats={str(k): v for k, v in q.items()},
                rot_deg={str(k): angle_deg(v) for k, v in q.items()},
                qoi=Q.qoi(dst, os.path.join(d, stem + ".xplt"), log) if "N O R M A L" in log else None,
                stderr=p.stderr[-2000:] if not log else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--config", action="append", required=True, help="NAME=FEBIO4[:SOLVER]")
    ap.add_argument("--models", default="49F_L45_hex,59F_L45_hybrid")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--out", default="finley_results.json")
    a = ap.parse_args()
    if not ZIPS:
        sys.exit("no models/finley2018/*.zip found")
    os.makedirs(a.workdir, exist_ok=True)
    src = unpack(a.workdir)
    configs = []
    for c in a.config:
        name, rest = c.split("=", 1)
        exe, _, solver = rest.partition(":")
        configs.append((name, os.path.abspath(exe), solver or "pardiso"))
    res = json.load(open(a.out)) if os.path.exists(a.out) else []
    for model in a.models.split(","):
        feb = prepare(src, a.workdir, model)
        ref = reference(src, model)
        for rep in range(a.reps):
            for name, exe, solver in configs:
                r = run(a.workdir, feb, name, exe, solver, a.threads, rep)
                r["reference_2017"] = ref
                res = [x for x in res if not (x["model"] == model and x["config"] == name and x["rep"] == rep)]
                res.append(r)
                json.dump(res, open(a.out, "w"), indent=1)
                rot = ", ".join(f"{k}:{v:.4f}" for k, v in r["rot_deg"].items() if v > 1e-9)
                ref = next((x for x in res if x["model"] == model and x["rep"] == rep and x["config"] == configs[0][0]), None)
                if ref is not None and ref is not r and ref.get("qoi") and r.get("qoi"):
                    r["qoi_diff_vs_" + configs[0][0]] = Q.compare(ref["qoi"], r["qoi"])
                    json.dump(res, open(a.out, "w"), indent=1)
                    rot += f" qoi_diff={r['qoi_diff_vs_' + configs[0][0]]}"
                print(f"{model:18s} {name:14s} ok={r['ok']} wall={r['wall']:7.1f}s total={r['total']} "
                      f"lin={r['linsolve']} it={r['iters']} refs={r['refs']} rot[deg]={rot}", flush=True)


if __name__ == "__main__":
    main()
