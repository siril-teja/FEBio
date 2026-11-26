#!/usr/bin/env python3
"""A/B two FEBio builds on the same models: wall-clock phase timings and bitwise output check.

usage: ab_builds.py --old PATH/febio4 --new PATH/febio4 --workdir DIR [--threads 4] [--reps 2]
"""
import argparse, json, os, shutil, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_bench import tsec, ival  # noqa: E402

PHASES = ["Total elapsed time", "time in linear solver", "model update", "evaluating stiffness",
          "evaluating residual", "reforming stiffness", "serialization", "IO-time (plot, dmp, data)"]
MODELS = [
    ("spine_n40_flex", "spine/gen_fsu.py", ["--mode", "flex", "--n", "40", "--nz", "12"]),
    ("spine_fric_n24_flex", "spine/gen_fsu.py", ["--mode", "flex", "--n", "24", "--nz", "8", "--friction", "0.1"]),
    ("hex20_n14_mr_newton", "gen_model.py", ["--n", "14", "--elem", "hex20", "--mat", "mr", "--newton"]),
    ("hex8_n20_newton", "gen_model.py", ["--n", "20", "--newton"]),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
    ap.add_argument("--workdir", required=True); ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--reps", type=int, default=2); ap.add_argument("--out", default="ab_builds.json")
    ap.add_argument("--solvers", default="pardiso,rust"); ap.add_argument("--models", default="all")
    a = ap.parse_args()
    a.old, a.new = os.path.abspath(a.old), os.path.abspath(a.new)
    res = []
    for name, gen, args in MODELS:
        if a.models != "all" and name not in a.models.split(","):
            continue
        for sol in a.solvers.split(","):
            stem = f"{name}_{sol}"
            for build in ("old", "new"):
                d = os.path.join(a.workdir, build); os.makedirs(d, exist_ok=True)
                subprocess.run([sys.executable, os.path.join(HERE, gen), *args, "--solver", sol,
                                "-o", os.path.join(d, stem + ".feb")], check=True, capture_output=True)
            for rep in range(a.reps):
                for build in ("old", "new"):
                    d = os.path.join(a.workdir, build)
                    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads), MKL_NUM_THREADS=str(a.threads))
                    subprocess.run([getattr(a, build), "-i", stem + ".feb", "-silent", "-noconfig"], cwd=d, env=env,
                                   capture_output=True)
                    log = open(os.path.join(d, stem + ".log")).read()
                    r = dict(model=name, solver=sol, build=build, rep=rep, ok="N O R M A L" in log,
                             iters=ival(log, "Total number of equilibrium iterations"))
                    for p in PHASES:
                        r[p] = tsec(log, p)
                    res.append(r)
                    json.dump(res, open(a.out, "w"), indent=1)
                    print(f"{name:22s} {sol:8s} {build} rep{rep} ok={r['ok']} iters={r['iters']} total={r[PHASES[0]]} "
                          f"lin={r[PHASES[1]]} update={r['model update']} stiff={r['evaluating stiffness']}", flush=True)
            same = open(os.path.join(a.workdir, "old", stem + "_disp.txt")).read() == \
                   open(os.path.join(a.workdir, "new", stem + "_disp.txt")).read()
            print(f"   -> final outputs identical: {same}", flush=True)

if __name__ == "__main__":
    main()
