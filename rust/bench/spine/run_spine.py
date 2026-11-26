#!/usr/bin/env python3
"""Spine FSU benchmark: PARDISO vs Rust solver on the L4-L5 surrogate.

Runs each (load case, mesh size, physics) with both solvers and thread counts, then reports
L4 rotation (ROM) at 7.5 N*m, iteration counts, solver/total time and the max nodal
displacement difference relative to PARDISO.

usage: run_spine.py --febio PATH/febio4 --workdir DIR [--cases ...] [--out results.json]
"""
import argparse, json, math, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from run_bench import tsec, ival, read_last_disp  # noqa: E402

LITERATURE = {  # Niosi et al. 2006 Eur Spine J 15:913, intact L3-L4, +-7.5 Nm, no preload, n=10
    "flex": (3.7, 1.5), "ext": (3.3, 1.5), "lb": (3.8, 1.4), "ar": (2.1, 0.9)}
AXIS = {"flex": 0, "ext": 0, "lb": 1, "ar": 2}

CASES = [
    # name, mode, n, nz, biphasic, solvers, threads
    ("fsu_n24_flex", "flex", 24, 8, False, [1, 4]),
    ("fsu_n24_ext", "ext", 24, 8, False, [1, 4]),
    ("fsu_n24_lb", "lb", 24, 8, False, [1, 4]),
    ("fsu_n24_ar", "ar", 24, 8, False, [1, 4]),
    ("fsu_n40_flex", "flex", 40, 12, False, [1, 4]),
    ("fsu_n40_ar", "ar", 40, 12, False, [1, 4]),
    ("fsu_n56_flex", "flex", 56, 16, False, [4]),
    ("fsu_bp_n24_flex", "flex", 24, 8, True, [1, 4]),
    ("fsu_bp_n40_flex", "flex", 40, 12, True, [4]),
]

def rom_deg(rb_file, mode):
    last = None
    for line in open(rb_file):
        if line.strip() and not line.startswith("*"):
            last = line.split()
    q = [float(v) for v in last[4:8]]
    return 2 * math.degrees(math.asin(q[AXIS[mode]]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--febio", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--cases", default="all")
    ap.add_argument("--out", default="spine_results.json")
    ap.add_argument("--solvers", default="pardiso,rust")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    res = json.load(open(a.out)) if os.path.exists(a.out) else []
    for name, mode, n, nz, bp, threads in CASES:
        if a.cases != "all" and name not in a.cases.split(","):
            continue
        ref = None
        refp = os.path.join(a.workdir, f"{name}_pardiso_disp.txt")
        if "pardiso" not in a.solvers.split(",") and os.path.exists(refp):
            ref = read_last_disp(refp)
        for sol in a.solvers.split(","):
            stem = f"{name}_{sol}"
            feb = os.path.join(a.workdir, stem + ".feb")
            g = [sys.executable, os.path.join(HERE, "gen_fsu.py"), "--mode", mode, "--n", str(n), "--nz", str(nz),
                 "--solver", sol, "-o", feb] + (["--biphasic"] if bp else [])
            subprocess.run(g, check=True, capture_output=True)
            for nt in threads:
                env = dict(os.environ, OMP_NUM_THREADS=str(nt), MKL_NUM_THREADS=str(nt))
                t0 = time.time()
                subprocess.run([a.febio, "-i", stem + ".feb", "-silent", "-noconfig"], cwd=a.workdir, env=env,
                               capture_output=True)
                wall = time.time() - t0
                log = open(feb[:-4] + ".log").read()
                r = dict(case=name, mode=mode, biphasic=bp, solver=sol, threads=nt, wall=wall,
                         ok="N O R M A L   T E R M I N A T I O N" in log,
                         neq=ival(log, "Nr of equations"), total=tsec(log, "Total elapsed time"),
                         linsolve=tsec(log, "time in linear solver"), iters=ival(log, "Total number of equilibrium iterations"),
                         refs=ival(log, "Total number of stiffness reformations"))
                if r["ok"]:
                    r["rom_deg"] = rom_deg(os.path.join(a.workdir, stem + "_rb.txt"), mode)
                    disp = read_last_disp(os.path.join(a.workdir, stem + "_disp.txt"))
                    if sol == "pardiso" and nt == max(threads):
                        ref = disp
                    if ref:
                        um = max(max(abs(v) for v in u) for u in ref.values())
                        r["max_diff_rel"] = max(max(abs(x - y) for x, y in zip(disp[k], ref[k])) for k in ref) / um
                res = [x for x in res if not (x["case"] == name and x["solver"] == sol and x["threads"] == nt)]
                res.append(r)
                json.dump(res, open(a.out, "w"), indent=1)
                print(f"{name:18s} {sol:8s} T={nt} neq={r['neq']} ok={r['ok']} total={r['total']} lin={r['linsolve']} "
                      f"iters={r['iters']} rom={r.get('rom_deg')} diff={r.get('max_diff_rel')}", flush=True)

if __name__ == "__main__":
    main()
