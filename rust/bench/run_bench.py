#!/usr/bin/env python3
"""FEBio linear-solver benchmark: runs each model with each solver / thread count,
parses timings + convergence from the FEBio log, and compares the final displacement
field against the reference solver (MKL PARDISO).

usage: run_bench.py --febio PATH/febio4 --workdir DIR [--quick]
"""
import argparse, json, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))

def tsec(log, label):
    m = re.search(re.escape(label) + r"[ .]*: [0-9:]+ \(([0-9.eE+-]+) sec\)", log)
    return float(m.group(1)) if m else None

def ival(log, label):
    m = re.search(re.escape(label) + r"[ .]*: *([0-9]+)", log)
    return int(m.group(1)) if m else None

def read_last_disp(path):
    """last step of a node_data logfile -> {node: (ux,uy,uz)}"""
    data, cur = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("*Step"):
                cur = {}
                data = cur
            elif line.startswith("*"):
                continue
            elif cur is not None and line.strip():
                p = line.strip().split(",")
                cur[int(p[0])] = tuple(float(x) for x in p[1:4])
    return data

def run(febio, feb, threads, env_extra=None):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    env["MKL_NUM_THREADS"] = str(threads)
    if env_extra: env.update(env_extra)
    t0 = time.time()
    p = subprocess.run([febio, "-i", os.path.basename(feb), "-silent", "-noconfig"],
                       cwd=os.path.dirname(feb), env=env, capture_output=True, text=True)
    wall = time.time() - t0
    log = open(feb[:-4] + ".log").read()
    ok = "N O R M A L   T E R M I N A T I O N" in log
    return dict(ok=ok, rc=p.returncode, wall=wall,
                total=tsec(log, "Total elapsed time"),
                linsolve=tsec(log, "time in linear solver"),
                stiff=tsec(log, "evaluating stiffness"),
                resid=tsec(log, "evaluating residual"),
                reform=tsec(log, "reforming stiffness"),
                iters=ival(log, "Total number of equilibrium iterations"),
                refs=ival(log, "Total number of stiffness reformations"),
                rhs=ival(log, "Total number of right hand evaluations"),
                neq=ival(log, "Nr of equations"))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--febio", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--cases", default="all")
    ap.add_argument("--out", default="results.json")
    ap.add_argument("--solvers", default="all", help="comma list to (re)run, e.g. rust")
    a = ap.parse_args()
    a.febio = os.path.abspath(a.febio)  # runs use the model directory as cwd
    os.makedirs(a.workdir, exist_ok=True)
    # (name, generator args, solvers, threads)
    cases = [
        ("hex8_n10_bfgs",   ["--n", "10"],                    ["skyline", "pardiso", "rust"], [1, 4]),
        ("hex8_n10_newton", ["--n", "10", "--newton"],        ["skyline", "pardiso", "rust"], [1, 4]),
        ("hex8_n20_bfgs",   ["--n", "20"],                    ["skyline", "pardiso", "rust"], [1, 4]),
        ("hex8_n20_newton", ["--n", "20", "--newton"],        ["pardiso", "rust"],            [1, 4]),
        ("hex8_n30_bfgs",   ["--n", "30"],                    ["pardiso", "rust"],            [1, 4]),
        ("hex8_n30_newton", ["--n", "30", "--newton"],        ["pardiso", "rust"],            [1, 4]),
        ("hex20_n10_newton",["--n", "10", "--elem", "hex20", "--newton"], ["pardiso", "rust"], [1, 4]),
        ("hex20_n14_mr_newton",["--n", "14", "--elem", "hex20", "--mat", "mr", "--newton"], ["pardiso", "rust"], [1, 4]),
        ("hex8_n40_newton", ["--n", "40", "--newton", "--steps", "3"], ["pardiso", "rust"],   [4]),
    ]
    if a.cases != "all":
        sel = a.cases.split(",")
        cases = [c for c in cases if c[0] in sel]
    results = []
    if os.path.exists(a.out):
        results = json.load(open(a.out))
    for name, gargs, solvers, threads in cases:
        ref = None
        refpath = os.path.join(a.workdir, f"{name}_pardiso_disp.txt")
        if os.path.exists(refpath):
            ref = read_last_disp(refpath)
        for sol in (["pardiso"] + [s for s in solvers if s != "pardiso"]):
            if sol not in solvers: continue
            if a.solvers != "all" and sol not in a.solvers.split(","): continue
            feb = os.path.join(a.workdir, f"{name}_{sol}.feb")
            subprocess.run([sys.executable, os.path.join(HERE, "gen_model.py"), *gargs, "--solver", sol, "-o", feb],
                           check=True, capture_output=True)
            for nt in threads:
                if sol == "skyline" and nt != 1:
                    continue  # skyline is single threaded
                r = run(a.febio, feb, nt)
                disp = read_last_disp(os.path.join(a.workdir, f"{name}_{sol}_disp.txt")) if r["ok"] else {}
                if sol == "pardiso" and nt == max(threads):
                    ref = disp
                if ref and disp:
                    umax = max(max(abs(v) for v in u) for u in ref.values())
                    dmax = max(max(abs(x - y) for x, y in zip(disp[k], ref[k])) for k in ref)
                    r["max_diff_rel"] = dmax / umax
                r.update(case=name, solver=sol, threads=nt)
                results = [x for x in results if not (x["case"] == name and x["solver"] == sol and x["threads"] == nt)]
                results.append(r)
                json.dump(results, open(a.out, "w"), indent=1)
                print(f"{name:22s} {sol:8s} T={nt}  neq={r['neq']}  total={r['total']}s  linsolve={r['linsolve']}s  "
                      f"iters={r['iters']} refs={r['refs']}  diff={r.get('max_diff_rel')}", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
