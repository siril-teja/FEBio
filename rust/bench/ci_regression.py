#!/usr/bin/env python3
"""CI regression gate for the Rust linear solver.

Runs a set of small FEBio models with MKL PARDISO (reference) and with the Rust solver,
and fails (exit 1) unless, for every model:
  * both runs terminate normally,
  * the number of equilibrium iterations and stiffness reformations are identical,
  * the final nodal displacement fields agree to max|du| / max|u| <= TOL,
  * determinism: repeating the run gives byte-identical output, for both solvers, and the
    Rust-solver output is byte-identical for 1 thread and for --threads threads.
Writes a JSON report and a markdown summary (for $GITHUB_STEP_SUMMARY).

usage: ci_regression.py --febio PATH/febio4 --workdir DIR [--threads 2] [--tol 1e-9]
"""
import argparse, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_bench import tsec, ival, read_last_disp  # noqa: E402

GEN = os.path.join(HERE, "gen_model.py")
FSU = os.path.join(HERE, "spine", "gen_fsu.py")
MODELS = [
    # name, generator, args, post-edit (sed-like replacement or None)
    ("hex8_bfgs", GEN, ["--n", "6"], None),
    ("hex8_newton_mr", GEN, ["--n", "6", "--newton", "--mat", "mr"], None),
    ("hex20_newton", GEN, ["--n", "4", "--elem", "hex20", "--newton"], None),
    ("hex8_nonsymmetric", GEN, ["--n", "6"],
     ("<symmetric_stiffness>symmetric</symmetric_stiffness>", "<symmetric_stiffness>non-symmetric</symmetric_stiffness>")),
    ("hex8_fiber_mixture", GEN, ["--n", "6", "--newton"],
     ('<material id="1" name="mat" type="neo-Hookean"><density>1</density><E>1</E><v>0.3</v></material>',
      '<material id="1" name="mat" type="solid mixture"><density>1</density>'
      '<solid type="neo-Hookean"><E>1</E><v>0.3</v></solid>'
      '<solid type="fiber-exp-pow"><ksi>2</ksi><alpha>0</alpha><beta>2.5</beta><fiber type="cylindrical">'
      '<center>0.5,0.5,0</center><axis>0,0,1</axis><vector>0,0.866,0.5</vector></fiber></solid>'
      '<solid type="fiber-exp-pow"><ksi>2</ksi><alpha>0</alpha><beta>2.5</beta><fiber type="vector">1,0,1</fiber></solid>'
      '</material>')),
    ("fsu_flex_contact", FSU, ["--mode", "flex", "--n", "12", "--nz", "4", "--steps", "5"], None),
    ("fsu_ar_contact", FSU, ["--mode", "ar", "--n", "12", "--nz", "4", "--steps", "5"], None),
    ("fsu_ar_sliding2", FSU, ["--mode", "ar", "--n", "12", "--nz", "4", "--steps", "5", "--contact", "sliding2"], None),
    ("fsu_friction_nonsym", FSU, ["--mode", "ar", "--n", "12", "--nz", "4", "--steps", "5", "--friction", "0.1"], None),
    ("fsu_biphasic", FSU, ["--mode", "flex", "--n", "12", "--nz", "4", "--steps", "5", "--biphasic"], None),
]

def run(febio, workdir, name, gen, args, edit, solver, threads, tag=""):
    stem = f"{name}_{solver}{tag}"
    feb = os.path.join(workdir, stem + ".feb")
    subprocess.run([sys.executable, gen, *args, "--solver", solver, "-o", feb], check=True, capture_output=True)
    if edit:
        txt = open(feb).read()
        assert edit[0] in txt
        open(feb, "w").write(txt.replace(edit[0], edit[1]))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    p = subprocess.run([febio, "-i", stem + ".feb", "-silent", "-noconfig"], cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=1800)
    logf = feb[:-4] + ".log"
    log = open(logf).read() if os.path.exists(logf) else ""
    r = dict(ok="N O R M A L   T E R M I N A T I O N" in log, rc=p.returncode,
             iters=ival(log, "Total number of equilibrium iterations"),
             refs=ival(log, "Total number of stiffness reformations"),
             linsolve=tsec(log, "time in linear solver"), total=tsec(log, "Total elapsed time"),
             neq=ival(log, "Nr of equations"))
    dfile = os.path.join(workdir, stem + "_disp.txt")
    r["disp"] = read_last_disp(dfile) if r["ok"] and os.path.exists(dfile) else {}
    r["raw"] = open(dfile, "rb").read() if os.path.exists(dfile) else b""
    if not r["ok"]:
        r["stderr"] = (p.stderr or "")[-2000:]
    return r

def fmt_diff(d):
    return 'n/a' if d is None else f'{d:.1e}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--febio", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--json", default="ci_regression.json")
    ap.add_argument("--md", default="ci_regression.md")
    a = ap.parse_args()
    a.febio = os.path.abspath(a.febio)  # runs use the model directory as cwd
    os.makedirs(a.workdir, exist_ok=True)
    rows, failed = [], []
    for name, gen, args, edit in MODELS:
        p = run(a.febio, a.workdir, name, gen, args, edit, "pardiso", a.threads)
        q = run(a.febio, a.workdir, name, gen, args, edit, "rust", a.threads)
        # determinism: repeat both runs, and the Rust run on one thread
        p2 = run(a.febio, a.workdir, name, gen, args, edit, "pardiso", a.threads, "_rep")
        q2 = run(a.febio, a.workdir, name, gen, args, edit, "rust", a.threads, "_rep")
        q1 = run(a.febio, a.workdir, name, gen, args, edit, "rust", 1, "_t1")
        diff = None
        if p["disp"] and q["disp"]:
            um = max(max(abs(v) for v in u) for u in p["disp"].values()) or 1.0
            diff = max(max(abs(x - y) for x, y in zip(q["disp"][k], p["disp"][k])) for k in p["disp"]) / um
        reasons = []
        if not p["ok"]: reasons.append("pardiso run failed")
        if not q["ok"]: reasons.append("rust run failed")
        if p["ok"] and q["ok"]:
            if (p["iters"], p["refs"]) != (q["iters"], q["refs"]):
                reasons.append(f"iterations/reformations differ {p['iters']}/{p['refs']} vs {q['iters']}/{q['refs']}")
            if diff is None or diff > a.tol:
                reasons.append(f"displacement difference {diff} > {a.tol}")
            if p2["raw"] != p["raw"]:
                reasons.append("PARDISO output differs between two identical runs")
            if q2["raw"] != q["raw"]:
                reasons.append("Rust output differs between two identical runs")
            if q1["raw"] != q["raw"]:
                reasons.append(f"Rust output differs between 1 and {a.threads} threads")
        det = (p2["raw"] == p["raw"], q2["raw"] == q["raw"], q1["raw"] == q["raw"])
        row = dict(model=name, neq=p["neq"], iters=(p["iters"], q["iters"]), refs=(p["refs"], q["refs"]), det=det,
                   linsolve=(p["linsolve"], q["linsolve"]), max_diff_rel=diff, passed=not reasons, reasons=reasons,
                   stderr=q.get("stderr") or p.get("stderr"))
        rows.append(row)
        if reasons:
            failed.append(name)
        print(f"{'PASS' if not reasons else 'FAIL'} {name:20s} neq={p['neq']} iters={row['iters']} det={det} "
              f"lin={row['linsolve']} diff={diff} {'; '.join(reasons)}", flush=True)
    for r in rows:
        r.pop("raw", None)
    json.dump(rows, open(a.json, "w"), indent=1)
    with open(a.md, "w") as f:
        f.write(f"### Rust solver regression ({a.threads} threads, tol {a.tol:g})\n\n")
        f.write("| model | dofs | iters P/R | refs P/R | lin. solver time P/R (s) | max rel. Δu | bitwise repeat P/R, Rust 1T=nT | result |\n|---|---:|---|---|---|---:|---|---|\n")
        for r in rows:
            lp, lr = r["linsolve"]
            f.write(f"| {r['model']} | {r['neq']} | {r['iters'][0]}/{r['iters'][1]} | {r['refs'][0]}/{r['refs'][1]} | "
                    f"{lp if lp is None else round(lp, 3)}/{lr if lr is None else round(lr, 3)} | "
                    f"{fmt_diff(r['max_diff_rel'])} | {'/'.join('yes' if d else 'NO' for d in r['det'])} | "
                    f"{'✅' if r['passed'] else '❌ ' + '; '.join(r['reasons'])} |\n")
    if failed:
        print("FAILED:", ", ".join(failed))
        return 1
    print("all regression models passed")
    return 0

if __name__ == "__main__":
    sys.exit(main())
