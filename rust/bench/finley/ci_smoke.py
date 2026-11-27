#!/usr/bin/env python3
"""CI smoke test on a published spine model (ROADMAP 2.3).

Runs the first --steps time steps of a Finley et al. (2018) model (default: 59F L4-L5
hybrid, a functional spinal unit) with MKL PARDISO and with the Rust solver, and checks:
  * both runs terminate normally with the same number of iterations and reformations;
  * the spine quantities of interest agree: range of motion within 1e-9, IDP and facet
    force within 1e-6 (relative to the group maximum; rust/bench/qoi.py);
  * with --repeat, the Rust run is reproducible (a repeat gives identical results; the
    regression suite already gates determinism, so CI leaves this off).

usage: ci_smoke.py --febio PATH/febio4 --workdir DIR [--model 59F_L45_hybrid]
                   [--steps 12] [--threads 2] [--repeat] [--md smoke.md]
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_finley as F  # noqa: E402
import qoi as Q  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--febio", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--model", default="59F_L45_hybrid")
    ap.add_argument("--steps", type=int, default=12)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--repeat", action="store_true")
    ap.add_argument("--md", default=None)
    a = ap.parse_args()
    if not F.ZIPS:
        sys.exit("no models/finley2018/*.zip found")
    os.makedirs(a.workdir, exist_ok=True)
    src = F.unpack(a.workdir)
    feb = F.prepare(src, a.workdir, a.model)
    # truncate the analysis to the first --steps time steps
    text = open(feb, encoding="latin-1").read()
    text, n = re.subn(r"<time_steps>\d+</time_steps>", f"<time_steps>{a.steps}</time_steps>", text)
    assert n == 1, "expected one <time_steps>"
    short = os.path.join(a.workdir, f"{a.model}_{a.steps}.feb")
    open(short, "w", encoding="latin-1").write(text)

    exe = os.path.abspath(a.febio)
    runs = {}
    todo = [("pardiso", 0), ("rust", 0)] + ([("rust", 1)] if a.repeat else [])
    for name, rep in todo:
        runs[(name, rep)] = F.run(a.workdir, short, name, exe, name, a.threads, rep)

    p, r = runs[("pardiso", 0)], runs[("rust", 0)]
    r2 = runs.get(("rust", 1), r)
    fails = []
    for x in (p, r, r2):
        if not x["ok"]:
            fails.append(f"{x['config']} rep {x['rep']} did not terminate normally: {x['stderr'][-300:]}")
    if not fails:
        if (p["iters"], p["refs"]) != (r["iters"], r["refs"]):
            fails.append(f"iterations/reformations differ: pardiso {p['iters']}/{p['refs']}, rust {r['iters']}/{r['refs']}")
        d = Q.compare(p["qoi"], r["qoi"])
        for k, tol in Q.TOL.items():
            if d.get(k) is not None and d[k] > tol:
                fails.append(f"QoI {k} differs by {d[k]:.3e} (tol {tol:g})")
        if r["quats"] != r2["quats"] or r["qoi"] != r2["qoi"]:
            fails.append("Rust run is not reproducible (results differ between repeats)")
    else:
        d = {}

    rom = max(r["rot_deg"].values()) if r.get("rot_deg") else float("nan")
    lines = [f"### Spine smoke test: Finley {a.model}, first {a.steps} steps ({a.threads} threads)", "",
             "| solver | iterations | reformations | linear solver (s) | total (s) |", "|---|---:|---:|---:|---:|"]
    for x in (p, r):
        lines.append(f"| {x['config']} | {x['iters']} | {x['refs']} | {x['linsolve']} | {x['total']} |")
    lines += ["", f"- max rotation (Rust): {rom:.4f}°", f"- QoI Δ (rel. to group max): {d}",
              f"- result: {'**PASS**' if not fails else '**FAIL**: ' + '; '.join(fails)}", ""]
    md = "\n".join(lines)
    print(md)
    if a.md:
        open(a.md, "w").write(md)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
