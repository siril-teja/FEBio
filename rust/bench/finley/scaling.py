#!/usr/bin/env python3
"""Scaling study (ROADMAP 2.4): native FEBio vs FEBio-Rust over thread counts.

For each model and thread count, runs native FEBio (PARDISO) and FEBio-Rust (Rust solver),
records total and linear-solver time, iterations, and the spine QoI differences, and
writes a JSON file plus a Markdown table.

usage: scaling.py --native PATH/febio4 --rust PATH/febio4 --workdir DIR
                  [--models 59F_L15_hybrid,49F_L15_hex] [--threads 4,8,16,32]
                  [--steps N] [--out scaling.json] [--md scaling.md]
--steps truncates the analysis to its first N time steps (quick checks).
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_finley as F  # noqa: E402
import qoi as Q  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", required=True)
    ap.add_argument("--rust", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--models", default="59F_L15_hybrid")
    ap.add_argument("--threads", default="4,8,16,32")
    ap.add_argument("--steps", type=int, default=0)
    ap.add_argument("--out", default="scaling.json")
    ap.add_argument("--md", default="scaling.md")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    src = F.unpack(a.workdir)
    res = json.load(open(a.out)) if os.path.exists(a.out) else []
    for model in a.models.split(","):
        feb = F.prepare(src, a.workdir, model)
        if a.steps:
            t = open(feb, encoding="latin-1").read()
            t = re.sub(r"<time_steps>\d+</time_steps>", f"<time_steps>{a.steps}</time_steps>", t)
            feb = os.path.join(a.workdir, f"{model}_{a.steps}.feb")
            open(feb, "w", encoding="latin-1").write(t)
        for nt in [int(x) for x in a.threads.split(",")]:
            row = dict(model=model, threads=nt, steps=a.steps or None)
            for name, exe, solver in (("native", a.native, "pardiso"), ("rust", a.rust, "rust")):
                r = F.run(a.workdir, feb, f"{name}_t{nt}", os.path.abspath(exe), solver, nt, 0)
                row[name] = {k: r[k] for k in ("ok", "wall", "total", "linsolve", "iters", "refs",
                                               "stiff", "resid", "update", "reform", "qoi")}
            if row["native"]["qoi"] and row["rust"]["qoi"]:
                row["qoi_diff"] = Q.compare(row["native"]["qoi"], row["rust"]["qoi"])
            res = [x for x in res if not (x["model"] == model and x["threads"] == nt and x.get("steps") == row["steps"])]
            res.append(row)
            json.dump(res, open(a.out, "w"), indent=1)
            n, r = row["native"], row["rust"]
            print(f"{model} T={nt}: native {n['total']} s, rust {r['total']} s, "
                  f"x{(n['total'] or 0) / (r['total'] or 1):.2f}, iters {n['iters']}/{r['iters']}, "
                  f"qoi {row.get('qoi_diff')}", flush=True)
    lines = ["| model | threads | native total (lin.) s | Rust total (lin.) s | speed-up | iterations | QoI Δ |",
             "|---|---:|---:|---:|---:|---|---|"]
    for x in sorted(res, key=lambda x: (x["model"], x["threads"])):
        n, r = x["native"], x["rust"]
        sp = (n["total"] / r["total"]) if n["total"] and r["total"] else float("nan")
        lines.append(f"| {x['model']}{' (' + str(x['steps']) + ' steps)' if x.get('steps') else ''} | {x['threads']} | "
                     f"{n['total']} ({n['linsolve']}) | {r['total']} ({r['linsolve']}) | {sp:.2f}× | "
                     f"{n['iters']} / {r['iters']} | {x.get('qoi_diff')} |")
    open(a.md, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
