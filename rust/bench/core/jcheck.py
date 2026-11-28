#!/usr/bin/env python3
"""Check the Rust FE core's element kinematics (ROADMAP 4.3) against C++ FEBio.

Runs the first --steps steps of a Finley et al. (2018) model with febio4, writing the
"relative volume" plot variable (mean J over each element's integration points). Then
the febio_core example `jcheck` recomputes J for every element from the reference mesh
and FEBio's nodal displacements, and reports the largest difference. The plot file stores
float32, so agreement is expected to about 1e-6.

usage: jcheck.py --febio PATH/febio4 --workdir DIR [--models 49F_L15_tet10,...] [--steps 2]
"""
import argparse
import os
import re
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "finley"))
sys.path.insert(0, os.path.join(ROOT, "python"))
import run_finley as F  # noqa: E402
from febio_rs import xplt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--febio", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--models", default="59F_L45_hybrid,49F_L45_hex,49F_L15_tet10")
    ap.add_argument("--steps", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    src = F.unpack(a.workdir)
    crate = os.path.join(ROOT, "rust", "febio_core")
    subprocess.run(["cargo", "build", "--release", "--example", "jcheck"], cwd=crate, check=True)
    exe = os.path.join(crate, "target", "release", "examples", "jcheck")
    worst = 0.0
    for model in a.models.split(","):
        feb = F.prepare(src, a.workdir, model)
        t = open(feb, encoding="latin-1").read()
        t = re.sub(r"<time_steps>\d+</time_steps>", f"<time_steps>{a.steps}</time_steps>", t)
        t = t.replace('<var type="displacement"/>', '<var type="displacement"/>\n\t\t\t<var type="relative volume"/>', 1)
        stem = f"{model}_j{a.steps}"
        open(os.path.join(a.workdir, stem + ".feb"), "w", encoding="latin-1").write(t)
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
        subprocess.run([os.path.abspath(a.febio), "-i", stem + ".feb", "-silent", "-noconfig"],
                       cwd=a.workdir, env=env, check=False, capture_output=True)
        p = xplt.read_xplt(os.path.join(a.workdir, stem + ".xplt"))
        d = p.last["displacement"][0]
        rv = p.last["relative volume"]
        dump = os.path.join(a.workdir, stem + ".jdump")
        with open(dump, "wb") as f:
            f.write(struct.pack("<q", len(d) // 3))
            f.write(struct.pack(f"<{len(d)}d", *d))
            f.write(struct.pack("<q", len(p.domains)))
            for k, dom in enumerate(p.domains):
                j = rv.get(k, [])
                f.write(struct.pack("<q", len(dom["elems"])))
                f.write(struct.pack(f"<{len(dom['elems'])}q", *dom["elems"]))
                f.write(struct.pack(f"<{len(j)}d", *j))
        out = subprocess.run([exe, os.path.join(a.workdir, stem + ".feb"), dump], capture_output=True, text=True)
        print(f"{model}: {out.stdout.strip()} {out.stderr.strip()[-300:]}")
        m = re.search(r"max rel diff ([0-9.eE+-]+)", out.stdout)
        worst = max(worst, float(m.group(1)) if m else float("inf"))
    print(f"worst relative difference in J: {worst:.3e}")
    sys.exit(0 if worst < 1e-5 else 1)


if __name__ == "__main__":
    main()
