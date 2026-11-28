#!/usr/bin/env python3
"""Check the Rust FE core's elements (ROADMAP 4.3) and materials (4.4) against C++ FEBio.

Runs the first --steps steps of a Finley et al. (2018) model with febio4 and writes, with
12 significant digits (FEBio log data records), the nodal displacements and, per element,
J and the Cauchy stress averaged over the integration points. The febio_core example
`jcheck` recomputes both from the reference mesh and those displacements (stress for the
materials the core supports: neo-Hookean, orthotropic elastic) and reports the largest
relative differences.

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
import run_finley as F  # noqa: E402


def last_block(path):
    """{item id: [values]} of the last step of a FEBio data-record file."""
    rows, cur = {}, None
    for line in open(path):
        if line.startswith("*Step"):
            cur = {}
            rows = cur
        elif line.startswith("*") or not line.strip():
            continue
        else:
            v = line.split()
            rows[int(v[0])] = [float(x) for x in v[1:]]
    return rows


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
        stem = f"{model}_j{a.steps}"
        rec = (f'\t\t\t<node_data data="ux;uy;uz" file="{stem}_u.txt"></node_data>\n'
               f'\t\t\t<element_data data="J;sx;sy;sz;sxy;syz;sxz" file="{stem}_e.txt"></element_data>\n')
        t, n = re.subn(r"(<logfile[^>]*>)", r"\1\n" + rec.replace("\\", "\\\\"), t, count=1)
        assert n == 1, "no <logfile> section"
        open(os.path.join(a.workdir, stem + ".feb"), "w", encoding="latin-1").write(t)
        env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
        subprocess.run([os.path.abspath(a.febio), "-i", stem + ".feb", "-silent", "-noconfig"],
                       cwd=a.workdir, env=env, check=False, capture_output=True)
        u = last_block(os.path.join(a.workdir, stem + "_u.txt"))
        e = last_block(os.path.join(a.workdir, stem + "_e.txt"))
        dump = os.path.join(a.workdir, stem + ".jdump")
        with open(dump, "wb") as f:
            f.write(struct.pack("<q", len(u)))
            for nid in sorted(u):
                f.write(struct.pack("<q3d", nid, *u[nid]))
            f.write(struct.pack("<q", len(e)))
            for eid in sorted(e):
                f.write(struct.pack("<q7d", eid, *e[eid]))
        out = subprocess.run([exe, os.path.join(a.workdir, stem + ".feb"), dump], capture_output=True, text=True)
        print(f"{model}: {out.stdout.strip()} {out.stderr.strip()[-300:]}")
        m = re.search(r"J max rel diff ([0-9.eE+-]+)", out.stdout)
        ms = re.search(r"stress max rel diff ([0-9.eE+-]+)", out.stdout)
        worst = max(worst, float(m.group(1)) if m else float("inf"), float(ms.group(1)) if ms else float("inf"))
    print(f"worst relative difference (J and stress): {worst:.3e}")
    sys.exit(0 if worst < 1e-5 else 1)


if __name__ == "__main__":
    main()
