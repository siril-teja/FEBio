"""python -m febio_rs run MODEL.feb [--solver rust] [--threads N] [--febio PATH] [--qoi]"""
import argparse
import json
import sys

from . import run


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m febio_rs")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run a model and print its summary as JSON")
    r.add_argument("feb")
    r.add_argument("--febio")
    r.add_argument("--solver")
    r.add_argument("--threads", type=int)
    r.add_argument("--nondeterministic", action="store_true")
    r.add_argument("--qoi", action="store_true", help="also print the spine quantities of interest")
    a = ap.parse_args(argv)
    res = run(a.feb, febio=a.febio, solver=a.solver, threads=a.threads, deterministic=not a.nondeterministic)
    out = res.summary()
    if a.qoi and res.ok:
        out["qoi"] = res.qoi()
    print(json.dumps(out, indent=1))
    return 0 if res.ok else 1


if __name__ == "__main__":
    sys.exit(main())
