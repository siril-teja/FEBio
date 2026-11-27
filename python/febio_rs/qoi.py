#!/usr/bin/env python3
"""Spine quantities of interest (QoI) from an FEBio run.

  * range of motion: rotation of every rigid body at the last logged state (log file,
    rigid_body_data qx;qy;qz;qw), in degrees
  * intradiscal pressure (IDP): per nucleus material (name contains "nucleus", any
    case), the mean hydrostatic pressure -(sxx+syy+szz)/3 over its elements at the
    last plot state, in the model's stress unit (MPa)
  * facet contact force: per contact surface, the magnitude of FEBio's total contact
    force vector ("contact force" plot variable, FEContactSurface::GetContactForce)

The plot file must contain "stress" (for IDP) and "contact force" (for facet forces).
Values come from the .xplt (float32), so equal runs give exactly equal QoIs.
"""
import math
import os
import re

from .xplt import read_xplt


def rotations_from_log(log):
    blocks = re.findall(r"Data = qx;qy;qz;qw\n(.*?)\n\n", log, re.S)
    out = {}
    if blocks:
        for line in blocks[-1].strip().splitlines():
            v = line.split()
            q = [float(x) for x in v[1:5]]
            out[v[0]] = 2.0 * math.degrees(math.atan2(math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2), q[3]))
    return out


def rotations_from_rb_file(path):
    """Rotation (deg) of each rigid body in the last record of a rigid_body_data file
    written with data="x;y;z;qx;qy;qz;qw" (lines: id x y z qx qy qz qw)."""
    out, last = {}, {}
    for line in open(path):
        v = line.split()
        if len(v) == 8 and not line.startswith("*"):
            last[v[0]] = [float(x) for x in v[4:8]]
    for k, q in last.items():
        out[k] = 2.0 * math.degrees(math.atan2(math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2), q[3]))
    return out


def material_names(feb_text):
    return {int(i): n for i, n in re.findall(r'<material id="(\d+)" name="([^"]*)"', feb_text)}


def qoi(feb_path, xplt_path, log_text=None):
    """QoIs of one run: dict(rom_deg={body: deg}, idp={material: p}, facet_force={surface: |F|})."""
    names = material_names(open(feb_path, encoding="latin-1").read())
    res = dict(rom_deg=rotations_from_log(log_text) if log_text else {}, idp={}, facet_force={})
    if not os.path.exists(xplt_path):
        return res
    p = read_xplt(xplt_path)
    last = p.last
    stress = last.get("stress", {})
    for k, dom in enumerate(p.domains, start=1):
        name = names.get(dom["part"], dom["name"])
        if "nucleus" not in name.lower() or k not in stress:
            continue
        s = stress[k]
        n = len(s) // 6
        if n:
            pm = sum(-(s[6 * e] + s[6 * e + 1] + s[6 * e + 2]) / 3.0 for e in range(n)) / n
            # several domains can share a material: average by element count
            acc = res["idp"].setdefault(name, [0.0, 0])
            acc[0] += pm * n
            acc[1] += n
    res["idp"] = {k: v[0] / v[1] for k, v in res["idp"].items()}
    cf = last.get("contact force", {})
    for surf in p.surfaces:
        v = cf.get(surf["id"])
        if v and len(v) >= 3:
            res["facet_force"][surf["name"]] = math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return res


def compare(a, b):
    """Largest difference per QoI group between two qoi() results, relative to the
    largest magnitude in that group (so a facet that carries no load cannot dominate).
    inf if a quantity exists in only one result; None if the group is empty."""
    out = {}
    for group in ("rom_deg", "idp", "facet_force"):
        keys = set(a[group]) | set(b[group])
        if not keys:
            out[group] = None
            continue
        scale = max([abs(v) for v in a[group].values()] + [abs(v) for v in b[group].values()] + [1e-30])
        worst = 0.0
        for k in keys:
            x, y = a[group].get(k), b[group].get(k)
            worst = float("inf") if (x is None or y is None) else max(worst, abs(x - y) / scale)
        out[group] = worst
    return out


# Tolerances used by the CI gate: rotations come from the log (12 digits), IDP and
# contact forces from the plot file (float32).
TOL = {"rom_deg": 1e-9, "idp": 1e-6, "facet_force": 1e-6}


if __name__ == "__main__":
    log = open(sys.argv[3], encoding="latin-1").read() if len(sys.argv) > 3 else None
    r = qoi(sys.argv[1], sys.argv[2], log)
    for g, d in r.items():
        print(g, {k: round(v, 6) for k, v in d.items()})
