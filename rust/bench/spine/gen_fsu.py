#!/usr/bin/env python3
"""Generate an L4-L5 functional spinal unit (FSU) surrogate model for FEBio 4.

This is a *surrogate*: an idealized geometry with literature-typical material
parameters, not the published Finley et al. (2018) model. It exercises the same
model ingredients that dominate solver cost in published FEBio spine models:

  * rigid vertebrae (L5 fixed, L4 loaded by a pure moment),
  * deformable intervertebral disc: nucleus (Mooney-Rivlin) + annulus
    (uncoupled Mooney-Rivlin ground matrix with two +-30 deg fiber families),
    bonded to the vertebrae via rigid node sets,
  * left/right facet joints: cartilage pads with sliding-elastic contact,
  * seven ligament groups as tension-only springs (ALL, PLL, LF, ISL, SSL, ITL, CL),
  * optional biphasic disc (unsymmetric system, transient creep),
  * optional facet friction (--friction mu, non-symmetric) and contact algorithm (--contact).

Coordinates in mm (x lateral/right, y anterior, z superior), forces in N, moments N*mm,
stress MPa.

usage: gen_fsu.py --mode flex|ext|lb|ar [--n 24] [--nz 8] [--solver pardiso] [--biphasic] -o out.feb
"""
import argparse
import math
import os

A, B, H = 25.0, 17.0, 10.0      # disc semi-axes (lateral, AP) and height
MOMENT = 7500.0                 # 7.5 N*m
FACET_INCLINATION = 30.0        # deg; facet normal tilted downward from the transverse plane
# (rigid-body rotation dof, sign) for each load case
MODES = {"flex": ("Ru", -1.0), "ext": ("Ru", 1.0), "lb": ("Rv", 1.0), "ar": ("Rw", 1.0)}

# representative linear ligament stiffnesses (N/mm, whole ligament), of the order
# reported for lumbar ligaments (e.g. Pintar et al. 1992); split over the fibers below
LIG_K = {"ALL": 33.0, "PLL": 20.0, "LF": 27.0, "ISL": 12.0, "SSL": 24.0, "ITL": 10.0, "CL": 34.0}


class Mesh:
    def __init__(self):
        self.nodes = []

    def add(self, p):
        self.nodes.append(p)
        return len(self.nodes)  # 1-based id


def disc_mesh(m, n, nz):
    """structured hex mesh of an elliptic cylinder (square-to-disk mapping)."""
    ids = {}
    for k in range(nz + 1):
        z = H * k / nz
        for j in range(n + 1):
            v = -1 + 2 * j / n
            for i in range(n + 1):
                u = -1 + 2 * i / n
                x = u * math.sqrt(1 - v * v / 2)
                y = v * math.sqrt(1 - u * u / 2)
                ids[(i, j, k)] = m.add((A * x, B * y, z))
    nuc, ann = [], []
    for k in range(nz):
        for j in range(n):
            for i in range(n):
                e = [ids[(i, j, k)], ids[(i + 1, j, k)], ids[(i + 1, j + 1, k)], ids[(i, j + 1, k)],
                     ids[(i, j, k + 1)], ids[(i + 1, j, k + 1)], ids[(i + 1, j + 1, k + 1)], ids[(i, j + 1, k + 1)]]
                cx = sum(m.nodes[q - 1][0] for q in e) / 8
                cy = sum(m.nodes[q - 1][1] for q in e) / 8
                (nuc if (cx / A) ** 2 + (cy / B) ** 2 < 0.55 ** 2 else ann).append(e)
    top = [ids[(i, j, nz)] for j in range(n + 1) for i in range(n + 1)]
    bot = [ids[(i, j, 0)] for j in range(n + 1) for i in range(n + 1)]
    return nuc, ann, top, bot


def box(m, origin, e1, e2, e3, L, ne):
    """hex block spanned by orthonormal right-handed axes e1,e2,e3 with lengths L, ne divisions.
    returns elements and a function to fetch node ids by (i,j,k)."""
    ids = {}
    for k in range(ne[2] + 1):
        for j in range(ne[1] + 1):
            for i in range(ne[0] + 1):
                a, b, c = L[0] * i / ne[0], L[1] * j / ne[1], L[2] * k / ne[2]
                p = tuple(origin[d] + a * e1[d] + b * e2[d] + c * e3[d] for d in range(3))
                ids[(i, j, k)] = m.add(p)
    el = []
    for k in range(ne[2]):
        for j in range(ne[1]):
            for i in range(ne[0]):
                el.append([ids[(i, j, k)], ids[(i + 1, j, k)], ids[(i + 1, j + 1, k)], ids[(i, j + 1, k)],
                           ids[(i, j, k + 1)], ids[(i + 1, j, k + 1)], ids[(i + 1, j + 1, k + 1)], ids[(i, j + 1, k + 1)]])
    return el, ids


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def facet_joint(m, side):
    """two cartilage pads (L4 inferior / L5 superior articular process) separated by a gap.
    returns (l4_elems, l5_elems, l4_back, l5_back, l4_face, l5_face)"""
    t, g, w, h = 1.5, 0.4, 12.0, 14.0        # pad thickness, gap, width, height
    ne = (2, 8, 8)
    ang = math.radians(45.0)                 # facet orientation: 45 deg to the sagittal plane
    inc = math.radians(FACET_INCLINATION)    # sagittal inclination (L4 pad above/behind L5 pad)
    nrm = (side * math.cos(ang) * math.cos(inc), math.sin(ang) * math.cos(inc), -math.sin(inc))
    tan = cross((0.0, 0.0, 1.0), nrm)
    tl = math.sqrt(sum(v * v for v in tan))
    tan = tuple(v / tl for v in tan)
    ez = cross(nrm, tan)                     # (nrm, tan, ez) orthonormal, right-handed
    c = (side * 15.0, -35.0, H / 2)
    def org(n0):
        return tuple(c[d] + n0 * nrm[d] - (w / 2) * tan[d] - (h / 2) * ez[d] for d in range(3))
    l4, i4 = box(m, org(-g / 2 - t), nrm, tan, ez, (t, w, h), ne)
    l5, i5 = box(m, org(g / 2), nrm, tan, ez, (t, w, h), ne)
    nn, nt, nv = ne
    l4_back = [i4[(0, j, k)] for k in range(nv + 1) for j in range(nt + 1)]
    l5_back = [i5[(nn, j, k)] for k in range(nv + 1) for j in range(nt + 1)]
    # contact faces: L4 face at i=nn (outward +nrm), L5 face at i=0 (outward -nrm)
    f4, f5 = [], []
    for k in range(nv):
        for j in range(nt):
            f4.append([i4[(nn, j, k)], i4[(nn, j + 1, k)], i4[(nn, j + 1, k + 1)], i4[(nn, j, k + 1)]])
            f5.append([i5[(0, j, k)], i5[(0, j, k + 1)], i5[(0, j + 1, k + 1)], i5[(0, j + 1, k)]])
    return l4, l5, l4_back, l5_back, f4, f5


def anchor(m, p, s=1.0):
    """tiny cube (one hex) used as a rigid attachment point; returns (elem, node at p)"""
    el, ids = box(m, (p[0] - s / 2, p[1] - s / 2, p[2] - s / 2), (1, 0, 0), (0, 1, 0), (0, 0, 1), (s, s, s), (1, 1, 1))
    return el[0], ids[(0, 0, 0)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="flex", choices=list(MODES))
    ap.add_argument("--n", type=int, default=24, help="disc in-plane divisions")
    ap.add_argument("--nz", type=int, default=8, help="disc divisions through height")
    ap.add_argument("--solver", default="pardiso")
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--biphasic", action="store_true")
    ap.add_argument("--contact", default="sliding-elastic", choices=["sliding-elastic", "sliding2"],
                    help="facet contact algorithm (sliding2 is the one the Finley 2018 models use)")
    ap.add_argument("--friction", type=float, default=0.0,
                    help="facet friction coefficient (sliding-elastic); > 0 gives a non-symmetric system")
    ap.add_argument("--newton", action="store_true")
    ap.add_argument("-o", required=True)
    a = ap.parse_args()
    if a.friction and a.contact == "sliding2":
        ap.error("--friction needs --contact sliding-elastic")
    stem = os.path.splitext(os.path.basename(a.o))[0]

    m = Mesh()
    nuc, ann, disc_top, disc_bot = disc_mesh(m, a.n, a.nz)
    facets = [facet_joint(m, s) for s in (+1, -1)]
    # rigid vertebral bodies (coarse blocks above/below the disc)
    v5, _ = box(m, (-A, -B, -25.0 - 1.0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2 * A, 2 * B, 25.0), (2, 2, 1))
    v4, _ = box(m, (-A, -B, H + 1.0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2 * A, 2 * B, 25.0), (2, 2, 1))
    # ligament attachment points: (name, [(L5 point, L4 point), ...])
    ligs = [
        ("ALL", [((x, B + 1.0, -3.0), (x, B + 1.0, H + 3.0)) for x in (-10.0, 0.0, 10.0)]),
        ("PLL", [((x, -B - 1.0, -3.0), (x, -B - 1.0, H + 3.0)) for x in (-6.0, 0.0, 6.0)]),
        ("LF",  [((x, -30.0, -2.0), (x, -30.0, H + 6.0)) for x in (-6.0, 6.0)]),
        ("ISL", [((0.0, -45.0, -4.0), (0.0, -45.0, H + 4.0))]),
        ("SSL", [((0.0, -58.0, -4.0), (0.0, -58.0, H + 4.0))]),
        ("ITL", [((x, -25.0, 0.0), (x, -25.0, H + 8.0)) for x in (-40.0, 40.0)]),
        ("CL",  [((s * 22.0, -38.0, -3.0), (s * 22.0, -38.0, H + 13.0)) for s in (-1, 1)]
               + [((s * 9.0, -33.0, -3.0), (s * 9.0, -33.0, H + 13.0)) for s in (-1, 1)]),
    ]
    anchors5, anchors4, springs = [], [], {}
    for name, pairs in ligs:
        springs[name] = []
        for p5, p4 in pairs:
            e5, n5 = anchor(m, p5)
            e4, n4 = anchor(m, p4)
            anchors5.append(e5)
            anchors4.append(e4)
            springs[name].append((n5, n4))

    L = []
    w = L.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w('<febio_spec version="4.0">')
    # sliding2 (FESlidingInterface2) is registered under FEBio 4's biphasic module
    bp_module = a.biphasic or a.contact == "sliding2"
    w(f'<Module type="{"biphasic" if bp_module else "solid"}"/>')
    w('<Control>')
    # a biphasic-module model without biphasic materials runs steady-state with the biphasic
    # solver; that is how FEBio 4 runs the Finley 2018 sliding2 models
    w(f'<analysis>{"TRANSIENT" if a.biphasic else "STEADY-STATE" if bp_module else "STATIC"}</analysis>')
    w(f'<time_steps>{a.steps}</time_steps><step_size>{1.0 / a.steps}</step_size>')
    w(f'<solver type="{"biphasic" if bp_module else "solid"}">')
    sym = a.friction == 0.0
    if not a.biphasic:
        w(f'<symmetric_stiffness>{"symmetric" if sym else "non-symmetric"}</symmetric_stiffness>')
    w('<max_refs>25</max_refs><dtol>0.001</dtol><etol>0.01</etol><rtol>0</rtol><lstol>0.9</lstol>')
    w(f'<qn_method type="BFGS"><max_ups>{0 if a.newton else 10}</max_ups></qn_method>')
    w(f'<linear_solver type="{a.solver}"/>')
    w('</solver>')
    w('<time_stepper type="default"><max_retries>10</max_retries><opt_iter>15</opt_iter>'
      f'<dtmin>0.001</dtmin><dtmax>{1.0 / a.steps}</dtmax></time_stepper>')
    w('</Control>')

    w('<Material>')
    w('<material id="1" name="L5" type="rigid body"><density>1e-6</density><E>1000</E><v>0.3</v></material>')
    w('<material id="2" name="L4" type="rigid body"><density>1e-6</density><E>1000</E><v>0.3</v></material>')
    nuc_solid = '<c1>0.12</c1><c2>0.03</c2><k>20</k>'
    fib = ('<solid type="fiber-exp-pow-uncoupled"><ksi>{ksi}</ksi><alpha>0</alpha><beta>2.5</beta>'
           '<fiber type="cylindrical"><center>0,0,0</center><axis>0,0,1</axis><vector>0,{c:.6f},{s:.6f}</vector></fiber></solid>')
    c30, s30 = math.cos(math.radians(30)), math.sin(math.radians(30))
    ann_mat = ('<solid type="Mooney-Rivlin"><c1>0.18</c1><c2>0.045</c2></solid>'
               + fib.format(ksi=8.0, c=c30, s=s30) + fib.format(ksi=8.0, c=c30, s=-s30))
    if a.biphasic:
        # biphasic disc: solid phase as above, isotropic constant permeability
        w('<material id="3" name="nucleus" type="biphasic"><phi0>0.2</phi0>'
          f'<solid type="Mooney-Rivlin">{nuc_solid}</solid>'
          '<permeability type="perm-const-iso"><perm>0.0003</perm></permeability></material>')
        w('<material id="4" name="annulus" type="biphasic"><phi0>0.3</phi0>'
          f'<solid type="uncoupled solid mixture"><k>20</k>{ann_mat}</solid>'
          '<permeability type="perm-const-iso"><perm>0.0002</perm></permeability></material>')
    else:
        w(f'<material id="3" name="nucleus" type="Mooney-Rivlin"><density>1e-6</density>{nuc_solid}</material>')
        w(f'<material id="4" name="annulus" type="uncoupled solid mixture"><density>1e-6</density><k>20</k>{ann_mat}</material>')
    w('<material id="5" name="cartilage" type="neo-Hookean"><density>1e-6</density><E>10</E><v>0.4</v></material>')
    w('</Material>')

    w('<Mesh><Nodes name="all">')
    for i, (x, y, z) in enumerate(m.nodes):
        w(f'<node id="{i + 1}">{x:.8g},{y:.8g},{z:.8g}</node>')
    w('</Nodes>')
    eid = [0]
    def elems(name, lst):
        w(f'<Elements type="hex8" name="{name}">')
        for e in lst:
            eid[0] += 1
            w(f'<elem id="{eid[0]}">' + ",".join(map(str, e)) + '</elem>')
        w('</Elements>')
    elems("nucleus", nuc)
    elems("annulus", ann)
    elems("cartilage", [e for f in facets for e in f[0] + f[1]])
    elems("L5body", v5 + anchors5)
    elems("L4body", v4 + anchors4)
    def nodeset(name, ids):
        w(f'<NodeSet name="{name}">' + ",".join(map(str, sorted(set(ids)))) + '</NodeSet>')
    nodeset("disc_top", disc_top)
    nodeset("disc_bot", disc_bot)
    nodeset("facet_L4_back", [q for f in facets for q in f[2]])
    nodeset("facet_L5_back", [q for f in facets for q in f[3]])
    fid = [0]
    for side, f in zip(("R", "L"), facets):
        for nm, faces in ((f"f4{side}", f[4]), (f"f5{side}", f[5])):
            w(f'<Surface name="{nm}">')
            for q in faces:
                fid[0] += 1
                w(f'<quad4 id="{fid[0]}">' + ",".join(map(str, q)) + '</quad4>')
            w('</Surface>')
        w(f'<SurfacePair name="facet{side}"><primary>f4{side}</primary><secondary>f5{side}</secondary></SurfacePair>')
    for name, pairs in springs.items():
        w(f'<DiscreteSet name="{name}">' + "".join(f'<delem>{p},{q}</delem>' for p, q in pairs) + '</DiscreteSet>')
    w('</Mesh>')

    w('<MeshDomains>')
    for part, mat in (("nucleus", "nucleus"), ("annulus", "annulus"), ("cartilage", "cartilage"),
                      ("L5body", "L5"), ("L4body", "L4")):
        w(f'<SolidDomain name="{part}" mat="{mat}"/>')
    w('</MeshDomains>')

    w('<Boundary>')
    w('<bc name="disc_L4" type="rigid" node_set="disc_top"><rb>L4</rb></bc>')
    w('<bc name="disc_L5" type="rigid" node_set="disc_bot"><rb>L5</rb></bc>')
    w('<bc name="facet_L4" type="rigid" node_set="facet_L4_back"><rb>L4</rb></bc>')
    w('<bc name="facet_L5" type="rigid" node_set="facet_L5_back"><rb>L5</rb></bc>')
    if a.biphasic:
        pass  # free-draining is the default for exposed biphasic surfaces (zero pressure not prescribed)
    w('</Boundary>')

    w('<Rigid>')
    w('<rigid_bc name="fixL5" type="rigid_fixed"><rb>L5</rb>'
      '<Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>1</Rz_dof><Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof></rigid_bc>')
    dof, sgn = MODES[a.mode]
    w(f'<rigid_load name="moment" type="rigid_moment"><rb>L4</rb><dof>{dof}</dof>'
      f'<value lc="1">{sgn * MOMENT}</value><relative>0</relative></rigid_load>')
    w('</Rigid>')

    w('<Contact>')
    for side in ("R", "L"):
        if a.contact == "sliding2":
            # FESlidingInterface2 (FEBioMix), frictionless; settings as in Finley et al. 2018
            w(f'<contact name="facet{side}" type="sliding2" surface_pair="facet{side}">'
              '<laugon>PENALTY</laugon><penalty>1</penalty><auto_penalty>1</auto_penalty><two_pass>1</two_pass>'
              '<symmetric_stiffness>1</symmetric_stiffness><search_tol>0.01</search_tol><search_radius>5</search_radius>'
              '<seg_up>0</seg_up></contact>')
            continue
        w(f'<contact name="facet{side}" type="sliding-elastic" surface_pair="facet{side}">'
          '<laugon>PENALTY</laugon><penalty>1</penalty><auto_penalty>1</auto_penalty><two_pass>0</two_pass>'
          f'<symmetric_stiffness>{1 if sym else 0}</symmetric_stiffness><search_tol>0.01</search_tol><search_radius>5</search_radius>'
          f'<tension>0</tension><fric_coeff>{a.friction:g}</fric_coeff></contact>')
    w('</Contact>')

    w('<Discrete>')
    for i, name in enumerate(springs):
        npair = len(springs[name])
        w(f'<discrete_material id="{i + 1}" name="{name}" type="tension-only linear spring">'
          f'<E>{LIG_K[name] / npair:.6g}</E></discrete_material>')
    for i, name in enumerate(springs):
        w(f'<discrete dmat="{i + 1}" discrete_set="{name}"/>')
    w('</Discrete>')

    w('<LoadData><load_controller id="1" type="loadcurve"><interpolate>LINEAR</interpolate>'
      '<points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>')
    w('<Output><plotfile type="febio"><var type="displacement"/><var type="stress"/></plotfile><logfile>'
      f'<rigid_body_data data="x;y;z;qx;qy;qz;qw" file="{stem}_rb.txt">2</rigid_body_data>'
      f'<node_data data="ux;uy;uz" delim="," file="{stem}_disp.txt"/>'
      '</logfile></Output>')
    w('</febio_spec>')
    open(a.o, "w").write("\n".join(L) + "\n")
    nfree = len(m.nodes) - len(set(disc_top) | set(disc_bot))
    print(f"{a.o}: {len(m.nodes)} nodes, {eid[0]} hex8 (disc {len(nuc) + len(ann)}), ~{3 * nfree} dofs")


if __name__ == "__main__":
    main()
