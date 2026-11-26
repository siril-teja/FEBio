#!/usr/bin/env python3
"""Generate parametric FEBio 4.0 benchmark models (hex8 / hex20 block, large-deformation).

usage: gen_model.py --n NX [--elem hex8|hex20] [--solver pardiso] [--newton] [--steps S] -o out.feb
Block [0,1]x[0,1]x[0,2], bottom face fixed, top face given a combined compression + shear.
"""
import argparse
import os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, help="elements per edge in x,y (z gets 2n)")
    ap.add_argument("--elem", default="hex8", choices=["hex8", "hex20"])
    ap.add_argument("--solver", default="pardiso")
    ap.add_argument("--newton", action="store_true", help="full Newton (refactor every iteration)")
    ap.add_argument("--steps", type=int, default=5)
    ap.add_argument("--mat", default="neo", choices=["neo", "mr"])
    ap.add_argument("-o", required=True)
    a = ap.parse_args()
    nx = ny = a.n; nz = 2 * a.n
    hx, hy, hz = 1.0 / nx, 1.0 / ny, 2.0 / nz
    quad = a.elem == "hex20"
    # node grid: for hex20 use the half-step grid, keep only serendipity nodes
    r = 2 if quad else 1
    NX, NY, NZ = r * nx + 1, r * ny + 1, r * nz + 1
    nid = {}
    nodes = []
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
                if quad and ((i % 2) + (j % 2) + (k % 2)) >= 2:
                    continue  # face/body centers not used by hex20
                nid[(i, j, k)] = len(nodes) + 1
                nodes.append((i * hx / r, j * hy / r, k * hz / r))
    elems = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                I, J, K = r * i, r * j, r * k
                c = [(0,0,0),(r,0,0),(r,r,0),(0,r,0),(0,0,r),(r,0,r),(r,r,r),(0,r,r)]
                e = [nid[(I+a_, J+b_, K+c_)] for a_, b_, c_ in c]
                if quad:
                    m = [(1,0,0),(2,1,0),(1,2,0),(0,1,0),(1,0,2),(2,1,2),(1,2,2),(0,1,2),
                         (0,0,1),(2,0,1),(2,2,1),(0,2,1)]
                    e += [nid[(I+a_, J+b_, K+c_)] for a_, b_, c_ in m]
                elems.append(e)
    bot = [nid[(i, j, 0)] for j in range(NY) for i in range(NX) if (i, j, 0) in nid]
    top = [nid[(i, j, NZ - 1)] for j in range(NY) for i in range(NX) if (i, j, NZ - 1) in nid]

    if a.mat == "neo":
        mat = '<material id="1" name="mat" type="neo-Hookean"><density>1</density><E>1</E><v>0.3</v></material>'
    else:
        mat = ('<material id="1" name="mat" type="Mooney-Rivlin"><density>1</density>'
               '<c1>0.2</c1><c2>0.05</c2><k>20</k></material>')
    qn = "<max_ups>0</max_ups>" if a.newton else "<max_ups>10</max_ups>"
    L = []
    w = L.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w('<febio_spec version="4.0">')
    w('<Module type="solid"/>')
    w('<Control><analysis>STATIC</analysis>')
    w(f'<time_steps>{a.steps}</time_steps><step_size>{1.0/a.steps}</step_size>')
    w('<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness>')
    w('<max_refs>50</max_refs><dtol>0.001</dtol><etol>0.01</etol><rtol>0</rtol><lstol>0.9</lstol>')
    w(f'<qn_method type="BFGS">{qn}</qn_method>')
    w(f'<linear_solver type="{a.solver}"/></solver>')
    w('<time_stepper type="default"><max_retries>5</max_retries><opt_iter>10</opt_iter>'
      f'<dtmin>0.001</dtmin><dtmax>{1.0/a.steps}</dtmax></time_stepper>')
    w('</Control>')
    w(f'<Material>{mat}</Material>')
    w('<Mesh><Nodes name="all">')
    for i, (x, y, z) in enumerate(nodes):
        w(f'<node id="{i+1}">{x:.10g},{y:.10g},{z:.10g}</node>')
    w('</Nodes>')
    w(f'<Elements type="{a.elem}" name="part">')
    for i, e in enumerate(elems):
        w(f'<elem id="{i+1}">' + ",".join(map(str, e)) + '</elem>')
    w('</Elements>')
    w('<NodeSet name="bottom">' + ",".join(map(str, bot)) + '</NodeSet>')
    w('<NodeSet name="top">' + ",".join(map(str, top)) + '</NodeSet>')
    w('</Mesh>')
    w('<MeshDomains><SolidDomain name="part" mat="mat"/></MeshDomains>')
    w('<Boundary>')
    w('<bc name="fix" type="zero displacement" node_set="bottom"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc>')
    w('<bc name="pz" type="prescribed displacement" node_set="top"><dof>z</dof><value lc="1">-0.5</value><relative>0</relative></bc>')
    w('<bc name="px" type="prescribed displacement" node_set="top"><dof>x</dof><value lc="1">0.3</value><relative>0</relative></bc>')
    w('<bc name="py" type="zero displacement" node_set="top"><y_dof>1</y_dof></bc>')
    w('</Boundary>')
    w('<LoadData><load_controller id="1" type="loadcurve"><interpolate>LINEAR</interpolate>'
      '<points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>')
    w('<Output><plotfile type="febio"><var type="displacement"/></plotfile>'
      f'<logfile><node_data data="ux;uy;uz" delim="," file="{os.path.splitext(os.path.basename(a.o))[0]}_disp.txt"/></logfile></Output>')
    w('</febio_spec>')
    open(a.o, "w").write("\n".join(L) + "\n")
    ndof = 3 * (len(nodes) - len(bot))
    print(f"{a.o}: {len(nodes)} nodes, {len(elems)} {a.elem}, ~{ndof} dofs")

if __name__ == "__main__":
    main()
