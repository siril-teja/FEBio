#!/usr/bin/env python3
"""Make the Finley et al. (2018) febio_spec 2.5 models readable by FEBio 4.

Five rewrites, all syntactic or scoping (plus one documented assumption); no value changes:

1. Fiber angles. FEBio 2.x gave fiber-exp-pow directions as <theta>/<phi> in the
   material: spherical angles in the local material axes,
   a = (cos t sin p, sin t sin p, cos p). FEBio 4 takes the same angles through
   a fiber generator, <fiber type="angles"><theta/><phi/></fiber>, evaluated with
   the same formula (FECore FESphericalAnglesVectorGenerator).
2. Material axes. FEBio 2.x put <mat_axis> on each fiber component of a
   "solid mixture"; FEBio 4 accepts it on the mixture. In these models every
   mixture is isotropic Holmes-Mow plus fiber families that all share one
   mat_axis, so moving it to the mixture is exact (FEMaterial::GetLocalCS
   composes parent*child). The script checks that assumption and refuses
   otherwise.
3. Contact name. FEBio 2.x "sliding-tension-compression" is FESlidingElasticInterface,
   registered as "sliding-elastic" since FEBio 3. Same class, same parameters.
4. Module scope. "sliding2" (FESlidingInterface2) is registered under FEBio 4's
   "biphasic" module, which builds on "solid"; FEBio 2.x had no module scoping.
   Files that use sliding2 therefore declare <Module type="biphasic"/>. They have
   no biphasic materials, so no pressure unknowns are created; the equation count
   equals the solid-module count. The biphasic solver defaults to a non-symmetric
   matrix, whereas FEBio 2.x's solid solver (which these models were built for) used
   a symmetric one and all their contacts set symmetric_stiffness=1; the converted
   file therefore sets <symmetric_stiffness>1</symmetric_stiffness> in Control.
5. Quasi-Newton method (assumption). 59F_L15_hybrid sets <qnmethod>2</qnmethod>.
   FEBio 4 reads 2 as JFNK, which needs an iterative linear solver; FEBio 2.5/2.6, which
   produced the authors' log, only had BFGS (0) and Broyden (1). BFGS reproduces the
   authors' 100 reformations (Broyden: 112) and is FEBio 2.x's default, so 2 -> 0. Both
   converge to the same rotation (19.571 vs 19.568 deg).

usage: convert_feb25.py in.feb out.feb
"""
import re
import sys

PAIR = re.compile(r"(?P<ind>[ \t]*)<theta>(?P<t>[^<]*)</theta>\s*\n[ \t]*<phi>(?P<p>[^<]*)</phi>")
MIX = re.compile(r'(<material [^>]*type="solid mixture"[^>]*>)(.*?)(</material>)', re.S)
AXIS = re.compile(r'[ \t]*<mat_axis\b[^>]*?(?:/>|>.*?</mat_axis>)[ \t]*\n', re.S)
ISOTROPIC_MATRIX = {"Holmes-Mow", "neo-Hookean"}


def fiber_angles(text):
    n = [0]

    def rep(m):
        n[0] += 1
        i = m.group("ind")
        return (f'{i}<fiber type="angles">\n{i}\t<theta>{m.group("t")}</theta>\n'
                f'{i}\t<phi>{m.group("p")}</phi>\n{i}</fiber>')

    out = PAIR.sub(rep, text)
    if re.search(r"<theta>", PAIR.sub("", text)):
        raise SystemExit("found <theta> without a following <phi>")
    return out, n[0]


def mixture_axes(text):
    n = [0]

    def rep(m):
        head, body, tail = m.groups()
        axes = AXIS.findall(body)
        if not axes:
            return m.group(0)
        if len({re.sub(r"\s+", "", a) for a in axes}) != 1:
            raise SystemExit(f"mixture with different fiber mat_axis values: {head}")
        for t in re.findall(r'<solid type="([^"]+)"', body):
            if t not in ISOTROPIC_MATRIX and not t.startswith("fiber-"):
                raise SystemExit(f"mixture component {t!r} may be anisotropic: {head}")
        n[0] += 1
        body = AXIS.sub("", body)
        ind = re.match(r"\n([ \t]*)", body).group(1) if body.startswith("\n") else "\t\t\t"
        return head + "\n" + ind + axes[0].strip() + body + tail

    return MIX.sub(rep, text), n[0]


def contact_names(text):
    old = '<contact type="sliding-tension-compression"'
    return text.replace(old, '<contact type="sliding-elastic"'), text.count(old)


def module_scope(text):
    if '<contact type="sliding2"' not in text:
        return text, 0
    if re.search(r'type="biphasic"|<permeability', text.split("<Geometry>")[0].split("<Module")[-1].split("</Material>")[0]):
        raise SystemExit("model already has biphasic materials")
    out, n = re.subn(r'<Module type="solid"\s*/>', '<Module type="biphasic"/>', text)
    control = re.sub(r"<!--.*?-->", "", out.split("</Control>")[0], flags=re.S)
    if n and "<symmetric_stiffness>" not in control:
        out = out.replace("</Control>", "\t<symmetric_stiffness>1</symmetric_stiffness>\n\t</Control>", 1)
    return out, n


def qn_method(text):
    out, n = re.subn(r"<qnmethod>\s*2\s*</qnmethod>", "<qnmethod>0</qnmethod>", text)
    return out, n


def convert(text):
    text, na = fiber_angles(text)
    text, nm = mixture_axes(text)
    text, nc = contact_names(text)
    text, nb = module_scope(text)
    text, nq = qn_method(text)
    return text, na, nm, nc + nb + nq


if __name__ == "__main__":
    src, dst = sys.argv[1:3]
    out, na, nm, nc = convert(open(src, encoding="latin-1").read())
    open(dst, "w", encoding="latin-1").write(out)
    print(f"{dst}: {na} fiber angle pairs, {nm} mixture mat_axis moves, {nc} contact/module rewrites", file=sys.stderr)
