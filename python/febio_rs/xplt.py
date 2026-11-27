#!/usr/bin/env python3
"""Minimal reader for FEBio .xplt plot files (format 3.x, as written by FEBio 4).

Reads what the spine quantities of interest need:
  * the dictionary (variable names, types, storage formats, classes)
  * the mesh: domains (part id, name, element ids) and surfaces (name, faces)
  * the last state: time and the data of every variable

Chunk layout (FEBioPlot/PltArchive): [id: u32][size: u32][payload]; branches nest.
Variable data (PLT_STATE_VAR_DATA) holds one chunk per region (node set / domain /
surface), whose id is the region index and whose payload is float32 values.

usage (library): p = read_xplt(path); p.last['contact force']['surface'][k] -> floats
"""
import struct
import zlib

ROOT, HEADER, HDR_COMPRESSION = 0x01000000, 0x01010000, 0x01010004
DICT, DIC_ITEM = 0x01020000, 0x01020001
DIC_TYPE, DIC_FMT, DIC_NAME = 0x01020002, 0x01020003, 0x01020004
DIC_CLASS = {0x01021000: "global", 0x01023000: "node", 0x01024000: "domain", 0x01025000: "surface", 0x01026000: "edge"}
MESH = 0x01040000
DOMAIN_SECTION, DOMAIN, DOMAIN_HDR = 0x01042000, 0x01042100, 0x01042101
DOM_ELEM_TYPE, DOM_PART_ID, DOM_ELEMS, DOM_NAME = 0x01042102, 0x01042103, 0x01032104, 0x01032105
DOM_ELEM_LIST, ELEMENT = 0x01042200, 0x01042201
SURFACE_SECTION, SURFACE, SURFACE_HDR = 0x01043000, 0x01043100, 0x01043101
SURFACE_ID, SURFACE_FACES, SURFACE_NAME = 0x01043102, 0x01043103, 0x01043104
STATE, STATE_HEADER, STATE_TIME, STATE_DATA = 0x02000000, 0x02010000, 0x02010002, 0x02020000
STATE_VARIABLE, VAR_ID, VAR_DATA = 0x02020001, 0x02020002, 0x02020003
DATA_CLASS = {0x02020100: "global", 0x02020300: "node", 0x02020400: "domain", 0x02020500: "surface", 0x02020600: "edge"}
VAR_SIZE = {0: 1, 1: 3, 2: 6, 3: 3, 4: 9, 5: 3, 6: 1, 7: 1}   # PLT_FLOAT, VEC3F, MAT3FS, MAT3FD, TENS4FS, MAT3F, ARRAY, ARRAY_VEC3F


def chunks(buf, pos=0, end=None):
    end = len(buf) if end is None else end
    while pos + 8 <= end:
        cid, size = struct.unpack_from("<II", buf, pos)
        yield cid, pos + 8, pos + 8 + size
        pos += 8 + size


def find(buf, cid, a, b):
    for i, s, e in chunks(buf, a, b):
        if i == cid:
            return s, e
    return None


def string(buf, s, e):
    raw = buf[s:e]
    if len(raw) >= 4:
        n = struct.unpack_from("<I", raw, 0)[0]
        if 0 < n <= len(raw) - 4 and raw[4:4 + n].isascii():   # length-prefixed
            return raw[4:4 + n].split(b"\0")[0].decode()
    return raw.split(b"\0")[0].decode(errors="replace")


def u32(buf, s):
    return struct.unpack_from("<I", buf, s)[0]


class Plot:
    def __init__(self):
        self.vars = []        # [(name, class, type, fmt)]
        self.domains = []     # [dict(name, part, elems=[element ids])]
        self.surfaces = []    # [dict(id, name, nfaces)]
        self.states = []      # [(time, {name: {region index: [floats]}})]

    @property
    def last(self):
        return self.states[-1][1]

    def var(self, name):
        for v in self.vars:
            if v[0] == name:
                return v
        raise KeyError(name)


def read_xplt(path, last_only=True):
    """Read dictionary, mesh and states (only the last state unless last_only=False)."""
    with open(path, "rb") as f:
        data = f.read()
    if u32(data, 0) != 0x00464542:
        raise ValueError(f"{path}: not an FEBio plot file")
    p = Plot()
    compressed = False
    pos = 4
    # Header and mesh are plain chunks. With compression on, everything after them is a
    # sequence of zlib streams, one per state, each holding a PLT_STATE chunk.
    top, states_raw = [], []
    for cid, s, e in chunks(data, pos):
        if cid not in (ROOT, MESH, STATE):
            break
        top.append((cid, s, e))
        pos = e
    if pos < len(data):
        rest = data[pos:]
        while rest:
            d = zlib.decompressobj()
            try:
                out = d.decompress(rest)
            except zlib.error:
                break
            for cid, s, e in chunks(out, 0):
                if cid == STATE:
                    states_raw.append(out[s:e])
                    if last_only and len(states_raw) > 1:
                        states_raw.pop(0)
            rest = d.unused_data
    if last_only:
        last_state = max((k for k, (cid, _, _) in enumerate(top) if cid == STATE), default=None)
        top = [c for k, c in enumerate(top) if c[0] != STATE or k == last_state]
        if states_raw:
            top = [c for c in top if c[0] != STATE]
    for cid, s, e in top + [(STATE, None, raw) for raw in states_raw]:
        if cid == ROOT:
            h = find(data, HEADER, s, e)
            c = find(data, HDR_COMPRESSION, *h) if h else None
            compressed = bool(c and u32(data, c[0]))
            d = find(data, DICT, s, e)
            for ccls, cs, ce in chunks(data, *d):
                for iid, is_, ie in chunks(data, cs, ce):
                    if iid != DIC_ITEM:
                        continue
                    t = find(data, DIC_TYPE, is_, ie); f = find(data, DIC_FMT, is_, ie); n = find(data, DIC_NAME, is_, ie)
                    p.vars.append((string(data, *n), DIC_CLASS.get(ccls, "?"), u32(data, t[0]), u32(data, f[0])))
        elif cid == MESH:
            ds = find(data, DOMAIN_SECTION, s, e)
            if ds:
                for did, s2, e2 in chunks(data, *ds):
                    if did != DOMAIN:
                        continue
                    hdr = find(data, DOMAIN_HDR, s2, e2)
                    part = find(data, DOM_PART_ID, *hdr); nm = find(data, DOM_NAME, *hdr)
                    el = find(data, DOM_ELEM_LIST, s2, e2)
                    ids = [u32(data, a) for _, a, _ in chunks(data, *el)] if el else []
                    p.domains.append(dict(part=u32(data, part[0]) if part else -1,
                                          name=string(data, *nm) if nm else "", elems=ids))
            ss = find(data, SURFACE_SECTION, s, e)
            if ss:
                for sid, s2, e2 in chunks(data, *ss):
                    if sid != SURFACE:
                        continue
                    hdr = find(data, SURFACE_HDR, s2, e2)
                    i = find(data, SURFACE_ID, *hdr); nf = find(data, SURFACE_FACES, *hdr); nm = find(data, SURFACE_NAME, *hdr)
                    p.surfaces.append(dict(id=u32(data, i[0]), nfaces=u32(data, nf[0]), name=string(data, *nm) if nm else ""))
        elif cid == STATE:
            payload = e if s is None else data[s:e]   # (None, bytes) for decompressed states
            sh = find(payload, STATE_HEADER, 0, len(payload))
            tm = find(payload, STATE_TIME, *sh) if sh else None
            t = struct.unpack_from("<f", payload, tm[0])[0] if tm else float("nan")
            vals = {}
            sd = find(payload, STATE_DATA, 0, len(payload))
            if sd:
                for ccls, cs, ce in chunks(payload, *sd):
                    cls = DATA_CLASS.get(ccls)
                    vars_cls = [v for v in p.vars if v[1] == cls]
                    for vid, vs, ve in chunks(payload, cs, ce):
                        if vid != STATE_VARIABLE:
                            continue
                        idc = find(payload, VAR_ID, vs, ve); dc = find(payload, VAR_DATA, vs, ve)
                        k = u32(payload, idc[0]) - 1     # 1-based within its class
                        if not (0 <= k < len(vars_cls)) or dc is None:
                            continue
                        regions = {}
                        for rid, rs, re_ in chunks(payload, *dc):
                            regions[rid] = list(struct.unpack_from(f"<{(re_ - rs) // 4}f", payload, rs))
                        vals[vars_cls[k][0]] = regions
            p.states.append((t, vals))
    return p


if __name__ == "__main__":
    import sys
    pl = read_xplt(sys.argv[1])
    print("vars:", [(v[0], v[1]) for v in pl.vars])
    print("domains:", len(pl.domains), [d["name"] for d in pl.domains[:5]])
    print("surfaces:", [(s["id"], s["name"], s["nfaces"]) for s in pl.surfaces])
    print("states:", len(pl.states), "last t =", pl.states[-1][0] if pl.states else None)
    for name, regions in pl.last.items():
        print(name, {r: len(v) for r, v in list(regions.items())[:4]})
