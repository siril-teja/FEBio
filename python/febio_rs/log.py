"""Parse an FEBio log file: run summary, timings and logged data blocks."""
import re

NORMAL = "N O R M A L   T E R M I N A T I O N"


def _int(log, label):
    m = re.search(re.escape(label) + r"[ .]*: *([0-9]+)", log)
    return int(m.group(1)) if m else None


def _sec(log, label):
    m = re.search(re.escape(label) + r"[ .]*: [0-9:]+ \(([0-9.eE+-]+) sec\)", log)
    return float(m.group(1)) if m else None


def summary(log):
    """Run summary of a log file's text: termination, counts and timings (s)."""
    return dict(
        ok=NORMAL in log,
        equations=_int(log, "Nr of equations"),
        iterations=_int(log, "Total number of equilibrium iterations"),
        reformations=_int(log, "Total number of stiffness reformations"),
        rhs_evaluations=_int(log, "Total number of right hand evaluations"),
        linear_solves=_int(log, "Total calls to linear solver"),
        time_total=_sec(log, "Total elapsed time"),
        time_linear_solver=_sec(log, "time in linear solver"),
        time_reform=_sec(log, "reforming stiffness"),
        time_stiffness=_sec(log, "evaluating stiffness"),
        time_residual=_sec(log, "evaluating residual"),
        time_update=_sec(log, "model update"),
    )


def data(log):
    """Logged data blocks (node_data, element_data, rigid_body_data, ... written to the
    log file): list of dict(time, data, rows), rows = {item id: [values]}, in log order."""
    out = []
    pat = re.compile(r"Data Record #\d+\n=+\nStep = \d+\nTime = ([0-9.eE+-]+)\nData = ([^\n]*)\n(.*?)(?:\n\n|\Z)", re.S)
    for m in pat.finditer(log):
        rows = {}
        for line in m.group(3).strip().splitlines():
            v = line.split()
            if v:
                try:
                    rows[int(v[0])] = [float(x) for x in v[1:]]
                except ValueError:
                    pass
        out.append(dict(time=float(m.group(1)), data=m.group(2).strip(), rows=rows))
    return out
