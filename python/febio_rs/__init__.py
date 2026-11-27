"""febio_rs: Python driver for the FEBio-Rust febio4 executable.

    import febio_rs as fr
    r = fr.run("model.feb", solver="rust", threads=4, check=True)
    r.summary()   # iterations, reformations, timings
    r.data()      # data blocks written to the log
    r.plot()      # plot file (last state): variables, domains, surfaces
    r.qoi()       # spine quantities of interest: range of motion, IDP, facet forces
    fr.set_values("model.feb", "stiff.feb", {"Control/time_steps": 20})

Standard library only. See python/README.md.
"""
from .feb import get_value, set_values
from .log import data, summary
from .qoi import compare, qoi
from .run import Result, find_febio, run
from .xplt import read_xplt

__all__ = ["run", "Result", "find_febio", "summary", "data", "read_xplt", "qoi", "compare",
           "set_values", "get_value"]
__version__ = "0.1.0"
