"""Moved to python/febio_rs/xplt.py; kept so the benchmark scripts' imports keep working."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
from febio_rs.xplt import *  # noqa: E402,F401,F403
