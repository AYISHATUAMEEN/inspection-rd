"""Load rdrobust and rddensity with the import-only plotnine stub on the path, and record versions."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "_stubs"))
from rdrobust import rdrobust, rdbwselect  # noqa: E402
from rddensity import rddensity  # noqa: E402
from importlib.metadata import version  # noqa: E402
VERSIONS = {p: version(p) for p in ("rdrobust", "rddensity", "lpdensity", "numpy", "pandas", "scipy")}
