"""Import-only stand-in for plotnine. rdrobust, rddensity and lpdensity import plotnine for their plotting
helpers; this project draws its figures with matplotlib and never calls those helpers. Any attempt to use a
plotnine object raises an error rather than failing silently."""
class _Unavailable:
    def __init__(self, name): self._name = name
    def __call__(self, *a, **k): raise RuntimeError(f"plotnine.{self._name} is not available (import-only stub)")
    def __getattr__(self, item): return _Unavailable(f"{self._name}.{item}")
def __getattr__(name):
    if name.startswith("__"): raise AttributeError(name)
    return _Unavailable(name)
__all__ = []
