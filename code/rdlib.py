"""rdlib.py - small, dependency-light regression-discontinuity and bunching tools (numpy only).

local_linear_rd : sharp or fuzzy local-linear RD with a triangular kernel at a fixed bandwidth and
                  cluster-robust (CR1) standard errors. Fixed bandwidths are pre-specified in the plan;
                  data-driven bandwidths and bias-corrected inference are done with rdrobust.
density_ratio   : local-linear estimate of the log ratio of the running-variable density at the cutoff from
                  binned counts (McCrary-style), with a delta-method standard error.
bunching        : excess mass at chosen points relative to a polynomial counterfactual fitted to integer
                  score counts outside an excluded window, with a parametric-bootstrap standard error.
"""
import numpy as np


def _wls(X, y, w, cl):
    XtW = X.T * w; b = np.linalg.solve(XtW @ X, XtW @ y); e = y - X @ b
    bread = np.linalg.inv(XtW @ X)
    order = np.argsort(cl, kind="stable"); cls, start = np.unique(cl[order], return_index=True)
    S = np.add.reduceat((X * (w * e)[:, None])[order], start, axis=0)
    G, n, k = len(cls), len(y), X.shape[1]
    V = bread @ (S.T @ S) @ bread * (G / (G - 1)) * ((n - 1) / (n - k))
    return b, V, e


def local_linear_rd(x, y, c, h, cluster=None, d=None, donut=0.0, order=1):
    """Jump in E[y|x] at c. x >= c is the 'above' side. If d is given, also returns the fuzzy (Wald) estimate
    of y on d with a delta-method SE. Observations with |x - c| < donut are excluded."""
    x = np.asarray(x, float); y = np.asarray(y, float); z = x - c
    keep = (np.abs(z) <= h) & (np.abs(z) >= donut) & np.isfinite(y)
    if d is not None: d = np.asarray(d, float); keep &= np.isfinite(d)
    z, yy = z[keep], y[keep]; cl = np.arange(len(z)) if cluster is None else np.asarray(cluster)[keep]
    cl = np.unique(cl, return_inverse=True)[1]
    T = (z >= 0).astype(float); w = 1 - np.abs(z) / h
    X = np.c_[np.ones_like(z), T, z, T * z]
    if order == 2: X = np.c_[X, z ** 2, T * z ** 2]
    b, V, e = _wls(X, yy, w, cl)
    out = {"n": int(len(z)), "n_below": int((T == 0).sum()), "n_above": int((T == 1).sum()), "clusters": int(cl.max() + 1),
           "h": h, "jump": float(b[1]), "se": float(np.sqrt(V[1, 1])), "left_limit": float(b[0])}
    if d is not None:
        dd = d[keep]; bd, Vd, ed = _wls(X, dd, w, cl)
        # joint covariance of the two jumps via stacked scores
        XtW = X.T * w; bread = np.linalg.inv(XtW @ X)
        order = np.argsort(cl, kind="stable"); _, start = np.unique(cl[order], return_index=True)
        Sy = np.add.reduceat((X * (w * e)[:, None])[order], start, axis=0); Sd = np.add.reduceat((X * (w * ed)[:, None])[order], start, axis=0)
        G, n, k = len(start), len(yy), X.shape[1]; adj = (G / (G - 1)) * ((n - 1) / (n - k))
        cov = (bread @ (Sy.T @ Sd) @ bread * adj)[1, 1]
        fs, rf = bd[1], b[1]; wald = rf / fs
        var = (V[1, 1] - 2 * wald * cov + wald ** 2 * Vd[1, 1]) / fs ** 2
        out.update({"first_stage": float(fs), "first_stage_se": float(np.sqrt(Vd[1, 1])), "wald": float(wald), "wald_se": float(np.sqrt(max(var, 0)))})
    return out


def density_ratio(x, c, h, binwidth):
    """Log ratio of the density just above c to just below c from binned counts. Bins are aligned to c."""
    x = np.asarray(x, float); z = x[(np.abs(x - c) <= h)] - c
    edges = np.arange(-h, h + binwidth / 2, binwidth); cnt, _ = np.histogram(z, bins=edges); mid = (edges[:-1] + edges[1:]) / 2
    res = {}
    for side, m in (("below", mid < 0), ("above", mid > 0)):
        xm, ym = mid[m], cnt[m].astype(float); w = 1 - np.abs(xm) / h
        X = np.c_[np.ones_like(xm), xm]; XtW = X.T * w; b = np.linalg.solve(XtW @ X, XtW @ ym)
        # Poisson variance of counts
        bread = np.linalg.inv(XtW @ X); V = bread @ ((X.T * (w ** 2 * np.maximum(ym, 1))) @ X) @ bread
        res[side] = (b[0], V[0, 0])
    (fb, vb), (fa, va) = res["below"], res["above"]
    theta = np.log(fa / fb); se = np.sqrt(va / fa ** 2 + vb / fb ** 2)
    return {"log_ratio": float(theta), "se": float(se), "z": float(theta / se), "n": int(len(z)), "bins": int(len(mid))}


def bunching(scores, lo, hi, excl_lo, excl_hi, degree=5, points=(59,), deficit=(60, 70), reps=500, seed=20261007):
    """Counterfactual integer-score distribution from a polynomial fitted to counts in [lo, hi] excluding
    [excl_lo, excl_hi]. Returns excess mass at `points` and missing mass in `deficit`, in observations and as a
    share of the counterfactual, with bootstrap standard errors (Poisson resampling of counts)."""
    s = np.round(np.asarray(scores, float)).astype(int); grid = np.arange(lo, hi + 1)
    cnt = np.array([(s == g).sum() for g in grid], float); inc = (grid < excl_lo) | (grid > excl_hi)
    def fit(cc):
        zc = (grid - grid.mean()) / grid.std(); Xp = np.vander(zc, degree + 1)
        b = np.linalg.lstsq(Xp[inc], np.log(np.maximum(cc[inc], 0.5)), rcond=None)[0]
        cf = np.exp(Xp @ b)
        ex = sum(cc[grid == p][0] - cf[grid == p][0] for p in points)
        dm = (grid >= deficit[0]) & (grid <= deficit[1]); miss = float((cf[dm] - cc[dm]).sum())
        return float(ex), miss, cf
    ex, miss, cf = fit(cnt); rng = np.random.default_rng(seed)
    bs = np.array([fit(rng.poisson(cnt).astype(float))[:2] for _ in range(reps)])
    cfp = sum(cf[grid == p][0] for p in points)
    return {"excess": ex, "excess_se": float(bs[:, 0].std()), "missing": miss, "missing_se": float(bs[:, 1].std()),
            "excess_minus_missing": ex - miss, "excess_minus_missing_se": float((bs[:, 0] - bs[:, 1]).std()),
            "excess_ratio": float(ex / cfp), "n_window": int(cnt.sum()), "grid": grid.tolist(), "observed": cnt.tolist(), "counterfactual": cf.tolist()}
