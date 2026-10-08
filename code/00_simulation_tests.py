#!/usr/bin/env python3
"""00_simulation_tests.py - validate rdlib on synthetic data with known answers. Exits non-zero on failure."""
import sys, numpy as np
import rdlib
rng = np.random.default_rng(1); fails = []
def check(name, ok, detail): print(("PASS " if ok else "FAIL ") + name + "  " + detail); (fails.append(name) if not ok else None)

# 1. sharp RD: coverage and unbiasedness with clustered data
est, cover = [], 0; R = 400; tau = 2.0
for r in range(R):
    n = 6000; cl = rng.integers(0, 2500, n); x = rng.uniform(70, 90, n); u = rng.normal(0, 3, 2500)[cl]
    y = 60 + 0.4 * (x - 79.5) + 0.01 * (x - 79.5) ** 2 + tau * (x >= 79.5) + u + rng.normal(0, 10, n)
    o = rdlib.local_linear_rd(x, y, 79.5, 5.0, cluster=cl); est.append(o["jump"]); cover += abs(o["jump"] - tau) <= 1.96 * o["se"]
check("sharp RD unbiased", abs(np.mean(est) - tau) < 0.15, f"mean {np.mean(est):.3f} (true {tau})")
check("sharp RD coverage", 0.92 <= cover / R <= 0.975, f"{cover / R:.3f}")
# 2. fuzzy RD
est, cover = [], 0; beta = -1.5
for r in range(R):
    n = 8000; x = rng.uniform(70, 90, n); v = rng.normal(0, 1, n)
    d = 1.2 + 0.02 * (x - 79.5) + 0.8 * (x >= 79.5) + 0.5 * v + rng.normal(0, 0.4, n)
    y = 80 + 0.3 * (x - 79.5) + beta * d + 2 * v + rng.normal(0, 8, n)
    o = rdlib.local_linear_rd(x, y, 79.5, 5.0, d=d); est.append(o["wald"]); cover += abs(o["wald"] - beta) <= 1.96 * o["wald_se"]
check("fuzzy RD unbiased", abs(np.median(est) - beta) < 0.15, f"median {np.median(est):.3f} (true {beta})")
check("fuzzy RD coverage", 0.92 <= cover / R <= 0.985, f"{cover / R:.3f}")
# 3. density test: size under no manipulation, power under manipulation
rej0 = rej1 = 0
for r in range(R):
    x = rng.normal(84, 12, 40000); x = x[(x > 40) & (x < 100)]
    rej0 += abs(rdlib.density_ratio(x, 79.5, 6, 0.5)["z"]) > 1.96
    m = (x >= 78.5) & (x < 79.5) & (rng.uniform(size=len(x)) < 0.25); x2 = x.copy(); x2[m] += 1.0
    rej1 += abs(rdlib.density_ratio(x2, 79.5, 6, 0.5)["z"]) > 1.96
check("density test size", rej0 / R <= 0.10, f"{rej0 / R:.3f}")
check("density test power", rej1 / R >= 0.90, f"{rej1 / R:.3f}")
# 4. bunching: recover a known relocation of mass from 60-70 to 59
errs = []
for r in range(100):
    s = np.clip(np.round(100 - rng.gamma(2.0, 6.0, 30000)), 0, 100); m = (s >= 60) & (s <= 70) & (rng.uniform(size=len(s)) < 0.5)
    moved = m.sum(); s[m] = 59; b = rdlib.bunching(s, 40, 85, 55, 72, reps=20); errs.append((b["excess"] - moved) / moved); last = b
check("bunching excess recovered", abs(np.mean(errs)) < 0.10, f"mean relative error {np.mean(errs):+.3f}")
check("bunching excess equals missing under pure relocation", abs(last["excess_minus_missing"]) < 4 * max(last["excess_minus_missing_se"], 1) + 0.1 * last["excess"], f"{last['excess_minus_missing']:.0f} (se {last['excess_minus_missing_se']:.0f})")
print("ALL SIMULATION TESTS PASSED" if not fails else "FAILED: " + ", ".join(fails)); sys.exit(1 if fails else 0)
