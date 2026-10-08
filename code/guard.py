"""guard.py - the registration gate. Confirmatory scripts refuse to read real outcomes until docs/REGISTRATION.json
records the public OSF registration. With --dry-run they run end to end on a simulated outcome instead."""
import json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def mode():
    if "--dry-run" in sys.argv: return "dry-run"
    p = os.path.join(ROOT, "docs", "REGISTRATION.json")
    if not os.path.exists(p):
        sys.exit("BLOCKED: no docs/REGISTRATION.json. Register the pre-analysis plan on OSF first (docs/OSF_GUIDE.md), "
                 "then record its URL and date there. Use --dry-run to test the code on a simulated outcome.")
    r = json.load(open(p))
    if not str(r.get("osf_url", "")).startswith("https://osf.io/") or not r.get("registered_on"):
        sys.exit("BLOCKED: docs/REGISTRATION.json must contain osf_url (https://osf.io/...) and registered_on (YYYY-MM-DD).")
    return "confirmatory"


def simulated_outcome(x, seed=20261007):
    """A smooth function of the index score plus noise, with NO discontinuity. Used only in dry runs."""
    rng = np.random.default_rng(seed)
    return np.clip(60 + 0.35 * (np.asarray(x, float) - 60) + rng.normal(0, 13.5, len(x)), 0, 100)
