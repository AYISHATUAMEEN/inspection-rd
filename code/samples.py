"""samples.py - the analysis samples, defined once and used by every script."""
import os
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); IN = os.path.join(ROOT, "data", "input")
RULE_START = pd.Timestamp("2001-01-08")      # effective date of 24 CFR 200.857 (65 FR 77230)
INDEX_END = pd.Timestamp("2005-12-31")       # leaves >= 46 months of follow-up before the history file ends (Nov 2009)
HISTORY_END = pd.Timestamp("2009-11-08")
CUTOFFS = {"80": 79.5, "90": 89.5}            # designation uses the score rounded half-up (89.5 -> 90)


def load():
    i = pd.read_csv(os.path.join(IN, "ahcp_inspections.csv.gz"), dtype={"inspection_id": str, "property_id": str}, parse_dates=["inspection_date"])
    return i.sort_values(["ahcp_property_id", "inspection_date", "ahcp_inspection_id"]).reset_index(drop=True)


def history_pairs(i):
    """Multifamily inspections from the 2011 history file, each with the property's next listed inspection.
    Outcome columns are attached here but are only read by scripts run after registration."""
    h = i[(i.program == "MF") & (i.first_vintage == 2011)].copy()
    g = h.groupby("ahcp_property_id", sort=False)
    h["next_date"] = g.inspection_date.shift(-1); h["next_score"] = g.inspection_score.shift(-1)
    h["prev_score_hist"] = g.inspection_score.shift(1); h["seq_hist"] = g.cumcount() + 1
    h["years_to_next"] = (h.next_date - h.inspection_date).dt.days / 365.25
    # fixed-horizon outcome: first listed inspection at least 21 months after the index inspection
    h = h.reset_index(drop=True)
    return h


def primary(h):
    """Index inspections for the confirmatory analysis."""
    s = h[(h.inspection_date >= RULE_START) & (h.inspection_date <= INDEX_END)].copy()
    s["has_next"] = s.next_date.notna().astype(int)
    return s
