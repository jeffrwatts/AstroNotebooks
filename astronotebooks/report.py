"""Summary tables."""

import numpy as np
import pandas as pd

PC_TO_LY = 3.2616

_GOOD, _BAD = "color: #006300; font-weight: bold", "color: #d03b3b; font-weight: bold"


def variable_star_summary(period, vsx_period, m, M, M_err, d_pc, d_err_pc, gaia,
                          M_gaia=None, M_gaia_err=None, band="V"):
    """Side-by-side table: what this notebook calculated vs. the catalogs.

    The Check column says whether the two agree within their combined
    1-sigma uncertainty (only for rows where both sides have one).
    """
    rows = []

    def add(quantity, ours, theirs, delta="", check=""):
        rows.append([quantity, ours, theirs, delta, check])

    def verdict(diff, err):
        return "agree" if abs(diff) <= err else "disagree"

    if vsx_period:
        add("Period", f"{period:.6f} d", f"{vsx_period:.6f} d (VSX)",
            f"{(period - vsx_period) * 1440:+.2f} min")
    else:
        add("Period", f"{period:.6f} d", "not in VSX")

    if gaia:
        add(f"Apparent mag m ({band})", f"{m:.3f}", f"{gaia['G_mag']:.3f} (Gaia G)",
            f"{m - gaia['G_mag']:+.3f}", "rough only: different band")
    else:
        add(f"Apparent mag m ({band})", f"{m:.3f}", "Gaia unavailable")

    if gaia and M_gaia is not None:
        add("Absolute mag M", f"{M:.3f} ± {M_err:.3f}", f"{M_gaia:.3f} ± {M_gaia_err:.3f} (Gaia)",
            f"{M - M_gaia:+.3f}", verdict(M - M_gaia, np.hypot(M_err, M_gaia_err)))
    else:
        add("Absolute mag M", f"{M:.3f} ± {M_err:.3f}", "Gaia unavailable")

    if gaia:
        g, ge = gaia["d_pc"], gaia["d_err_pc"]
        diff = d_pc - g
        check = verdict(diff, np.hypot(d_err_pc, ge))
        add("Distance (pc)", f"{d_pc:.0f} ± {d_err_pc:.0f}", f"{g:.0f} ± {ge:.0f} (Gaia)",
            f"{diff:+.0f} ({100 * diff / g:+.1f}%)", check)
        add("Distance (light years)", f"{d_pc * PC_TO_LY:.0f} ± {d_err_pc * PC_TO_LY:.0f}",
            f"{g * PC_TO_LY:.0f} ± {ge * PC_TO_LY:.0f} (Gaia)", "", check)
    else:
        add("Distance (pc)", f"{d_pc:.0f} ± {d_err_pc:.0f}", "Gaia unavailable")

    table = pd.DataFrame(rows, columns=["Quantity", "This notebook", "Catalog", "Difference", "Check"])

    def color(val):
        return {"agree": _GOOD, "disagree": _BAD}.get(val, "")

    return table.style.hide(axis="index").map(color, subset=["Check"])


def observer_fit_table(res, band, my_observer=None, max_rows=15, excluded=None):
    """Per-observer summary of how closely their points follow the model in one band.

    offset:  median (data - model) in mag; negative = brighter than the model
    scatter: standard deviation of their residuals around their own offset
    |pull|:  median of |data - model| / quoted error bar (about 1 if the error bars are honest)

    Long tables are trimmed to the `max_rows` observers with the most points
    (plus `my_observer`). Set-aside observers (`excluded`, from fit_residuals)
    are added at the bottom, labeled "(set aside)".
    """
    d = res[res["FILT"] == band]
    t = _observer_stats(d)
    if len(t) > max_rows:   # the busiest observers, plus yours wherever it ranks
        keep = t.index[:max_rows].union(t.index[t["observer"] == my_observer])
        print(f"showing the {max_rows} observers with the most {band}-band points "
              f"(of {len(t)}){', plus ' + my_observer if my_observer in set(t['observer'][max_rows:]) else ''}")
        t = t.loc[keep]
    if excluded is not None:
        x = _observer_stats(excluded[excluded["FILT"] == band])
        x["observer"] = x["observer"] + " (set aside)"
        t = pd.concat([t, x], ignore_index=True)
    t.columns = ["Observer", "Points", "First day", "Last day", "Offset (mag)",
                 "Scatter (mag)", "Quoted error (mag)", "Typical |pull|"]

    def highlight(row):
        if row["Observer"].endswith("(set aside)"):
            return ["color: #898781; font-style: italic"] * len(row)
        mine = row["Observer"] == my_observer
        return ["font-weight: bold; background-color: #fff3c4" if mine else ""] * len(row)

    return (t.style.hide(axis="index")
             .apply(highlight, axis=1)
             .format({"First day": "{:+.1f}", "Last day": "{:+.1f}", "Offset (mag)": "{:+.3f}",
                      "Scatter (mag)": "{:.3f}", "Quoted error (mag)": "{:.3f}",
                      "Typical |pull|": "{:.1f}"}, na_rep="–"))


def _observer_stats(d):
    return (d.groupby("observer")
           .agg(points=("resid", "size"),
                first_day=("phase", "min"), last_day=("phase", "max"),
                offset=("resid", "median"),
                scatter=("resid", "std"),
                typical_error=("MERR", "median"),
                abs_pull=("pull", lambda p: np.median(np.abs(p))))
           .sort_values("points", ascending=False)
           .reset_index())
