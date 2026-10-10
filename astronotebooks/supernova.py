"""Helpers for fitting a SALT2 light curve with sncosmo."""

import numpy as np
import pandas as pd
from astropy.table import Table

# AID band name -> sncosmo's built-in Bessell bandpass. Bessell (1990)
# curves are a close, but not exact, stand-in for amateur Johnson/Cousins filters.
SALT2_BANDS = {"Johnson U": "bessellux", "Johnson B": "bessellb", "Johnson V": "bessellv",
               "Cousins R": "bessellr", "Cousins I": "besselli"}

# Days around the guessed peak to keep. SALT2 is trained over about -15 to +45 days.
PHASE_WINDOW = (-10.0, 50.0)

ZP, ZPSYS = 25.0, "vega"   # flux zero-point for the magnitude -> flux conversion


def build_photometry_table(df, t0_guess, phase_window=PHASE_WINDOW,
                           min_points=5, min_span_days=10.0):
    """Turn AAVSO magnitudes (an AID download) into the (time, band, flux, fluxerr) table that
    `sncosmo.fit_lc` expects, using every band SALT2 understands.

    Keeps points inside `phase_window` days of `t0_guess`, and skips any band
    with too few points (or too short a time span) to say anything about the
    light-curve shape.
    """
    rows = []
    for band, bandpass in SALT2_BANDS.items():
        d = df.loc[df["band"] == band].dropna(subset=["jd", "mag", "uncertainty"])
        phase = d["jd"] - t0_guess
        d = d.loc[(phase >= phase_window[0]) & (phase <= phase_window[1])]
        if len(d) == 0:
            continue
        span = float(d["jd"].max() - d["jd"].min())
        if len(d) < min_points or span < min_span_days:
            print(f"  skipping {band}: {len(d)} point(s) over {span:.1f} d "
                  f"(need {min_points}+ points over {min_span_days:.0f}+ d)")
            continue

        flux = 10 ** (-0.4 * (d["mag"] - ZP))
        fluxerr = flux * np.log(10) / 2.5 * d["uncertainty"]
        rows += [(t, bandpass, f, fe, ZP, ZPSYS) for t, f, fe in zip(d["jd"], flux, fluxerr)]
    return Table(rows=rows, names=["time", "band", "flux", "fluxerr", "zp", "zpsys"])


def check_fit(photdata, fitted_model, bounds, t0_tolerance_days=5.0):
    """Print a warning for the two signatures of a failed SALT2 fit:
    a parameter stuck at its bound, or a peak time with no data near it."""
    ok = True
    t_lo, t_hi = float(np.min(photdata["time"])), float(np.max(photdata["time"]))
    t0 = fitted_model.get("t0")
    if t0 < t_lo - t0_tolerance_days or t0 > t_hi + t0_tolerance_days:
        print(f"  WARNING: fitted peak t0 = {t0:.1f} is far outside the data ({t_lo:.1f} to {t_hi:.1f}). "
              "The fit probably failed; look for an observer whose data disagrees with the rest.")
        ok = False
    for name, (lo, hi) in bounds.items():
        val = fitted_model.get(name)
        pad = 0.02 * (hi - lo)
        if val <= lo + pad or val >= hi - pad:
            print(f"  WARNING: {name} = {val:.3f} is stuck at its limit {bounds[name]}, "
                  "so it isn't a real measurement.")
            ok = False
    if ok:
        print("  fit looks healthy (peak inside the data, no parameter at its limit)")
    return ok


def fit_residuals(df, fit, photdata, t0_guess, phase_window=PHASE_WINDOW):
    """How far every fitted observation sits from the SALT2 model, in magnitudes.

    Uses the same points the fit used (same bands, same time window). Returns a
    DataFrame with one row per observation:
        band, observer, jd, mag, uncertainty
        phase      days from the fitted peak t0
        model_mag  the model's magnitude at that time, in that band
        resid      mag - model_mag (negative = brighter than the model)
        pull       resid / uncertainty (the miss, in units of the point's own error bar)
    """
    fitted = set(photdata["band"])
    parts = []
    for band, bandpass in SALT2_BANDS.items():
        if bandpass not in fitted:
            continue
        d = df.loc[df["band"] == band].dropna(subset=["jd", "mag", "uncertainty"])
        phase = d["jd"] - t0_guess
        d = d.loc[(phase >= phase_window[0]) & (phase <= phase_window[1])].copy()
        d["model_mag"] = fit.bandmag(bandpass, ZPSYS, d["jd"].to_numpy())
        parts.append(d)

    res = pd.concat(parts)
    if "observer" not in res.columns:
        res["observer"] = "unknown"
    res["phase"] = res["jd"] - fit.get("t0")
    res["resid"] = res["mag"] - res["model_mag"]
    res["pull"] = res["resid"] / res["uncertainty"]
    return res[["band", "observer", "jd", "phase", "mag", "uncertainty", "model_mag", "resid", "pull"]]
