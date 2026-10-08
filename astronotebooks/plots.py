"""Every plot the notebooks draw."""

import numpy as np
import matplotlib.pyplot as plt

# Colors: one hue per job, with marker shapes as a second cue for color-blind readers.
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#898781"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
GOOD, BAD = "#0ca30c", "#d03b3b"

# Supernova bands, ordered blue -> red by wavelength.
BAND_STYLES = {
    "U": ("#4a3aa7", "v"),
    "B": ("#2a78d6", "o"),
    "V": ("#008300", "s"),
    "R": ("#eda100", "^"),
    "I": ("#e34948", "D"),
}


def _style(ax, xlabel, ylabel, title, invert_y=False):
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if invert_y:
        ax.invert_yaxis()   # astronomers' convention: brighter (smaller mag) at the top


# --------------------------------------------------------------------------
# Variable stars
# --------------------------------------------------------------------------

def plot_light_curve(times, mags, title):
    """Raw magnitude vs. time for one band."""
    jd0 = np.floor(times.min())
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.scatter(times - jd0, mags, s=10, color=BLUE, alpha=0.6)
    _style(ax, f"JD - {jd0:.0f}", "magnitude", title, invert_y=True)
    fig.tight_layout()
    return fig, ax


def plot_periodogram(periods, power, best_period, title, literature_period=None):
    """Lomb-Scargle power vs. period."""
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(periods, power, lw=1, color=BLUE)
    ax.axvline(best_period, color=INK, lw=1.5, label=f"recovered: {best_period:.6f} d")
    if literature_period:
        ax.axvline(literature_period, color=ORANGE, ls="--", lw=1.5,
                   label=f"VSX catalog: {literature_period:.6f} d")
    _style(ax, "period (days)", "Lomb-Scargle power", title)
    ax.legend(frameon=False)
    fig.tight_layout()
    return fig, ax


def plot_folded(phase, mags, bin_centers, bin_mags, mean_mag, title):
    """Light curve folded on the period, shown over two cycles."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for shift in (0, 1):
        first = shift == 0
        ax.scatter(phase + shift, mags, s=8, color=GRAY, alpha=0.4,
                   label="observations" if first else None)
        ax.plot(bin_centers + shift, bin_mags, "o-", color=BLUE, ms=5, lw=2,
                label="phase-bin average" if first else None)
    ax.axhline(mean_mag, ls="--", color=ORANGE, lw=1.5,
               label=f"intensity mean m = {mean_mag:.3f}")
    _style(ax, "phase (cycles)", "magnitude", title, invert_y=True)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    return fig, ax


def plot_pl_lookup(rel, period, M, M_err, title=""):
    """The period-luminosity line over its calibrated range, showing how a
    measured period is read off as an absolute magnitude: a dotted line up
    from log10(P) to the relation, then across to M."""
    p_lo, p_hi = rel["period_range"]
    p_lo, p_hi = min(p_lo, period / 1.5), max(p_hi, period * 1.5)   # keep the star on the plot
    logP_grid = np.linspace(np.log10(p_lo), np.log10(p_hi), 200)
    M_line = rel["a"] + rel["b"] * logP_grid
    sig = np.sqrt(rel["a_err"] ** 2 + (logP_grid * rel["b_err"]) ** 2 + rel["sigma_intrinsic"] ** 2)
    logP = np.log10(period)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.fill_between(logP_grid, M_line - sig, M_line + sig, color=BLUE, alpha=0.15, lw=0,
                    label="± uncertainty (incl. star-to-star scatter)")
    ax.plot(logP_grid, M_line, color=BLUE, lw=2,
            label=f"{rel['label']}:  M = {rel['a']:.3f} {'−' if rel['b'] < 0 else '+'} {abs(rel['b']):.3f} log10(P)")

    # Bottom of an inverted y-axis is its largest value.
    y_bottom = M_line.max() + sig.max() + 0.3
    x_left = logP_grid[0]
    ax.plot([logP, logP], [y_bottom, M], ls=":", color=INK, lw=1.5)
    ax.plot([logP, x_left], [M, M], ls=":", color=INK, lw=1.5)
    ax.errorbar([logP], [M], yerr=[M_err], fmt="o", color=ORANGE, ms=8, capsize=4, zorder=3,
                label=f"this star: log10(P) = {logP:.3f}  →  M = {M:.2f} ± {M_err:.2f}")

    _style(ax, "log10(period / days)", "absolute magnitude M", title)
    ax.set_xlim(x_left, logP_grid[-1])
    ax.set_ylim(y_bottom, M_line.min() - sig.max() - 0.3)   # brighter (smaller M) at the top
    ax.legend(frameon=False, loc="upper left", fontsize=9)   # the line runs bottom-left to top-right
    fig.tight_layout()
    return fig, ax


def plot_pl_relation(rel, period, M_check=None, M_check_err=None, check_label="Gaia", title=""):
    """The period-luminosity line, zoomed to `period`, with its uncertainty
    bands, plus an optional independent M (e.g. from Gaia) to compare."""
    logP = np.linspace(np.log10(period) - 0.15, np.log10(period) + 0.15, 200)
    M_line = rel["a"] + rel["b"] * logP
    sig_cal = np.sqrt(rel["a_err"] ** 2 + (logP * rel["b_err"]) ** 2)
    sig_tot = np.sqrt(sig_cal ** 2 + rel["sigma_intrinsic"] ** 2)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.fill_between(logP, M_line - sig_tot, M_line + sig_tot, color=BLUE, alpha=0.12, lw=0,
                    label=f"total uncertainty (incl. {rel['sigma_intrinsic']:.2f} mag star-to-star scatter)")
    ax.fill_between(logP, M_line - sig_cal, M_line + sig_cal, color=BLUE, alpha=0.25, lw=0,
                    label="calibration uncertainty only")
    ax.plot(logP, M_line, color=BLUE, lw=2, label=rel["label"])
    ax.axvline(np.log10(period), color=GRAY, ls=":", lw=1)

    if M_check is not None:
        ax.errorbar([np.log10(period)], [M_check], yerr=[M_check_err], fmt="s", color=ORANGE,
                    ms=8, capsize=4, elinewidth=1.5,
                    label=f"{check_label}: M = {M_check:.2f} ± {M_check_err:.2f}")

    _style(ax, "log10(period / days)", "absolute magnitude M", title)
    M_mid = rel["a"] + rel["b"] * np.log10(period)
    half = 1.6 * np.interp(np.log10(period), logP, sig_tot)
    if M_check is not None:
        half = max(half, abs(M_check - M_mid) + M_check_err + 0.1)
    ax.set_xlim(logP[0], logP[-1])
    ax.set_ylim(M_mid + half, M_mid - half)   # brighter at the top
    ax.legend(frameon=False, loc="lower left", fontsize=9)
    fig.tight_layout()
    return fig, ax


# --------------------------------------------------------------------------
# Supernovae
# --------------------------------------------------------------------------

def plot_multiband(df, title, excluded=None, t0_guess=None, phase_window=None):
    """Raw light curve in every band SALT2 uses (U, B, V, R, I).

    Points in `excluded` (a DataFrame of set-aside observations) are drawn
    hollow. If `t0_guess` is given, the dashed line marks it and the gray
    shading marks data outside `phase_window` (dropped from the fit).
    """
    bands = list(BAND_STYLES)
    use = df[df["FILT"].isin(bands)].dropna(subset=["DATE", "MAG", "MERR"])
    out = (excluded[excluded["FILT"].isin(bands)].dropna(subset=["DATE", "MAG", "MERR"])
           if excluded is not None else use.iloc[0:0])
    jd0 = np.floor(min(use["DATE"].min(), out["DATE"].min() if len(out) else np.inf))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for band, (color, marker) in BAND_STYLES.items():
        b_use, b_out = use[use["FILT"] == band], out[out["FILT"] == band]
        if len(b_use):
            ax.scatter(b_use["DATE"] - jd0, b_use["MAG"], s=24, marker=marker, color=color,
                       alpha=0.8, edgecolors="white", linewidths=0.4, label=f"{band} ({len(b_use)})")
        if len(b_out):
            ax.scatter(b_out["DATE"] - jd0, b_out["MAG"], s=24, marker=marker,
                       facecolors="none", edgecolors=color, linewidths=1.0)
    if len(out):
        who = ", ".join(sorted(out["observer"].unique()))
        ax.scatter([], [], s=24, marker="o", facecolors="none", edgecolors=INK_2,
                   label=f"hollow = set aside ({who})")

    if t0_guess is not None and phase_window is not None:
        x_min, x_max = ax.get_xlim()
        lo, hi = t0_guess - jd0 + phase_window[0], t0_guess - jd0 + phase_window[1]
        shade = dict(color=GRAY, alpha=0.18, lw=0, zorder=0)
        ax.axvspan(x_min, lo, label="outside fit window", **shade)
        ax.axvspan(hi, x_max, **shade)
        ax.axvline(t0_guess - jd0, color=INK_2, ls="--", lw=1, label="peak guess")
        ax.set_xlim(x_min, x_max)

    _style(ax, f"JD - {jd0:.0f}", "magnitude", title, invert_y=True)
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    fig.tight_layout()
    return fig, ax


def plot_distance_check(d_mpc, d_err_mpc, ned, title):
    """Our distance vs. NED-D's published distances ("all" and "refined").

    Gray band = range of the published values; gray diamond = their average;
    colored dot = our measurement. PASS/FAIL is written out, not just colored.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharex=True)
    for ax, key, name in [(axes[0], "all", "every published distance"),
                          (axes[1], "refined", "refined (Cepheid/TRGB/SBF/PNLF)")]:
        ref = ned[key]
        passed = abs(d_mpc - ref["d_mpc"]) <= np.hypot(d_err_mpc, ref["d_err_mpc"])
        color = GOOD if passed else BAD
        ax.axvspan(*ref["range"], ymin=0.15, ymax=0.85, color=GRAY, alpha=0.2, lw=0,
                   label="range of published values")
        ax.errorbar([ref["d_mpc"]], [0.15], xerr=[ref["d_err_mpc"]], fmt="D", color=INK_2,
                    ms=7, capsize=4, label=f"NED-D ({ref['n']} papers)")
        ax.errorbar([d_mpc], [-0.15], xerr=[d_err_mpc], fmt="o", color=color,
                    ms=9, capsize=4, elinewidth=2, label="this notebook")
        ax.set_ylim(-0.7, 0.7)
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_xlabel("distance (Mpc)")
        verdict = "PASS" if passed else "FAIL"
        ax.set_title(f"vs. {name}\n{d_mpc:.2f} ± {d_err_mpc:.2f}  vs.  "
                     f"{ref['d_mpc']:.2f} ± {ref['d_err_mpc']:.2f} Mpc   {verdict}",
                     fontsize=10, color=INK)
        ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig, axes
