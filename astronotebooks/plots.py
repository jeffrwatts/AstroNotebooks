"""Every plot the notebooks draw."""

import numpy as np
import matplotlib.pyplot as plt

# Colors: one hue per job, with marker shapes as a second cue for color-blind readers.
BLUE, ORANGE, GRAY = "#2a78d6", "#eb6834", "#898781"
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
GOOD, WARN, BAD = "#0ca30c", "#c98500", "#d03b3b"

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
    colored dot = our measurement. The verdict ("agrees" within 1 sigma,
    "consistent" within 2, "disagrees" beyond) is written out, not just colored.
    """
    from .report import agreement

    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharex=True)
    for ax, key, name in [(axes[0], "all", "every published distance"),
                          (axes[1], "refined", "refined (Cepheid/TRGB/SBF/PNLF)")]:
        ref = ned[key]
        n_sigma, tier = agreement(d_mpc, d_err_mpc, ref["d_mpc"], ref["d_err_mpc"])
        color = {"agrees": GOOD, "consistent": WARN, "disagrees": BAD}[tier]
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
        ax.set_title(f"vs. {name}\n{d_mpc:.2f} ± {d_err_mpc:.2f}  vs.  "
                     f"{ref['d_mpc']:.2f} ± {ref['d_err_mpc']:.2f} Mpc:  "
                     f"{tier.upper()} ({n_sigma:.1f}σ)", fontsize=10, color=INK)
        ax.legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig, axes


def plot_band_fit(res, fit, band, my_observer=None, title=None):
    """One band of the SALT2 fit, in detail.

    Top: the model light curve and every observation used in the fit.
    Bottom: each point's distance from the curve (data - model, in mag).
    Points from `my_observer` are drawn larger, in black, so they stand out.
    """
    from .supernova import SALT2_BANDS, ZPSYS

    color, marker = BAND_STYLES[band]
    d = res[res["FILT"] == band]
    if len(d) == 0:
        print(f"No {band}-band points were used in the fit for this supernova.")
        return None
    mine = d[d["observer"] == my_observer] if my_observer else d.iloc[0:0]
    others = d.drop(mine.index)

    t0 = fit.get("t0")
    grid = np.linspace(max(-10, d["phase"].min() - 3), 50, 400)
    with np.errstate(all="ignore"):
        model = fit.bandmag(SALT2_BANDS[band], ZPSYS, t0 + grid)
    peak = np.nanmin(model)
    peak_day = grid[np.nanargmin(model)]

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(10, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.4], "hspace": 0.08})

    # --- light curve ---
    ax.plot(grid, model, color=color, lw=2, zorder=3, label="SALT2 best fit")
    ax.errorbar(others["phase"], others["MAG"], yerr=others["MERR"], fmt=marker, ms=5,
                color=color, alpha=0.55, elinewidth=0.8, capsize=0, zorder=2,
                label=f"other observers ({len(others)} points)")
    if len(mine):
        ax.errorbar(mine["phase"], mine["MAG"], yerr=mine["MERR"], fmt=marker, ms=9,
                    color=INK, mec="white", mew=0.8, elinewidth=1.2, capsize=3, zorder=4,
                    label=f"{my_observer} ({len(mine)} points)")
    ax.axvline(0, color=GRAY, ls=":", lw=1)
    ax.plot([peak_day], [peak], marker="*", ms=12, color=color, mec=INK, mew=0.6, zorder=5,
            ls="none", label=f"model peak: {peak:.2f} mag on day {peak_day:+.1f}")
    _style(ax, "", f"{band} magnitude", title or f"{band} band: data vs. SALT2 model", invert_y=True)
    lo = np.nanmin([d["MAG"].min(), peak]) - 0.25
    hi = np.nanmax([d["MAG"].max(), np.nanmax(model)]) + 0.25
    ax.set_ylim(hi, lo)
    ax.legend(frameon=False, loc="upper right", fontsize=9)

    # --- residuals ---
    axr.axhspan(-0.05, 0.05, color=GRAY, alpha=0.15, lw=0, label="±0.05 mag")
    axr.axhline(0, color=color, lw=1.5)
    axr.errorbar(others["phase"], others["resid"], yerr=others["MERR"], fmt=marker, ms=4,
                 color=color, alpha=0.55, elinewidth=0.8, capsize=0)
    if len(mine):
        axr.errorbar(mine["phase"], mine["resid"], yerr=mine["MERR"], fmt=marker, ms=8,
                     color=INK, mec="white", mew=0.8, elinewidth=1.2, capsize=3, zorder=4)
        for _, r in mine.iterrows():
            axr.annotate(f"{r['resid']:+.3f}", (r["phase"], r["resid"]), xytext=(6, 0),
                         textcoords="offset points", fontsize=8, color=INK, va="center")
    _style(axr, "days from B-band peak (t0)", "data − model\n(mag)", "", invert_y=True)
    lim = max(0.15, np.nanpercentile(np.abs(d["resid"]), 98) * 1.15,
              np.abs(mine["resid"]).max() * 1.25 if len(mine) else 0)
    axr.set_ylim(lim, -lim)   # inverted: brighter than the model is up, as in the top panel
    axr.text(0.005, 0.95, "brighter than model ↑", transform=axr.transAxes, fontsize=8,
             color=INK_2, va="top")
    axr.legend(frameon=False, loc="upper right", fontsize=8)
    fig.align_ylabels([ax, axr])
    return fig, (ax, axr)


# --------------------------------------------------------------------------
# Bayestar19 dust map
# --------------------------------------------------------------------------

from matplotlib.colors import LinearSegmentedColormap, LogNorm

# One hue, light -> dark: pale = little dust, dark brown = a lot. Gray = no data.
DUST_CMAP = LinearSegmentedColormap.from_list(
    "dust", ["#fdf4ee", "#f6c3a4", "#eb6834", "#b6431a", "#4f1b07"])
DUST_CMAP.set_bad("#cfcec7")
AV_NORM = LogNorm(vmin=0.02, vmax=10)   # A_V spans ~3 decades across the sky


def _wrap_l(l_deg):
    """Galactic longitude to Mollweide x (radians): l = 0 in the middle,
    increasing to the left, as on every published all-sky map."""
    return -np.radians((np.asarray(l_deg) + 180) % 360 - 180)


def _mollweide_ticks(ax):
    xt = np.radians(np.arange(-120, 121, 60))
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{(-np.degrees(x)) % 360:.0f}°" for x in xt], fontsize=7, color=INK)
    ax.tick_params(axis="y", labelsize=7, labelcolor=INK_2)
    ax.grid(color=INK_2, lw=0.3, alpha=0.4)


def plot_sky_av(l_deg, b_deg, av_maps, target=None, target_label=None):
    """All-sky A_V maps (Galactic Mollweide), one panel per distance.

    `av_maps` maps distance_pc -> A_V grid from dust.sky_grid_av. All panels
    share one color scale, so you can watch the dust build up with distance.
    """
    n = len(av_maps)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 3.6 * nrow + 0.6),
                             subplot_kw={"projection": "mollweide"}, squeeze=False)
    # pcolormesh needs monotonic x, so sort longitude into Mollweide order once
    x = _wrap_l(l_deg)
    order = np.argsort(x)
    X, Y = np.meshgrid(x[order], np.radians(b_deg))
    for ax, (dist, av) in zip(axes.flat, av_maps.items()):
        shown = np.where(np.isnan(av), np.nan, np.maximum(av, AV_NORM.vmin))   # 0 is not "no data"
        mesh = ax.pcolormesh(X, Y, shown[:, order], cmap=DUST_CMAP, norm=AV_NORM,
                             shading="nearest", rasterized=True)
        if target is not None:
            g = target.galactic
            ax.plot(_wrap_l(g.l.deg), g.b.radian, marker="o", ms=9, mfc="none",
                    mec=BLUE, mew=2, ls="none")
            ax.annotate(target_label or "target", (_wrap_l(g.l.deg), g.b.radian),
                        xytext=(7, 5), textcoords="offset points", fontsize=8, color=INK)
        _mollweide_ticks(ax)
        label = f"{dist / 1000:g} kpc" if dist >= 1000 else f"{dist:g} pc"
        ax.set_title(f"dust between us and {label}", fontsize=10, color=INK, pad=8)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    cbar = fig.colorbar(mesh, ax=axes, orientation="horizontal", fraction=0.04,
                        pad=0.06, aspect=50, extend="both")
    cbar.set_label("V-band extinction $A_V$ (mag, log scale)    ·    gray: south of Dec −30°, not mapped")
    return fig, axes


def plot_line_of_sight(los, target_label, distance_pc=None):
    """Dust along one sightline: cumulative A_V (top) and where it piles up (bottom).

    Gray shading marks distances the map considers unreliable for this
    sightline (too few stars nearer/farther to constrain the dust).
    """
    d = los["distance_pc"]
    lo, med, hi = los["A_V"][:, 0], los["A_V"][:, 1], los["A_V"][:, 2]
    step = np.diff(med)

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.6], "hspace": 0.08})
    ax.fill_between(d, lo, hi, color=BLUE, alpha=0.2, lw=0,
                    label=f"{los['pct'][0]}th–{los['pct'][2]}th percentile")
    ax.plot(d, med, color=BLUE, lw=2, label="median")
    axd.bar(np.sqrt(d[:-1] * d[1:]), step, width=np.diff(d) * 0.85, color=BLUE, lw=0)

    r_lo, r_hi = los["reliable_pc"]
    for a in (ax, axd):
        if r_lo > d[0]:
            a.axvspan(d[0], r_lo, color=GRAY, alpha=0.15, lw=0)
        if r_hi < d[-1]:
            a.axvspan(r_hi, d[-1], color=GRAY, alpha=0.15, lw=0)
    ax.axvspan(np.nan, np.nan, color=GRAY, alpha=0.15, lw=0, label="outside reliable range")

    if distance_pc is not None:
        av_at = np.interp(np.log10(distance_pc), np.log10(d), med)
        for a in (ax, axd):
            a.axvline(distance_pc, color=INK, ls="--", lw=1)
        ax.plot([distance_pc], [av_at], "o", ms=8, color=INK, mec="white", mew=1.5, zorder=4)
        ax.annotate(f"$A_V$ = {av_at:.2f} mag at {distance_pc:.0f} pc", (distance_pc, av_at),
                    xytext=(8, -14), textcoords="offset points", fontsize=9, color=INK)

    ax.set_xscale("log")
    ax.set_xlim(d[0], d[-1])
    ax.set_ylim(bottom=0)
    title = f"Dust toward {target_label}"
    if not los["converged"]:
        title += "  (map's fit did not converge here: treat with caution)"
    _style(ax, "", "$A_V$ out to distance (mag)", title)
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    _style(axd, "distance (pc, log scale)", "$A_V$ added\nper bin (mag)", "")
    fig.align_ylabels([ax, axd])
    return fig, (ax, axd)


def plot_field_av(offsets_deg, av, center_label, distance_pc):
    """A_V out to one distance, in a square patch of sky around a target."""
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    mesh = ax.pcolormesh(offsets_deg, offsets_deg, av, cmap=DUST_CMAP,
                         shading="nearest", rasterized=True)
    ax.plot(0, 0, marker="o", ms=12, mfc="none", mec=BLUE, mew=2)
    ax.set_aspect("equal")
    ax.invert_xaxis()   # east to the left, as on the sky
    ax.set(xlabel="RA offset (deg, east ←)", ylabel="Dec offset (deg)")
    ax.set_title(f"{center_label}: dust out to {distance_pc:.0f} pc", fontsize=11, color=INK)
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("$A_V$ (mag)")
    fig.tight_layout()
    return fig, ax


def plot_plane_density(l_deg, r_edges_pc, density, b_deg=0.0, target=None,
                       target_label=None, distance_pc=None):
    """Bird's-eye view of the Galactic disk: dust density around the Sun.

    The Sun is at the center, the Galactic center straight up (l = 0), and
    longitude increases counter-clockwise, as seen from the north Galactic pole.
    """
    max_pc = r_edges_pc[-1]
    pos = density[np.isfinite(density) & (density > 0)]
    vmax = np.percentile(pos, 99) if pos.size else 1

    step = np.diff(l_deg).mean()   # pcolormesh needs cell edges in longitude too
    th = np.radians(np.append(l_deg - step / 2, l_deg[-1] + step / 2))

    fig, ax = plt.subplots(figsize=(8, 8.4), subplot_kw={"projection": "polar"})
    mesh = ax.pcolormesh(th, r_edges_pc, density.T, cmap=DUST_CMAP, vmin=0, vmax=vmax,
                         rasterized=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(1)
    ax.set_thetagrids(range(0, 360, 30), [f"{a}°" for a in range(0, 360, 30)],
                      fontsize=8, color=INK_2)
    ax.set_rlim(0, max_pc)
    ax.set_rlabel_position(105)
    ax.tick_params(axis="y", labelsize=8, labelcolor=INK_2)
    ax.grid(color=INK_2, lw=0.3, alpha=0.5)
    ax.plot(0, 0, marker="*", ms=14, color="#f2c94c", mec=INK, mew=0.6)
    ax.annotate("Sun", (0, 0), xytext=(8, -12), textcoords="offset points", fontsize=8, color=INK)

    if target is not None and distance_pc is not None:
        l_t = target.galactic.l.radian
        ax.plot(l_t, min(distance_pc, max_pc), marker="o", ms=10, mfc="none", mec=BLUE, mew=2)
        ax.annotate(target_label or "target", (l_t, min(distance_pc, max_pc)),
                    xytext=(9, 4), textcoords="offset points", fontsize=9, color=INK)
        if abs(target.galactic.b.deg - b_deg) > 3:
            ax.annotate(f"({target_label or 'target'} is at b = {target.galactic.b.deg:+.0f}°: "
                        "shown at its longitude, above/below this slab)",
                        (0.5, -0.07), xycoords="axes fraction", ha="center", fontsize=8, color=INK_2)

    ax.set_title(f"Dust in the Galactic plane (b = {b_deg:g}°), seen from above\n"
                 "↑ Galactic center (l = 0°)  ·  rings: distance from the Sun (pc)  ·  "
                 "gray: unreliable or not mapped", fontsize=10, color=INK, pad=22)
    cbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", fraction=0.04, pad=0.09,
                        aspect=40, extend="max")
    cbar.set_label("dust density: $A_V$ per kpc of path (mag/kpc)")
    return fig, ax
