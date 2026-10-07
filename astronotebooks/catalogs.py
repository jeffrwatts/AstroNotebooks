"""Online catalog lookups used by the notebooks.

Every function needs an internet connection. Each one returns None (and
prints why) on failure instead of raising, so one slow service doesn't
stop the rest of a notebook.
"""

import contextlib
import io
import urllib.parse

import numpy as np
import pandas as pd
import requests

_HEADERS = {"User-Agent": "Mozilla/5.0 (AstroNotebooks)"}


# --------------------------------------------------------------------------
# Stars: VSX, SIMBAD, Gaia
# --------------------------------------------------------------------------

def vsx_period(name):
    """Literature period (days) for a variable star, from the AAVSO VSX
    catalog (via its VizieR mirror). Name matching is case-insensitive."""
    params = urllib.parse.urlencode({"-source": "B/vsx/vsx", "-out": "Name,Period"})
    url = (f"https://vizier.cds.unistra.fr/viz-bin/asu-tsv?{params}"
           f"&Name=~{urllib.parse.quote(name)}")
    try:
        r = requests.get(url, headers=_HEADERS, timeout=20)
        r.raise_for_status()
    except Exception as e:
        print(f"  VSX lookup failed for {name!r}: {e}")
        return None

    lines = [l for l in r.text.splitlines() if not l.startswith("#")]
    header = next((i for i, l in enumerate(lines) if l.split("\t")[0] == "Name"), None)
    if header is None:
        return None
    data = [l for l in lines[header + 3:] if l.strip()]   # skip header, units, dashes
    if not data:
        return None
    try:
        return float(data[0].split("\t")[1])
    except (IndexError, ValueError):
        return None


def simbad_coords(name):
    """Sky position (astropy SkyCoord) for an object name, via SIMBAD."""
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    try:
        res = Simbad.query_object(name)
    except Exception as e:
        print(f"  SIMBAD lookup failed for {name!r}: {e}")
        return None
    if res is None or len(res) == 0:
        print(f"  SIMBAD could not resolve {name!r}")
        return None
    ra_col = "ra" if "ra" in res.colnames else "RA"
    dec_col = "dec" if "dec" in res.colnames else "DEC"
    ra, dec = res[ra_col][0], res[dec_col][0]
    try:
        return SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
    except (TypeError, ValueError):
        return SkyCoord(str(ra), str(dec), unit=(u.hourangle, u.deg))


def gaia_distance(name, radius_arcsec=3.0):
    """Gaia DR3 parallax and Bailer-Jones geometric distance for a star.

    Returns a dict with d_pc / d_err_pc (Bailer-Jones; preferred over 1/parallax,
    which is biased when parallax is noisy), parallax, and Gaia G magnitude.
    """
    with contextlib.redirect_stdout(io.StringIO()):   # hide the archive's import banner
        from astroquery.gaia import Gaia

    coord = simbad_coords(name)
    if coord is None:
        return None
    adql = f"""
    SELECT g.source_id, g.parallax, g.parallax_error, g.phot_g_mean_mag,
           bj.r_med_geo, bj.r_lo_geo, bj.r_hi_geo
    FROM gaiadr3.gaia_source AS g
    LEFT JOIN external.gaiaedr3_distance AS bj ON g.source_id = bj.source_id
    WHERE 1=CONTAINS(POINT('ICRS', g.ra, g.dec),
                     CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {radius_arcsec / 3600.0}))
    ORDER BY g.phot_g_mean_mag ASC
    """
    try:
        t = Gaia.launch_job(adql).get_results()
    except Exception as e:
        print(f"  Gaia query failed for {name!r}: {e}")
        return None
    if len(t) == 0:
        print(f"  no Gaia source within {radius_arcsec}\" of {name!r}")
        return None

    row = t[0]   # brightest source in the cone: the target itself
    return {
        "source_id": int(row["source_id"]),
        "parallax_mas": float(row["parallax"]),
        "parallax_err_mas": float(row["parallax_error"]),
        "d_pc": float(row["r_med_geo"]),
        "d_err_pc": (float(row["r_hi_geo"]) - float(row["r_lo_geo"])) / 2,
        "G_mag": float(row["phot_g_mean_mag"]),
    }


# --------------------------------------------------------------------------
# Galaxies: NED redshift, NED-D distances, Milky Way dust
# --------------------------------------------------------------------------

def _ned_table(endpoint, name):
    from astropy.io.votable import parse_single_table

    url = f"https://ned.ipac.caltech.edu/NED::API/{endpoint}"
    r = requests.get(url, params={"TARGET": name}, timeout=30, headers=_HEADERS)
    r.raise_for_status()
    df = parse_single_table(io.BytesIO(r.content)).to_table().to_pandas()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].apply(lambda x: x.decode() if isinstance(x, bytes) else x)
    return df


def ned_redshift(name):
    """Spectroscopic redshift of a galaxy, from NED."""
    try:
        df = _ned_table("OverviewOfObject", name)
    except Exception as e:
        print(f"  NED lookup failed for {name!r}: {e}")
        return None
    if len(df) == 0 or "z" not in df.columns or pd.isna(df["z"].iloc[0]):
        print(f"  NED has no redshift for {name!r}")
        return None
    return float(df["z"].iloc[0])


def milky_way_ebv(name):
    """Milky Way foreground reddening E(B-V) toward an object, from the
    NASA/IPAC dust map (Schlafly & Finkbeiner 2011)."""
    from astroquery.ipac.irsa.irsa_dust import IrsaDust

    try:
        table = IrsaDust.get_extinction_table(name)
    except Exception as e:
        print(f"  IRSA dust lookup failed for {name!r}: {e}")
        return None
    rows = table[table["Filter_name"] == "CTIO B"]
    if len(rows) == 0:
        return None
    # A_B / (A_B / E(B-V)) = E(B-V), at this exact position
    return float(rows[0]["A_SandF"]) / float(rows[0]["A_over_E_B_V_SandF"])


PRIMARY_METHODS = ("Cepheids", "TRGB", "SBF", "PNLF")


def ned_distances(host, methods=PRIMARY_METHODS, sigma=3.0):
    """Published distances to a galaxy from NED-D, boiled down two ways.

    "all":     every paper on file, any method, unfiltered.
    "refined": only primary distance-ladder methods (Cepheids, TRGB, SBF,
               PNLF), outliers removed with a median/MAD sigma clip. This
               is the number a new measurement is judged against.

    Tully-Fisher is left out of "refined" because its NED-D entries mostly
    assume an old H0 (~87 km/s/Mpc) and read ~20-30% low. Supernova-based
    distances are left out because they may be measurements of the very
    supernova we are checking.

    In both, a paper that reports several distances for the same galaxy is
    averaged down to one value, so it counts as one vote.

    Returns a dict, or None if the lookup fails.
    """
    from astropy.stats import sigma_clip

    try:
        rows = _ned_table("DistancesOfObject", host)
    except Exception as e:
        print(f"  NED-D lookup failed for {host!r}: {e}")
        return None
    if len(rows) == 0:
        print(f"  NED-D has no distances for {host!r}")
        return None

    rows["method"] = rows["method"].astype(str).str.strip()
    rows["refcode"] = rows["refcode"].fillna("unknown")

    def one_per_paper(df):
        return df.groupby("refcode", as_index=False).agg(
            distance=("distance", "mean"), dist_mod_err=("dist_mod_err", "mean"))

    # --- all: every paper, any method ---
    all_papers = one_per_paper(rows)
    all_d = all_papers["distance"].astype(float).to_numpy()

    # --- refined: primary methods only, then clip outliers ---
    primary = one_per_paper(rows[rows["method"].isin(methods)])
    if len(primary) == 0:
        print(f"  NED-D has no {'/'.join(methods)} distances for {host!r}")
        return None
    d = primary["distance"].astype(float).to_numpy()
    err_mag = pd.to_numeric(primary["dist_mod_err"], errors="coerce").to_numpy()
    keep = np.ones(len(d), bool)
    if len(d) >= 3:
        keep = ~sigma_clip(d, sigma=sigma, maxiters=5, cenfunc="median", stdfunc="mad_std").mask
    d, err_mag = d[keep], err_mag[keep]

    d_mpc = float(d.mean())
    # Uncertainty: the spread between papers, combined with the typical
    # precision each paper quotes for itself (converted from mag to Mpc).
    quoted = err_mag[~np.isnan(err_mag)]
    quoted_mpc = d_mpc * np.log(10) / 5 * quoted.mean() if len(quoted) else 0.0
    spread = d.std() if len(d) > 1 else 0.0
    d_err = float(np.hypot(spread, quoted_mpc)) or np.nan

    return {
        "refined": {"d_mpc": d_mpc, "d_err_mpc": d_err,
                    "range": (float(d.min()), float(d.max())), "n": len(d),
                    "n_clipped": int((~keep).sum())},
        "all": {"d_mpc": float(all_d.mean()),
                "d_err_mpc": float(all_d.std()) if len(all_d) > 1 else 0.0,
                "range": (float(all_d.min()), float(all_d.max())), "n": len(all_d)},
        "method_counts": rows["method"].value_counts().to_dict(),
        "n_rows": len(rows),
    }
