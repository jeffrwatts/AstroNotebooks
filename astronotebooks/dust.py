"""Bayestar19 3D dust map (Green et al. 2019).

Bayestar19 tells you how much dust lies between us and a given *distance*
along a line of sight, not just the total in that direction. It needs:

- the `dustmaps` package (installed automatically on macOS/Linux; it isn't
  available on Windows), and
- a one-time ~700 MB data download: `download_bayestar()`.

Without either, `bayestar_av` returns (None, None) and the notebook carries
on with no extinction correction.
"""

import contextlib
import io
import os
from pathlib import Path

import requests

DEFAULT_DATA_DIR = ".dustmaps_data"

# Converts Bayestar19's reddening units into Johnson V extinction A_V.
# Green et al. (2019) only publish Pan-STARRS/2MASS coefficients; this one
# is derived through 2MASS Ks using the same R_V = 3.3 reddening law.
BAYESTAR_V_COEFF = 2.7476

_DATAVERSE_DOI = "10.7910/DVN/2EJ9TX"   # Bayestar19 on Harvard Dataverse
_query = None


def _bayestar_file(data_dir):
    return Path(data_dir) / "bayestar" / "bayestar2019.h5"


def download_bayestar(data_dir=DEFAULT_DATA_DIR):
    """Download the Bayestar19 data file (~700 MB, one time only).

    Uses a browser-like User-Agent: Harvard Dataverse rejects the default
    Python one, which is why `dustmaps.bayestar.fetch()` can fail.
    """
    dest = _bayestar_file(data_dir)
    if dest.exists():
        print(f"Already downloaded: {dest}")
        return dest

    headers = {"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"}
    meta = requests.get(
        f"https://dataverse.harvard.edu/api/datasets/:persistentId?persistentId=doi:{_DATAVERSE_DOI}",
        headers=headers, timeout=30).json()
    file_id = next(f["dataFile"]["id"] for f in meta["data"]["latestVersion"]["files"]
                   if f["dataFile"]["filename"] == "bayestar2019.h5")

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    print(f"Downloading Bayestar19 (~700 MB) to {dest} ...")
    with requests.get(f"https://dataverse.harvard.edu/api/access/datafile/{file_id}",
                      headers=headers, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    os.replace(tmp, dest)
    print("Done.")
    return dest


def load_bayestar(data_dir=DEFAULT_DATA_DIR):
    """Load the map once per session (takes ~20-30 s).

    Returns a `dustmaps` BayestarQuery, or None if the map is unavailable.
    """
    global _query
    if _query is not None:
        return _query
    try:
        from dustmaps.config import config
        config["data_dir"] = str(data_dir)
        from dustmaps.bayestar import BayestarQuery
    except ImportError:
        print("  `dustmaps` isn't installed, so no 3D dust map is available")
        return None
    if not _bayestar_file(data_dir).exists():
        print("  Bayestar19 data isn't downloaded yet (see download_bayestar())")
        return None
    print("  loading Bayestar19 dust map (first time only, ~30 s) ...")
    with contextlib.redirect_stdout(io.StringIO()):   # hide the loader's progress log
        _query = BayestarQuery(version="bayestar2019")
    return _query


def bayestar_av(coord, distance_pc, data_dir=DEFAULT_DATA_DIR):
    """V-band extinction A_V (mag) between us and `distance_pc` along `coord`.

    Returns (A_V, A_V_err), with A_V_err the 1-sigma spread of the map's own
    estimate, or (None, None) if the map is unavailable or the star is
    outside its footprint (declination below -30 deg).
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    q = load_bayestar(data_dir)
    if q is None or coord is None:
        return None, None
    if coord.dec.deg <= -30.0:
        print(f"  Dec {coord.dec.deg:.1f} is outside Bayestar19's sky coverage (Dec > -30 only)")
        return None, None

    c = SkyCoord(ra=coord.ra, dec=coord.dec, distance=distance_pc * u.pc)
    try:
        median = float(q(c, mode="median"))
        lo, hi = q(c, mode="percentile", pct=[16, 84])
    except Exception as e:
        print(f"  Bayestar19 query failed: {e}")
        return None, None
    return median * BAYESTAR_V_COEFF, (float(hi) - float(lo)) / 2 * BAYESTAR_V_COEFF


# --------------------------------------------------------------------------
# Whole-map queries, for the dust map notebook. All return A_V in magnitudes
# and NaN where the map has no data (declination below -30 deg).
# --------------------------------------------------------------------------

def sky_grid_av(q, l_deg, b_deg, distance_pc):
    """A_V out to `distance_pc` on a Galactic (l, b) grid.

    `l_deg` and `b_deg` are 1-D; returns an array of shape (len(b), len(l)).
    """
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord

    L, B = np.meshgrid(l_deg, b_deg)
    c = SkyCoord(l=L.ravel() * u.deg, b=B.ravel() * u.deg,
                 distance=distance_pc * u.pc, frame="galactic")
    return q(c, mode="median").reshape(L.shape) * BAYESTAR_V_COEFF


def field_grid_av(q, center, half_width_deg, distance_pc, n=241):
    """A_V out to `distance_pc` on an n x n grid of RA/Dec offsets around `center`.

    Returns (offsets_deg, A_V) with A_V of shape (n, n), rows = Dec offset.
    """
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord

    off = np.linspace(-half_width_deg, half_width_deg, n)
    dRA, dDec = np.meshgrid(off, off)
    dec = center.dec.deg + dDec
    ra = center.ra.deg + dRA / np.cos(np.radians(dec))   # equal-angle steps on the sky
    c = SkyCoord(ra=ra.ravel() * u.deg, dec=np.clip(dec, -90, 90).ravel() * u.deg,
                 distance=distance_pc * u.pc)
    return off, q(c, mode="median").reshape(dRA.shape) * BAYESTAR_V_COEFF


def line_of_sight_av(q, coord, pct=(16, 50, 84)):
    """The full dust profile along one line of sight.

    Returns a dict with the map's distance nodes (pc), A_V percentiles at each
    node (shape (n_dist, len(pct))), and the distance range (pc) the map
    considers reliable for this sightline.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=coord.ra, dec=coord.dec)
    av, flags = q(c, mode="percentile", pct=list(pct), return_flags=True)
    return {
        "distance_pc": q.distances.to(u.pc).value,
        "pct": tuple(pct),
        "A_V": av * BAYESTAR_V_COEFF,
        "converged": bool(flags["converged"]),
        "reliable_pc": (10 ** (flags["min_reliable_distmod"] / 5 + 1),
                        10 ** (flags["max_reliable_distmod"] / 5 + 1)),
    }


def plane_density(q, l_deg, b_deg=0.0, b_halfwidth_deg=1.5, max_pc=3000, dr_pc=25):
    """Dust density (A_V per kpc) in a thin slab of the Galactic disk.

    For each longitude, averages the sightlines from b - b_halfwidth to
    b + b_halfwidth (a single sightline is too noisy to read), resamples the
    cumulative A_V onto even `dr_pc` steps out to `max_pc`, and takes the
    A_V gained per step. Stretches beyond a sightline's reliable distance,
    or outside the map, come back NaN.

    Returns (r_edges_pc, density) with density of shape (len(l), n_r).
    """
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord

    l_deg = np.asarray(l_deg, dtype=float)
    bs = np.arange(b_deg - b_halfwidth_deg, b_deg + b_halfwidth_deg + 1e-9, 0.5)
    L, B = np.meshgrid(l_deg, bs, indexing="ij")
    c = SkyCoord(l=L.ravel() * u.deg, b=B.ravel() * u.deg, frame="galactic")
    av, flags = q(c, mode="median", return_flags=True)
    av = av * BAYESTAR_V_COEFF

    nodes = np.log10(q.distances.to(u.pc).value)
    r = np.arange(0, max_pc + dr_pc, dr_pc, dtype=float)
    logr = np.log10(np.maximum(r, 1.0))
    cum = np.array([np.interp(logr, nodes, row, left=0.0) for row in av])   # A_V(r), 0 inside 63 pc
    r_max = 10 ** (flags["max_reliable_distmod"] / 5 + 1)
    cum[r[None, :] > r_max[:, None]] = np.nan
    cum[np.isnan(av).any(axis=1)] = np.nan

    with np.errstate(invalid="ignore"):
        density = np.diff(cum, axis=1) / (dr_pc / 1000.0)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # all-NaN columns south of Dec -30
        density = np.nanmean(density.reshape(len(l_deg), len(bs), -1), axis=1)
    return r, density
