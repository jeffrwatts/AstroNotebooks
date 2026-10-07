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


def _get_query(data_dir):
    """Load the map once per session (takes ~20-30 s)."""
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

    q = _get_query(data_dir)
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
