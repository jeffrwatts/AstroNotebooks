"""Read AAVSO photometry files.

Two formats are supported, auto-detected, and normalized onto the same
columns (NAME, DATE, MAG, MERR, FILT) so the notebooks don't care which
one they were given:

- AAVSO Extended Format: the report you export from VPhot or ASTAP after
  reducing your own images (transformed or not).
- AAVSO International Database (AID) CSV download: existing observations
  of any star, from https://www.aavso.org/data-download
"""

from pathlib import Path

import pandas as pd

# AID downloads use full band names ("Johnson V"); Extended Format uses short
# codes ("V"). Map AID names onto the short codes (matched case-insensitively).
# Anything not listed passes through unchanged.
AID_BAND_CODES = {
    "johnson u": "U", "johnson b": "B", "johnson v": "V",
    "cousins r": "R", "cousins i": "I",
    "sloan u'": "SU", "sloan g'": "SG", "sloan r'": "SR",
    "sloan i'": "SI", "sloan z'": "SZ",
    "tri-color green": "TG", "tri-color blue": "TB", "tri-color red": "TR",
    "unfiltered with v zeropoint": "CV", "unfiltered with r zeropoint": "CR",
    "visual": "VISUAL",
}

# Column spellings seen across AID exports (keyed lowercase).
AID_COLUMNS = {
    "jd": "DATE",
    "mag": "MAG", "magnitude": "MAG",
    "uncertainty": "MERR",
    "band": "FILT",
    "target": "NAME", "star name": "NAME",
}


def load_aavso_report(path):
    """Load an AAVSO Extended Format report or AID CSV download.

    Returns a DataFrame with (at least) NAME, DATE (JD), MAG, MERR, FILT.
    AID 'fainter than' rows (upper limits, not measurements) are dropped.
    """
    path = Path(path)
    first_line = next((l for l in path.read_text().splitlines() if l.strip()), "")
    # Extended Format always starts with a '#KEY=VALUE' line (#TYPE=EXTENDED).
    # An AID CSV can also start with '#' (its row-counter column is named '#'),
    # so test for the '#KEY=' shape rather than just a leading '#'.
    key, sep, _ = first_line.partition("=")
    if sep and key.startswith("#") and key[1:].strip().isalpha():
        return _load_extended_format(path)
    return _load_aid_download(path)


def _load_extended_format(path):
    delim = ","
    header = None
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#DELIM="):
            delim = line.split("=", 1)[1].strip()
        elif line.startswith("#NAME" + delim):
            header = line.lstrip("#").split(delim)
        elif line.startswith("#"):
            # Metadata, or (in a TransformApplier report) the original
            # untransformed row, commented out above its replacement.
            continue
        else:
            rows.append(line.split(delim))

    if header is None:
        raise ValueError(f"Could not find the #NAME,DATE,... header line in {path}")

    df = pd.DataFrame(rows, columns=header)
    for col in ["DATE", "MAG", "MERR"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _load_aid_download(path):
    df = pd.read_csv(path, low_memory=False)
    df = df.rename(columns={c: AID_COLUMNS.get(c.strip().lower(), c) for c in df.columns})
    df["FILT"] = df["FILT"].map(lambda b: AID_BAND_CODES.get(str(b).strip().lower(), b))
    for col in ["DATE", "MAG", "MERR"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if "fainterthan" in df.columns:
        faint = df["fainterthan"].astype(str).str.strip().str.lower() == "true"
        if faint.any():
            print(f"  dropped {int(faint.sum())} 'fainter than' (non-detection) rows")
            df = df.loc[~faint]
    return df


def split_by_observer(df, exclude):
    """Split `df` into (kept, excluded) rows by AAVSO observer code.

    `exclude` is an iterable of observer codes (or a {code: reason} dict).
    """
    exclude = list(exclude or [])
    if not exclude or "observer" not in df.columns:
        return df, df.iloc[0:0]
    mask = df["observer"].isin(exclude)
    if mask.any():
        print(f"  set aside {int(mask.sum())} point(s) from excluded observer(s): {', '.join(exclude)}")
    return df.loc[~mask], df.loc[mask]


def get_band_series(df, band):
    """Time (JD) and magnitude arrays for one filter, sorted by time.

    Points with no magnitude or no reported uncertainty are dropped.
    """
    d = df.loc[df["FILT"] == band].dropna(subset=["DATE", "MAG"]).sort_values("DATE")
    no_err = d["MERR"].isna()
    if no_err.any():
        print(f"  dropped {int(no_err.sum())} {band}-band point(s) with no reported uncertainty")
        d = d.loc[~no_err]
    return d["DATE"].to_numpy(), d["MAG"].to_numpy()
