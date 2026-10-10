"""Read AAVSO photometry files.

Two formats, each loaded with its own function and kept in its own schema,
so the column names match what AAVSO documents (and what other observers
call them):

- AAVSO Extended Format: the report you export from VPhot or ASTAP after
  reducing your own images. `load_extended_report()` -> NAME, DATE, MAG,
  MERR, FILT, ... with short band codes ("V").
- AAVSO International Database (AID) CSV download: existing observations
  of any star, from https://www.aavso.org/data-download.
  `load_aid_download()` -> target, jd, mag, uncertainty, band, observer, ...
  with full band names ("Johnson V").

When one notebook needs to treat both alike, `aid_as_extended()` renames an
AID download onto the Extended Format columns and band codes.
"""

from pathlib import Path

import pandas as pd

# AID band names -> Extended Format band codes (matched case-insensitively).
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

# AID column -> the Extended Format column holding the same thing.
AID_TO_EXTENDED = {
    "target": "NAME",
    "jd": "DATE",
    "mag": "MAG",
    "uncertainty": "MERR",
    "band": "FILT",
}


def load_extended_report(path):
    """Load an AAVSO Extended Format report (VPhot / ASTAP export).

    Columns come from the report's own #NAME,DATE,... header line.
    """
    path = Path(path)
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
        raise ValueError(f"Could not find the #NAME,DATE,... header line in {path} "
                         "(is it an AID download? use load_aid_download)")

    df = pd.DataFrame(rows, columns=header)
    for col in ["DATE", "MAG", "MERR"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df.loc[df["MERR"] <= 0, "MERR"] = float("nan")   # a 0.000 uncertainty means "not reported"
    return df


def load_aid_download(path):
    """Load an AAVSO International Database CSV download.

    'Fainter than' rows (upper limits, not measurements) are dropped.
    """
    df = pd.read_csv(path, low_memory=False)
    if "jd" not in df.columns:
        raise ValueError(f"No 'jd' column in {path} "
                         "(is it an Extended Format report? use load_extended_report)")
    for col in ["jd", "mag", "uncertainty"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df.loc[df["uncertainty"] <= 0, "uncertainty"] = float("nan")   # a 0.000 uncertainty means "not reported"

    faint = df["fainterthan"].astype(str).str.strip().str.lower() == "true"
    if faint.any():
        print(f"  dropped {int(faint.sum())} 'fainter than' (non-detection) rows")
        df = df.loc[~faint]
    return df


def aid_as_extended(df):
    """An AID download renamed onto the Extended Format schema.

    jd/mag/uncertainty/band/target become DATE/MAG/MERR/FILT/NAME, and band
    names become Extended Format codes ("Johnson V" -> "V"). The values are
    otherwise untouched; other columns (observer, airmass, ...) keep their names.
    """
    out = df.rename(columns=AID_TO_EXTENDED)
    out["FILT"] = out["FILT"].map(lambda b: AID_BAND_CODES.get(str(b).strip().lower(), b))
    return out


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
    Expects Extended Format columns (use `aid_as_extended` on an AID download).
    """
    d = df.loc[df["FILT"] == band].dropna(subset=["DATE", "MAG"]).sort_values("DATE")
    no_err = d["MERR"].isna()
    if no_err.any():
        print(f"  dropped {int(no_err.sum())} {band}-band point(s) with no reported uncertainty")
        d = d.loc[~no_err]
    return d["DATE"].to_numpy(), d["MAG"].to_numpy()
