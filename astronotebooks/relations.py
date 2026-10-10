"""Published calibration constants used by the notebooks.

Each entry is copied from the cited paper. Never mix constants from
different entries: each set was fit together, against one sample.
"""

# Period-luminosity relations: M_V = a + b * log10(P [days])
PL_RELATIONS = {
    "SXPhe_CohenSarajedini2012": {
        "label": "SX Phe (Cohen & Sarajedini 2012)",
        "a": -1.640, "a_err": 0.110,      # zero-point
        "b": -3.389, "b_err": 0.090,      # slope
        "sigma_intrinsic": 0.20,          # star-to-star scatter around the line (mag)
        "period_range": (0.02, 0.20),     # days the relation is calibrated over
        "reference": "Cohen & Sarajedini (2012), MNRAS 419, 342",
    },
    "Cepheid_Benedict2007": {
        "label": "Classical Cepheid (Benedict et al. 2007)",
        # Published as M_V = -4.05 - 2.43 (logP - 1), from HST parallaxes of
        # 10 Galactic Cepheids. Rewritten here around logP = 0 (a = -4.05 + 2.43).
        "a": -1.620, "a_err": 0.123,
        "b": -2.430, "b_err": 0.120,
        "sigma_intrinsic": 0.10,
        "period_range": (3.7, 35.6),
        "reference": "Benedict et al. (2007), AJ 133, 1810, Table 12",
    },
}

# SALT2 + Tripp estimator calibration for Type Ia supernovae:
#   mu = m_B + alpha * x1 - beta * c - M_B0
TRIPP = {
    "label": "SALT2-B22 + Pantheon+/SH0ES",
    "alpha": 0.148, "alpha_err": 0.004,   # Brout et al. (2022), Table 2 (BS21 baseline)
    "beta": 3.09, "beta_err": 0.04,       # Brout et al. (2022), Table 2 (BS21 baseline)
    "M_B0": -19.253, "M_B0_err": 0.027,   # Riess et al. (2022), Cepheid-anchored
    "rms": 0.171,                         # Brout et al. (2022): scatter of real SNe around the relation
    "x1_range": (-3.0, 3.0),              # Brout et al. (2022): stretch and color cuts on the sample
    "c_range": (-0.3, 0.3),               #   these constants were fit to; outside them they're extrapolation
    "H0": 73.04,                          # Riess et al. (2022): the Hubble constant this implies
    "reference": "Brout et al. (2022), ApJ 938, 110; Riess et al. (2022), ApJ 934, L7",
}
