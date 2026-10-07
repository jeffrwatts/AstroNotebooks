"""Helper library for the AstroNotebooks distance notebooks.

The notebooks hold the science (each step in building a distance); this
package holds the plumbing they share: reading AAVSO files, querying online
catalogs, the 3D dust map, published calibration constants, and plotting.

Modules:
    aavso      - read AAVSO photometry files (VPhot/ASTAP reports, AID downloads)
    catalogs   - online lookups: VSX, SIMBAD, Gaia, NED, IRSA
    dust       - Bayestar19 3D dust map (download + query)
    relations  - published period-luminosity and SN Ia calibration constants
    supernova  - helpers for the SALT2 light-curve fit
    plots      - every plot the notebooks draw
    report     - summary tables
"""

__version__ = "0.1.0"

import importlib.util as _ilu
import sys as _sys

_REQUIRED = ["numpy", "pandas", "scipy", "matplotlib", "astropy", "astroquery",
             "requests", "sncosmo", "iminuit", "jinja2"]


def check_environment():
    """Fail early, with instructions, if Jupyter is running a Python that
    doesn't have this project's packages (e.g. a system or Anaconda Jupyter
    instead of the one installed in the project's .venv)."""
    missing = [m for m in _REQUIRED if _ilu.find_spec(m) is None]
    if missing:
        raise ModuleNotFoundError(
            f"This notebook is running on {_sys.executable}, which is missing: {', '.join(missing)}.\n"
            "Jupyter was probably started from a different Python than this project's environment.\n"
            "Fix: close Jupyter, then in a terminal in the AstroNotebooks folder run\n"
            "    source .venv/bin/activate      (Windows: .venv\\Scripts\\activate)\n"
            "    jupyter lab\n"
            "If you haven't installed yet, see the Installation section of README.md."
        )
    print("environment OK")
