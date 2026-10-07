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
