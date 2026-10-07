# AstroNotebooks

Measure the distance to stars and galaxies using **your own photometry** (or
public data from the AAVSO), one step at a time, in Jupyter notebooks.

| Notebook | What it measures | Method | Checked against |
|---|---|---|---|
| [`variable_star_distance.ipynb`](variable_star_distance.ipynb) | Distance to a pulsating star (hundreds of parsecs) | Period–luminosity relation | Gaia parallax |
| [`sn1a_distance.ipynb`](sn1a_distance.ipynb) | Distance to a galaxy (millions of parsecs) | Type Ia supernova light curve (SALT2 + Tripp) | Published distances (NED-D) |

Both start from the same idea, the **distance modulus**: if you know how bright
something *really* is (absolute magnitude $M$) and you measure how bright it
*looks* (apparent magnitude $m$), the difference gives its distance:

$$ d = 10^{\,(m - M + 5)/5} \text{ parsecs} $$

Each notebook is a different way of working out $M$.

---

## Installation

You need **Python 3.10 or newer** ([python.org](https://www.python.org/downloads/)).
To check, run `python3 --version` (macOS/Linux) or `py --version` (Windows).

### 1. Get the code

```bash
git clone https://github.com/jeffrwatts/AstroNotebooks.git
cd AstroNotebooks
```

(No git? On GitHub, click **Code → Download ZIP** and unzip it.)

### 2. Create an environment and install

A *virtual environment* keeps this project's packages separate from everything
else on your computer.

**macOS / Linux**
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

**Windows** (Command Prompt or PowerShell)
```bat
py -m venv .venv
.venv\Scripts\activate
pip install -e .
```

`pip install -e .` reads `pyproject.toml` and installs everything: the
science packages, JupyterLab, and this project's own helper library.

> **Using conda?** `conda create -n astro python=3.11`, then `conda activate astro`,
> then `pip install -e .`
>
> **Using [uv](https://docs.astral.sh/uv/)?** `uv venv` then `uv pip install -e .`

### 3. Start Jupyter

```bash
jupyter lab
```

Open one of the notebooks and choose **Run → Run All Cells**. Next time, you only
need to activate the environment (`source .venv/bin/activate`, or
`.venv\Scripts\activate` on Windows) and run `jupyter lab`.

### 4. (Optional) Download the 3D dust map

The variable star notebook corrects for interstellar dust using the
[Bayestar19](http://argonaut.skymaps.info/) 3D dust map, which is a **one-time
~700 MB download**. In Step 4 of the notebook, uncomment this line and run it once:

```python
dust.download_bayestar(DUST_DIR)
```

Without the map, the notebook still runs but skips the dust correction (and
reports a distance about 5–15% too large). The `dustmaps` package isn't
available on Windows, so Windows users always skip this step.

**Internet required.** Both notebooks look things up online (VSX, SIMBAD, Gaia,
NED, IRSA). If a service is slow or down, the notebook prints a message. Re-run
the cell later.

---

## The variable star notebook

Some pulsating stars follow a **period–luminosity relation**: the longer the
period, the brighter the star really is. Henrietta Leavitt discovered it in 1912.

1. **Period**: find how often the star pulses (Lomb–Scargle periodogram).
2. **Absolute magnitude $M$**: from the period, using the relation for the star's class.
3. **Apparent magnitude $m$**: the average brightness over one cycle (fold the
   light curve on the period, then average in flux, not in magnitudes).
4. **Dust $A_V$**: how much interstellar dust dims the star (3D dust map).
5. **Distance**: put $m$, $M$ and $A_V$ into the distance modulus.
6. **Check**: compare with Gaia's parallax distance.

**Sample results** (from the data in `photometry/`):

| Star | Type | Data | This notebook | Gaia |
|---|---|---|---|---|
| CY Aqr | SX Phe | VPhot report, one night | 425 ± 50 pc | 413 ± 4 pc |
| YZ Boo | SX Phe / δ Scuti | VPhot report, one night | 585 ± 66 pc | 580 ± 5 pc |
| XX Cyg | SX Phe | AAVSO download, one night | 1105 ± 123 pc | 1126 ± 20 pc |
| DY Peg | SX Phe | AAVSO download, one night | 367 ± 43 pc | 400 ± 7 pc |
| δ Cep | Cepheid | AAVSO download, 2 years | 280 ± 23 pc | 281 ± 10 pc |
| η Aql | Cepheid | AAVSO download, 1 year | 286 ± 26 pc | 272 ± 13 pc |

To switch stars, change `STAR = "CY Aqr"` near the top of the notebook.

## The Type Ia supernova notebook

Type Ia supernovae (exploding white dwarfs) all peak at *nearly* the same
brightness, and the differences can be corrected: slow-fading ones are brighter,
and red ones are dimmer.

1. **Milky Way dust** toward the host galaxy (NASA/IPAC dust map).
2. **Redshift** of the host galaxy (NED).
3. **Light-curve fit**: fit the SALT2 model to every observation in every band,
   measuring the peak time, *stretch* $x_1$ and *color* $c$.
4. **Tripp estimator**: correct the standard peak brightness for stretch and color
   to get $M$.
5. **Distance**: the distance modulus again.
6. **Check**: compare with published distances to the host galaxy (NED-D).

**Sample results:**

| Supernova | Host | This notebook | NED-D (refined) | |
|---|---|---|---|---|
| SN 2011fe | M101 | 6.52 ± 0.52 Mpc | 7.10 ± 0.56 Mpc | |
| SN 2011by | NGC 3972 | 23.3 ± 1.9 Mpc | 20.8 ± 0.7 Mpc | only one published distance |
| SN 2012fr | NGC 1365 | 19.8 ± 1.6 Mpc | 18.6 ± 1.0 Mpc | |
| SN 2014J | M82 | 2.93 ± 0.24 Mpc | 3.67 ± 0.27 Mpc | fails: M82's unusual dust |
| SN 2025rbs | NGC 7331 | 14.8 ± 1.2 Mpc | 14.9 ± 0.9 Mpc | |
| SN 2026aaiv | NGC 7331 | 15.1 ± 1.2 Mpc | 14.9 ± 0.9 Mpc | |

To switch supernovae, change `SUPERNOVA = "SN 2026aaiv"` near the top of the notebook.

---

## Using your own data

Both notebooks read either kind of AAVSO file:

- **Your own report**: the AAVSO Extended Format file exported from **VPhot** or
  **ASTAP** (transformed or not).
- **An AAVSO download**: from the [AAVSO International Database](https://www.aavso.org/data-download).
  Choose CSV.

Put the file in `photometry/` and add a line to the `TARGETS` table at the top of
the notebook. For a variable star, pick the period–luminosity relation that matches
its class. For a supernova, give its host galaxy.

## What's in this repository

```
variable_star_distance.ipynb   the variable star notebook
sn1a_distance.ipynb            the supernova notebook
photometry/                    sample AAVSO reports and downloads
astronotebooks/                helper library used by the notebooks
    aavso.py                       read AAVSO files
    catalogs.py                    online lookups: VSX, SIMBAD, Gaia, NED, IRSA
    dust.py                        Bayestar19 3D dust map
    relations.py                   published calibration constants (with citations)
    supernova.py                   SALT2 fitting helpers
    plots.py, report.py            plots and summary tables
pyproject.toml                 package list for `pip install -e .`
```

The notebooks hold the science: each step's calculation is written out in the
notebook. The library handles the plumbing (file formats, web queries,
plotting), so you can read a notebook from top to bottom without wading through
it.

## References

- Leavitt & Pickering (1912), *Periods of 25 Variable Stars in the Small Magellanic Cloud*, Harvard College Observatory Circular 173
- Cohen & Sarajedini (2012), MNRAS 419, 342: SX Phe period–luminosity relation
- Benedict et al. (2007), AJ 133, 1810: Cepheid period–luminosity relation from HST parallaxes
- Green et al. (2019), ApJ 887, 93: Bayestar19 3D dust map
- Bailer-Jones et al. (2021), AJ 161, 147: Gaia distances
- Guy et al. (2007), A&A 466, 11: the SALT2 light-curve model
- Tripp (1998), A&A 331, 815: the Tripp estimator
- Brout et al. (2022), ApJ 938, 110: Pantheon+ (SALT2-B22, α, β)
- Riess et al. (2022), ApJ 934, L7: SH0ES ($M_B^0$, $H_0$)
- Schlafly & Finkbeiner (2011), ApJ 737, 103: Milky Way dust map
- Steer et al. (2017), AJ 153, 37: NED-D distances
- Amanullah et al. (2014), ApJL 788, L21: dust toward SN 2014J

Photometry courtesy of the [AAVSO International Database](https://www.aavso.org/)
and its observers.
