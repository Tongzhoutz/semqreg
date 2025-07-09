# semqreg

**Non-Parametric Earnings-Dynamics Toolkit**  
MIT-Licensed • Maintainer: Tong Zhou ( tongzhoutz [at] gmail )

---

## Why this repo matters

The Panel Study of Income Dynamics (PSID) underpins U.S. fiscal, monetary, and social-insurance models at the Federal Reserve, CBO, SSA, and other agencies.  
`semqreg` supplies the **first fully non-parametric estimator** that:

1. **Separates persistent vs transitory risk** in earnings without functional-form assumptions;  
2. **Captures heavy-tailed, ARCH-type volatility** documented in PSID data;  
3. Reproduces Table 2 of Zhou (2024) with < 2 min runtime on a laptop.

Accurate persistence and volatility inputs flow directly into **FRB/US stress tests**, **CBO dynamic tax scoring**, and Social-Security actuarial balance calculations—central arguments in the accompanying National-Interest Waiver filing.

---

## Features

| Module | Purpose | Key Output |
|--------|---------|------------|
| `np_id.py` | Spectral operator identification of a hidden-Markov earnings process | Persistent shock `rhô` and transitory variance `σ̂²ₜ` |
| `psid_loader.py` | One-click PSID 1968-2017 extract (792 balanced households) | Clean long panel `earnings.csv` |
| `forecast_tools.py` | 5-year Gini-path & mobility-elasticity projection | `gini_proj.csv`, `mobility_proj.csv` |
| `demo_notebook.ipynb` | Replicates Figures 5–9 in Zhou (2024) | Plots + LaTeX tables |

---

## Quick start

```bash
# clone & enter
git clone https://github.com/Tongzhoutz/semqreg.git
cd semqreg

# create environment
conda env create -f environment.yml
conda activate semqreg

# run demo replication
python np_id.py             # estimate HMM
python forecast_tools.py    # 5-year projections
jupyter notebook demo_notebook.ipynb
