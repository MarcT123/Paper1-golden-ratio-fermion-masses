# Golden Ratio Structure in Fermion Mass Ratios

**Companion code for:**
M. Trömel, "Statistical evidence for golden ratio structure in elementary fermion
mass ratios with QCD regime dependence"

## Overview

This repository provides all code to reproduce the statistical results,
figures, and falsification tests reported in the paper.

Starting from the electron mass alone, a recursive bootstrap algorithm
reproduces 11 Standard Model particle masses using five Φ-based operators
(mean error 2.5%). The operator exponents {3, 4, 11} are Lucas numbers
(L₂, L₃, L₅). Among standard mathematical constants, only Φ achieves
full reproduction at 10% tolerance; all successful bases in a 500-base
scan are powers of Φ.

## Quick Start

```bash
# Clone and install dependencies
git clone https://github.com/MarcT123/Paper1-golden-ratio-fermion-masses.git
cd Paper1-golden-ratio-fermion-masses
pip install -r requirements.txt

# For bitwise-identical p-values as reported in the paper:
pip install -r requirements-pinned.txt

# Run all statistical tests (produces terminal output)
python scripts/master_verification_v10.py

# Run regime hypothesis analysis
python scripts/regime_hypothesis_v9.py

# Generate all paper figures
cd figures
python make_figures.py
```

## Repository Structure

```
├── README.md
├── requirements.txt
├── requirements-pinned.txt
├── LICENSE
├── paper/
│   └── paper1.tex                  # Full LaTeX source
├── scripts/
│   ├── master_verification_v10.py  # All Level A, B, D tests (FINAL)
│   └── regime_hypothesis_v9.py     # Regime hypothesis, ROC, chiral
└── figures/
    └── make_figures.py             # Generates Fig. 1–3 as PDF
```

**Note on script versions:** `master_verification_v10.py` is the final,
authoritative version. Earlier versions (v9 and below) contained a
directional error in the D1 holdout test (see Version History below)
and should not be used.

## What the Scripts Reproduce

### master_verification_v10.py

| Test | Description | Expected p-value |
|------|-------------|-----------------|
| A1 | Genesis chain (MC, 3 leptons) | p ≈ 0.003 |
| A2 | Phase uniformity (Greenwood, 6 ladders) | p ≈ 0.007 |
| A3 | Φ-universality (200k random K) | p ≈ 0.001 |
| A4 | Regime hypothesis (Fisher exact, 9 hadrons) | p = 0.008 |
| A5 | Chain connectivity (MC, 5000 pseudo-spectra) | p ≈ 0.006 |
| A6 | Convention-independent subset (7 particles) | p ≈ 0.017 |
| A7 | Base scan (500 bases, Table 5) | Φ unique @ 10% |
| B1 | Ω⁻ = s × Φ⁶ (MC, Bonferroni) | p_raw ≈ 0.009, p_Bonf ≈ 0.028 |
| B2 | Decuplet amplification ≈ Φ (Bonferroni) | p_raw ≈ 0.024, p_Bonf ≈ 0.071 |
| B3 | ΔM(u→b) ≈ τ × Φ² (4 baryon pairs) | p ≈ 0.0007 |
| Fisher | Combined Level B (raw p-values, χ², 6 dof) | p ≈ 0.00002 |
| D1 | Holdout: 167 PDG hadrons (proximity test) | p ≈ 0.055 (not significant) |
| D2 | Gauge bosons (W, Z, H; computed live) | all p > 0.5 |

**Note on D1 (holdout test):** The test computes, for 167 PDG 2024
hadron masses, the mean distance to the nearest Φ-grid point (using
all 9 elementary fermion masses as seeds). The p-value is defined as
the fraction of Monte Carlo random spectra whose mean grid distance
is **≤** the observed value (i.e., as close or closer than the real
hadrons). A small p-value would indicate that hadrons are unusually
close to the Φ-grid; the observed p ≈ 0.055 is marginal and not
significant at the 5% level. A robustness check with only 4 seeds
(e, u, μ, τ) yields p ≈ 0.91 — no proximity whatsoever.

The 167 hadron masses are embedded directly in the script (39 light
mesons, 13 strange mesons, 10 charm mesons, 6 bottom mesons,
12 charmonium, 9 bottomonium, 58 light baryons, 12 charm baryons,
8 bottom baryons). All p-values are computed live via Monte Carlo —
no hardcoded values.

### regime_hypothesis_v9.py

- Fisher exact test on 9-hadron sample with QCD dominance parameter R
- ROC analysis identifying optimal threshold R_crit = 0.5
- Chiral derivation: Σ_s − Σ_ud = Δm(MS) × (Φ−1) ≈ 56 MeV

### make_figures.py

- **Figure 1:** Predicted vs. observed masses (log-log) → `fig1_pred_vs_obs.pdf`
- **Figure 2:** Φ-universality scan over K ∈ [1.2, 2.5] → `fig2_universality.pdf`
- **Figure 3:** Regime hypothesis — perfect separation at R = 0.5 → `fig3_regime.pdf`

All figures are saved to the `figures/` directory using relative paths.

## Verification Checkpoints

After running `master_verification_v10.py`, verify these approximate values:

| Checkpoint | Expected |
|-----------|----------|
| Φ-universality p | ≈ 0.001 |
| Phase uniformity (Greenwood) p | ≈ 0.007 |
| Chain connectivity (5%) p | ≈ 0.006 |
| Regime hypothesis (Fisher) p | 0.008 (exact) |
| Fisher Combined (Level B) p | ≈ 0.00002 |
| Holdout 167 hadrons p | ≈ 0.055 (not significant) |
| Holdout 4-seed robustness p | ≈ 0.91 (no proximity) |
| Bosons (W, Z, H) p | all > 0.5 |

## Reproducibility

All random number generators are seeded (`seed=42`). The code requires
Python 3.8+ with NumPy and SciPy. No custom packages or GPU required.
Total runtime < 3 minutes on a standard laptop.

Monte Carlo p-values are reproducible within ±0.002 across platforms
with identical NumPy/SciPy versions. Minor variations may occur between
NumPy versions due to RNG implementation differences. The Fisher exact
test (A4) produces an exact p-value independent of RNG.

## PDG Data

All particle masses are from PDG 2024 (Review of Particle Physics).
Mass conventions follow Table 1 of the paper:
- Pole masses for leptons and top quark
- MS-bar at μ = 2 GeV for light quarks (u, d, s)
- MS-bar at μ = m_q for charm and bottom
- Physical masses for pion and proton

## Version History

**v10 (current, March 2026)** — Final version as submitted.
- Complete 167-hadron D1 holdout test with all masses embedded in script
- D1 p-value direction corrected: counts MC spectra with mean grid
  distance **≤** observed (proximity test). This is the correct
  one-sided test: small p = hadrons unusually close to Φ-grid.
- All Level A, B, D tests in a single script

**v9 (deprecated)** — Development version. Contains a directional error
in the D1 holdout test: the Monte Carlo comparison used `>=` instead of
`<=`, effectively computing 1 − p instead of p. This produced p ≈ 0.83
(anti-proximity) instead of the correct p ≈ 0.055 (proximity). The error
affected only the D1 test; all other tests (A1–A7, B1–B3, D2) were
unaffected. **Do not use v9 for D1 results.**

## License

MIT License. See LICENSE file.

## Contact

Marc Trömel — marc.troe@gmail.com
Independent Researcher, Nürtingen, Germany
