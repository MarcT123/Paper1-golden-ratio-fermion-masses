# Golden Ratio Structure in Fermion Mass Ratios

**Companion code for:**
M. Tromel, "A Golden Ratio Structure in the Fermion Mass Spectrum: Statistical Evidence and Regime Dependence"

## Overview

This repository provides all code to reproduce the statistical results,
figures, and falsification tests reported in the paper.

Starting from the electron mass alone, a recursive bootstrap algorithm
reproduces 11 Standard Model particle masses using five Φ-based operators
(mean error 2.5%). The operator exponents {3, 4, 11} are consecutive
Lucas numbers. A universality scan over 500 alternative bases shows Φ is
the only mathematical constant achieving complete reproduction at 10%
tolerance.

## Quick Start

```bash
# Clone and install dependencies
git clone https://github.com/MarcT123/Paper1-golden-ratio-fermion-masses.git
cd Paper1-golden-ratio-fermion-masses
pip install -r requirements.txt

# Run all statistical tests (produces terminal output)
python scripts/master_verification_v9.py

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
├── LICENSE
├── paper/
│   └── paper1.tex              # Full LaTeX source
├── scripts/
│   ├── master_verification_v9.py   # All Level A, B, D tests (742 lines)
│   └── regime_hypothesis_v9.py     # Regime hypothesis, ROC, chiral (282 lines)
└── figures/
    └── make_figures.py             # Generates Fig. 1–3 as PDF
```

## What the Scripts Reproduce

### master_verification_v9.py

| Test | Description | Expected p-value |
|------|-------------|-----------------|
| A1 | Genesis chain (MC, 3 leptons) | p = 0.003 |
| A2 | Phase uniformity (Greenwood, 6 ladders) | p = 0.007 |
| A3 | Φ-universality (200k random K) | p = 0.001 |
| A4 | Regime hypothesis (Fisher exact, 9 hadrons) | p = 0.008 |
| B1 | Ω⁻ = s × Φ⁶ (MC, Bonferroni) | p = 0.028 |
| B2 | Decuplet amplification ≈ Φ (Bonferroni) | p = 0.071 |
| B3 | ΔM(u→b) ≈ τ × Φ² (4 baryon pairs) | p = 0.001 |
| Fisher | Combined Level B (χ², 6 dof) | p = 0.0002 |
| D1 | Holdout: 167 PDG hadrons | p = 0.83 (falsified) |
| D2 | Gauge bosons (W, Z, H) | all p > 0.5 |

### regime_hypothesis_v9.py

- Fisher exact test on 9-hadron sample with QCD dominance parameter R
- ROC analysis identifying optimal threshold R_crit = 0.5
- Chiral derivation: Σ_s − Σ_ud = Δm(MS) × (Φ−1) ≈ 56 MeV

### make_figures.py

- **Figure 1:** Predicted vs. observed masses (log-log)
- **Figure 2:** Regime hypothesis — perfect separation at R = 0.5
- **Figure 3:** Φ-universality scan over K ∈ [1.2, 2.5]

## Verification Checkpoints

After running `master_verification_v9.py`, verify these values:

| Checkpoint | Expected |
|-----------|----------|
| Genesis R2 mean error | 1.81% |
| Φ-universality p | 0.001 |
| Phase uniformity (Greenwood) p | 0.007 |
| Regime hypothesis (Fisher) p | 0.008 |
| Fisher Combined (Level B) p | 0.0002 |
| Holdout 167 hadrons p | 0.83 |
| Bosons (W, Z, H) p | all > 0.5 |

## Reproducibility

All random number generators are seeded (`seed=42`) for exact
bitwise reproduction. The code requires only Python 3.8+ with NumPy
and SciPy. No custom packages or GPU required. Total runtime < 2 minutes
on a standard laptop.

## PDG Data

All particle masses are from PDG 2024 (Review of Particle Physics).
Mass conventions follow Table 1 of the paper.

## License

MIT License. See LICENSE file.

## Contact

Marc Tromel — marc.troe@gmail.com
Independent Researcher, Nürtingen, Germany
