# Statistical evidence for golden ratio structure in elementary fermion mass ratios

This repository contains the paper, verification code, and figures for:

> M. Trömel, *Statistical evidence for golden ratio structure in elementary fermion mass ratios with QCD regime dependence* (2026).

## Repository structure

```
Paper1-golden-ratio-fermion-masses/
├── LICENSE
├── README.md
├── requirements.txt
├── requirements-pinned.txt
├── releases/
│   └── v1.0.0-original.md
├── paper/
│   ├── paper1.tex          # LaTeX source
│   └── paper1.pdf          # Compiled paper (16 pages)
├── scripts/
│   ├── master_verification_v10.py   # All statistical tests (1171 lines)
│   └── regime_hypothesis_v9.py      # Regime hypothesis analysis (288 lines)
└── figures/
    ├── fig1_pred_vs_obs.pdf         # Genesis predictions vs. PDG masses
    ├── fig2_universality.pdf        # Φ-universality scan
    └── fig3_regime.pdf              # QCD regime separation
```

## Quick start

```bash
pip install -r requirements.txt
cd scripts
python master_verification_v10.py
```

All random number generators are seeded (`seed=42`) for reproducibility.

## Key results

| Test | Result | p-value |
|------|--------|---------|
| A1: Genesis chain | 11 masses from electron, mean error 2.5% | — |
| A2: Phase uniformity | 6 Φ-ladders, Greenwood test | 0.007 |
| A3: Φ-universality | No alternative K in [1.2, 2.5] | 0.001 |
| A4: Regime hypothesis | Perfect 4/4 vs 0/5 separation at R = 0.5 | 0.008 |
| A5: Chain connectivity | 10/11 at 5%, strengthens at tighter tolerance | 0.006 |
| B1–B3: Hadron traces | Fisher combined | 2.2 × 10⁻⁵ |
| D1: 167-hadron holdout | No general Φ-grid proximity | 0.055 (n.s.) |
| D2: Gauge bosons | No Φ-signal | all p > 0.5 |

## Compiling the paper

```bash
cd paper
pdflatex paper1.tex
pdflatex paper1.tex   # second pass for references
```

Figures are loaded from `../figures/` via `\graphicspath`.

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

M. Trömel — marc.troe@gmail.com
