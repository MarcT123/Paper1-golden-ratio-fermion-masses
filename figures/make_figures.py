#!/usr/bin/env python3
"""
Generate publication-quality figures for Paper 1.

Companion code for:
M. Trömel, "Statistical evidence for golden ratio structure in
elementary fermion mass ratios with QCD regime dependence"

Output (figure numbers match paper):
    fig1_pred_vs_obs.pdf   — Paper Figure 1
    fig2_universality.pdf  — Paper Figure 2
    fig3_regime.pdf        — Paper Figure 3

Usage:
    cd figures
    python make_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHI = (1 + np.sqrt(5)) / 2

# Output directory = directory containing this script
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# PDG 2024
PDG = {
    'e': 0.51099895, 'u': 2.16, 'd': 4.67, 's': 93.4,
    'mu': 105.6583755, 'pi': 139.57, 'p': 938.272,
    'c': 1270.0, 'tau': 1776.86, 'b': 4180.0, 't': 172760.0
}

PRED = {
    'u': PDG['e']*PHI**3, 'd': PDG['u']*PHI**PHI,
    's': np.sqrt(PDG['u']*PDG['b']), 'mu': PDG['e']*PHI**11,
    'pi': np.sqrt(PDG['d']*PDG['b']), 'p': PDG['pi']*PHI**4,
    'c': PDG['p']*PHI**(1/PHI), 'tau': PDG['c']*PHI**(1/PHI),
    'b': np.sqrt(PDG['s']*PDG['t']), 't': PDG['p']*PHI**11
}

# ========================================
# FIGURE 1: Predicted vs Observed masses
# ========================================
fig1, ax1 = plt.subplots(1, 1, figsize=(5.5, 5.5))

names = ['u','d','s','mu','pi','p','c','tau','b','t']
labels = ['$u$','$d$','$s$','$\\mu$','$\\pi$','$p$','$c$','$\\tau$','$b$','$t$']
pdg_vals = [PDG[n] for n in names]
pred_vals = [PRED[n] for n in names]

ax1.plot([0.3, 3e5], [0.3, 3e5], 'k--', lw=0.8, alpha=0.5,
         label='$M_\\mathrm{pred} = M_\\mathrm{PDG}$')
xx = np.logspace(-0.5, 5.5, 100)
ax1.fill_between(xx, xx*0.95, xx*1.05, alpha=0.08, color='blue',
                 label='$\\pm 5\\%$ band')
ax1.scatter(pdg_vals, pred_vals, s=60, c='#1F3864', zorder=5,
            edgecolors='white', linewidths=0.5)

for i, (name, x, y) in enumerate(zip(labels, pdg_vals, pred_vals)):
    offset = (8, 8)
    if name == '$\\mu$': offset = (-15, 8)
    if name == '$\\pi$': offset = (8, -12)
    if name == '$s$': offset = (-12, 8)
    ax1.annotate(name, (x, y), textcoords="offset points",
                 xytext=offset, fontsize=10)

ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('$M_\\mathrm{PDG}$ (MeV)', fontsize=12)
ax1.set_ylabel('$M_\\mathrm{pred}$ (MeV)', fontsize=12)
ax1.set_title('Genesis predictions vs. PDG 2024', fontsize=13, pad=10)
ax1.legend(fontsize=9, loc='upper left')
ax1.set_xlim(0.3, 3e5); ax1.set_ylim(0.3, 3e5)
ax1.set_aspect('equal')
ax1.grid(True, alpha=0.2, which='both')
fig1.tight_layout()

path1 = os.path.join(OUT_DIR, 'fig1_pred_vs_obs.pdf')
fig1.savefig(path1, bbox_inches='tight', dpi=300)
print(f"Figure 1 saved: {path1}")

# ========================================
# FIGURE 2: Phi-universality (K-scan)
# ========================================
fig2, ax2 = plt.subplots(1, 1, figsize=(6, 3.5))

np.random.seed(42)
relations = [
    ('$m_u/m_e$', PDG['u']/PDG['e']),
    ('$m_\\mu/m_e$', PDG['mu']/PDG['e']),
    ('$m_\\tau/m_e$', PDG['tau']/PDG['e']),
    ('$M_{\\Omega}/m_s$', 1672.45/PDG['s']),
    ('$\\Delta M/\\Delta m$', 146.8 / (PDG['s'] - (PDG['u']+PDG['d'])/2)),
]

K_vals = np.linspace(1.2, 2.5, 2000)
max_deltas = []
for K in K_vals:
    logK = np.log(K)
    deltas = []
    for _, ratio in relations:
        n = np.log(ratio) / logK
        deltas.append(abs(n - round(n)))
    max_deltas.append(max(deltas))

ax2.plot(K_vals, max_deltas, 'k-', lw=0.6, alpha=0.6)
ax2.axhline(y=0.080, color='#C00000', ls='--', lw=1, alpha=0.7,
            label='$\\Phi$ threshold (0.080)')
ax2.axvline(x=PHI, color='#1F3864', ls='-', lw=2, alpha=0.8,
            label=f'$\\Phi = {PHI:.3f}$')
ax2.plot(PHI, 0.080, 'o', color='#1F3864', markersize=8, zorder=5)

ax2.set_xlabel('Constant $K$', fontsize=12)
ax2.set_ylabel('$\\max|\\Delta n|(K)$', fontsize=12)
ax2.set_title('$\\Phi$-universality: no alternative $K$ matches all 5 relations',
              fontsize=11, pad=8)
ax2.legend(fontsize=9, loc='upper right')
ax2.set_ylim(0, 0.5); ax2.set_xlim(1.2, 2.5)
ax2.grid(True, alpha=0.2)
fig2.tight_layout()

path2 = os.path.join(OUT_DIR, 'fig2_universality.pdf')
fig2.savefig(path2, bbox_inches='tight', dpi=300)
print(f"Figure 2 saved: {path2}")

# ========================================
# FIGURE 3: Regime hypothesis
# ========================================
fig3, ax3 = plt.subplots(1, 1, figsize=(6, 4))

hadrons = [
    ("$\\Delta(1232)$", 1232, 8.99, True),
    ("$\\Sigma^*(1385)$", 1383.7, 97.73, True),
    ("$\\Xi^*(1532)$", 1531.8, 188.96, True),
    ("$\\Omega^-$", 1672.45, 280.2, True),
    ("$\\Lambda_c$", 2286.46, 1276.83, False),
    ("$\\Lambda_b$", 5619.60, 4186.83, False),
    ("$\\Omega_b^-$", 6046.0, 4366.8, False),
    ("$\\Omega_{ccc}$", 4761.0, 3810.0, False),
    ("$\\Omega_{bbb}$", 14371.0, 12540.0, False),
]

for name, M, sumq, phi in hadrons:
    R = (M - sumq) / M
    color = '#1F3864' if phi else '#C00000'
    marker = 'o' if phi else 'x'
    ms = 10 if phi else 9
    lw = 1 if phi else 2
    ax3.plot(R, 0.5, marker, color=color, markersize=ms, markeredgewidth=lw)
    yoff = 0.08 if phi else -0.08
    ax3.annotate(name, (R, 0.5+yoff), fontsize=8.5, ha='center', va='center')

ax3.axvline(x=0.5, color='gray', ls='--', lw=1.2, alpha=0.7)
ax3.annotate('$R_{\\mathrm{crit}} = 0.5$', (0.5, 0.85), fontsize=10,
             ha='center', color='gray')
ax3.text(0.25, 0.25, 'Quark-dominated\n(no $\\Phi$)', fontsize=10,
         ha='center', color='#C00000', style='italic')
ax3.text(0.80, 0.25, 'QCD-dominated\n($\\Phi$ confirmed)', fontsize=10,
         ha='center', color='#1F3864', style='italic')
ax3.annotate('', xy=(0.83, 0.65), xytext=(0.44, 0.65),
            arrowprops=dict(arrowstyle='<->', color='#666666', lw=1.5))
ax3.text(0.635, 0.69, '$\\Delta R = 0.39$', fontsize=9,
         ha='center', color='#666666')

ax3.set_xlabel('QCD dominance parameter $R = E_\\mathrm{bind}/M$', fontsize=12)
ax3.set_xlim(0.05, 1.05); ax3.set_ylim(0.1, 0.95)
ax3.set_yticks([])
ax3.set_title('Regime hypothesis: perfect separation ($p = 0.008$)',
              fontsize=12, pad=10)
ax3.grid(True, axis='x', alpha=0.2)
fig3.tight_layout()

path3 = os.path.join(OUT_DIR, 'fig3_regime.pdf')
fig3.savefig(path3, bbox_inches='tight', dpi=300)
print(f"Figure 3 saved: {path3}")

print("\nAll figures generated successfully.")
