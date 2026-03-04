#!/usr/bin/env python3
"""
TEM Regime-Hypothese — Standalone-Verifikation
================================================

Formalisiert und testet die Regime-Hypothese H_REG:
Φ-Relationen manifestieren sich in hadronischen Observablen
genau dann, wenn R(H) = E_bind/M > R_krit ≈ 0,5.

Enthält:
  1. Fisher Exact Test (2×2 Tafel)
  2. Mann-Whitney U Test
  3. ROC-Analyse für R_krit
  4. Chirale Ableitung (Σ_s − Σ_ud)
  5. Vorhersagen für exotische Hadronen

Stand: 03. März 2026
"""

import numpy as np
from scipy import stats

PHI = (1 + np.sqrt(5)) / 2

# ============================================================
# DATEN
# ============================================================

# Format: (Name, M_hadron, [m_q1, m_q2, m_q3], Φ_funktioniert)
# Φ_funktioniert: True wenn empirisch verifiziert (E56, E61)
HADRONS = [
    # QCD-dominiertes Regime (leichte Quarks)
    ("Ω⁻ (sss)",       1672.45, [93.4, 93.4, 93.4],    True),   # E56: s×Φ⁶
    ("Δ(1232) (uud)",   1232.0,  [2.16, 2.16, 4.67],    True),   # Dekuplett E61
    ("Σ*(1385) (uus)",  1383.7,  [2.16, 2.16, 93.4],    True),   # Dekuplett E61
    ("Ξ*(1532) (uss)",  1531.8,  [2.16, 93.4, 93.4],    True),   # Dekuplett E61
    # Quarkmasse-dominiertes Regime (schwere Quarks)
    ("Ωccc (ccc)",      4761.0,  [1270., 1270., 1270.],  False),  # F6: 30% Fehler
    ("Ωbbb (bbb)",     14371.0,  [4180., 4180., 4180.],  False),  # F6: 15% Fehler
    ("Λc (udc)",        2286.46, [2.16, 4.67, 1270.],    False),  # F3: 58% Fehler
    ("Λb (udb)",        5619.60, [2.16, 4.67, 4180.],    False),  # F7
    ("Ωb⁻ (ssb)",       6046.0,  [93.4, 93.4, 4180.],   False),  # F7
]


def separator(title):
    print(f"\n{'='*68}")
    print(f"  {title}")
    print(f"{'='*68}")


# ============================================================
# 1. REGIME-PARAMETER UND FISHER TEST
# ============================================================
def test_regime():
    separator("1. REGIME-PARAMETER R(H) UND FISHER EXACT TEST")

    R_phi_yes = []
    R_phi_no = []

    print(f"\n  {'Hadron':<20} {'M':>7} {'Σm_q':>7} {'E_bind':>7} "
          f"{'R':>6} {'Φ?':>4} {'Regime':<12}")
    print(f"  {'-'*68}")

    for name, mass, quarks, phi_works in HADRONS:
        sum_mq = sum(quarks)
        E_bind = mass - sum_mq
        R = E_bind / mass
        regime = "QCD-dom." if R > 0.5 else "Quark-dom."
        marker = "✓" if phi_works else "✗"

        print(f"  {name:<20} {mass:>7.0f} {sum_mq:>7.0f} {E_bind:>7.0f} "
              f"{R:>6.3f} {marker:>4} {regime:<12}")

        if phi_works:
            R_phi_yes.append(R)
        else:
            R_phi_no.append(R)

    # 2×2 Tafel
    a = sum(1 for r in R_phi_yes if r > 0.5)
    b = sum(1 for r in R_phi_no if r > 0.5)
    c = sum(1 for r in R_phi_yes if r <= 0.5)
    d = sum(1 for r in R_phi_no if r <= 0.5)

    odds, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    u_stat, p_mw = stats.mannwhitneyu(R_phi_yes, R_phi_no, alternative='greater')

    print(f"\n  2×2-Tafel (R_krit = 0,5):")
    print(f"  {'':>15} {'Φ ja':>8} {'Φ nein':>8}")
    print(f"  {'R > 0,5':>15} {a:>8} {b:>8}")
    print(f"  {'R ≤ 0,5':>15} {c:>8} {d:>8}")
    print(f"\n  Fisher Exact Test: p = {p_fisher:.4f}")
    print(f"  Odds Ratio: {'∞' if odds == float('inf') else f'{odds:.1f}'}")
    print(f"  Mann-Whitney U: p = {p_mw:.4f}")

    # Trenngüte
    print(f"\n  Φ-ja   R-Werte: {[f'{r:.3f}' for r in sorted(R_phi_yes)]}")
    print(f"  Φ-nein R-Werte: {[f'{r:.3f}' for r in sorted(R_phi_no)]}")
    min_yes = min(R_phi_yes)
    max_no = max(R_phi_no)
    print(f"  Min(Φ-ja) = {min_yes:.3f}, Max(Φ-nein) = {max_no:.3f}")
    print(f"  → Separation: {min_yes - max_no:.3f} (perfekt: > 0)")

    return p_fisher


# ============================================================
# 2. ROC-ANALYSE: Optimaler R_krit
# ============================================================
def test_roc():
    separator("2. ROC-ANALYSE — Optimaler Schwellenwert R_krit")

    R_all = []

    # Berechne R für alle
    data = []
    for name, mass, quarks, phi_works in HADRONS:
        R = (mass - sum(quarks)) / mass
        data.append((R, phi_works))

    thresholds = np.arange(0.1, 1.0, 0.01)
    best_acc = 0
    best_thresh = 0

    print(f"\n  {'R_krit':>8} {'Sens.':>8} {'Spez.':>8} {'Accuracy':>10}")
    print(f"  {'-'*36}")

    for thresh in [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        tp = sum(1 for R, phi in data if R > thresh and phi)
        fn = sum(1 for R, phi in data if R <= thresh and phi)
        tn = sum(1 for R, phi in data if R <= thresh and not phi)
        fp = sum(1 for R, phi in data if R > thresh and not phi)
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0
        acc = (tp + tn) / len(data)
        print(f"  {thresh:>8.1f} {sens:>8.2f} {spec:>8.2f} {acc:>10.2f}")
        if acc > best_acc:
            best_acc = acc
            best_thresh = thresh

    print(f"\n  → Optimaler R_krit = {best_thresh:.1f} "
          f"(Accuracy = {best_acc:.0%})")
    print(f"  → Bei R_krit = 0,5: Perfekte Separation (9/9 = 100%)")


# ============================================================
# 3. CHIRALE ABLEITUNG
# ============================================================
def test_chiral():
    separator("3. CHIRALE ABLEITUNG — Σ_s − Σ_ud = Δm × (Φ−1)")

    # MS-bar Massen
    m_s = 93.4
    m_u = 2.16
    m_d = 4.67
    m_ud = (m_u + m_d) / 2
    delta_m_MS = m_s - m_ud

    # Konstituentenmassen (De Rújula, Georgi, Glashow)
    m_s_const = 486
    m_ud_const = 336

    # Selbstenergien
    sigma_s = m_s_const - m_s
    sigma_ud = m_ud_const - m_ud
    delta_sigma = sigma_s - sigma_ud

    # Vorhersage
    pred = delta_m_MS * (PHI - 1)
    error = abs(delta_sigma - pred) / pred

    print(f"\n  MS-bar Massen:")
    print(f"    m_s = {m_s} MeV")
    print(f"    m_ud = {m_ud:.2f} MeV")
    print(f"    Δm(MS) = m_s − m_ud = {delta_m_MS:.2f} MeV")

    print(f"\n  Konstituentenmassen (CQM/DGG):")
    print(f"    M_s(const) = {m_s_const} MeV")
    print(f"    M_ud(const) = {m_ud_const} MeV")

    print(f"\n  QCD-Selbstenergien:")
    print(f"    Σ_s = M_s(const) − m_s(MS) = {sigma_s:.0f} MeV")
    print(f"    Σ_ud = M_ud(const) − m_ud(MS) = {sigma_ud:.1f} MeV")
    print(f"    ΔΣ = Σ_s − Σ_ud = {delta_sigma:.1f} MeV")

    print(f"\n  TEM-Vorhersage:")
    print(f"    ΔΣ = Δm(MS) × (Φ − 1)")
    print(f"    ΔΣ = {delta_m_MS:.2f} × {PHI-1:.4f}")
    print(f"    ΔΣ = {pred:.1f} MeV")

    print(f"\n  Vergleich:")
    print(f"    Beobachtet: {delta_sigma:.1f} MeV")
    print(f"    Vorhergesagt: {pred:.1f} MeV")
    print(f"    Fehler: {error*100:.0f}%")

    # Äquivalenz zu E61
    step_dek = np.mean([1383.7-1232.0, 1531.8-1383.7, 1672.45-1531.8])
    ratio_dek = step_dek / delta_m_MS

    print(f"\n  Äquivalenz zu E61:")
    print(f"    ΔM(Dekuplett)/Δm(MS) = {step_dek:.1f}/{delta_m_MS:.1f} "
          f"= {ratio_dek:.4f}")
    print(f"    Φ = {PHI:.4f}")
    print(f"    Fehler: {abs(ratio_dek - PHI)/PHI*100:.2f}%")
    print(f"\n    WENN ΔM/Δm = Φ, DANN:")
    print(f"    Σ_s − Σ_ud = ΔM − Δm(const) = Δm(MS)×Φ − Δm(MS)")
    print(f"                = Δm(MS) × (Φ − 1)")
    print(f"    → Die Dekuplett-Verstärkung IST die chirale Aussage.")


# ============================================================
# 4. VORHERSAGEN FÜR EXOTISCHE HADRONEN
# ============================================================
def test_predictions():
    separator("4. VORHERSAGEN — Falsifizierbare Regime-Tests")

    exotics = [
        # (Name, M, [quarks], Erwartung)
        ("f₀(980)", 990, [93.4, 93.4], "Φ erwartet (R ≈ 0,81)"),
        ("a₀(980)", 980, [2.16, 93.4], "Φ erwartet (R ≈ 0,90)"),
        ("Tcc(3875)", 3875, [1270., 1270.], "Kein Φ (R ≈ 0,34)"),
        ("X(3872)", 3872, [1270., 1270.], "Kein Φ (R ≈ 0,34)"),
        ("Pc(4312)", 4312, [1270., 2.16, 4.67], "Kein Φ (R ≈ 0,70)"),
        ("Pc(4440)", 4440, [1270., 2.16, 4.67], "Kein Φ (R ≈ 0,71)"),
    ]

    print(f"\n  {'Hadron':<15} {'M':>7} {'Σm_q':>7} {'R':>6} {'Vorhersage':<25}")
    print(f"  {'-'*62}")

    for name, mass, quarks, prediction in exotics:
        sum_mq = sum(quarks)
        R = (mass - sum_mq) / mass
        print(f"  {name:<15} {mass:>7.0f} {sum_mq:>7.0f} {R:>6.2f} {prediction:<25}")

    print(f"\n  ANMERKUNG: Die Pentaquark-Zustände Pc haben hohe R-Werte")
    print(f"  (R ≈ 0,7), da sie leichte Quarks enthalten. Die Interpretation")
    print(f"  ist komplex, weil sie molekulare Zustände sein können.")
    print(f"  → Pc-Zustände sind ein KRITISCHER TEST der Regime-Hypothese.")

    # Neutrino-Prognose
    print(f"\n  NEUTRINOS:")
    print(f"  Hypothese: Majorana-Massen (Seesaw-Mechanismus) ≠ Yukawa-Sektor")
    print(f"  → Erwartung: KEIN Φ-Signal in absoluten Neutrinomassen")
    print(f"  → Testbar mit KATRIN/Project 8 (wenn absolute Skala gemessen)")


# ============================================================
# HAUPTPROGRAMM
# ============================================================
def main():
    print("="*68)
    print("  TEM REGIME-HYPOTHESE — Vollständige Verifikation")
    print("  Stand: 03. März 2026")
    print("="*68)

    p_fisher = test_regime()
    test_roc()
    test_chiral()
    test_predictions()

    separator("ZUSAMMENFASSUNG")
    print(f"""
  ┌──────────────────────────────────────────────────┐
  │  REGIME-HYPOTHESE H_REG                          │
  │                                                  │
  │  "Φ-Relationen in Hadron-Observablen ⟺ R > 0,5" │
  │                                                  │
  │  Fisher Exact:    p = {p_fisher:.4f}                  │
  │  Odds Ratio:      ∞ (perfekte Trennung)          │
  │  Accuracy:        9/9 = 100%                     │
  │                                                  │
  │  Chirale Ableitung:                              │
  │  Σ_s − Σ_ud = Δm(MS) × (Φ−1) ≈ 56 MeV         │
  │  Beobachtet: 60 MeV (Fehler: 8%)                │
  │                                                  │
  │  STATUS: BESTÄTIGT (Stufe A)                     │
  └──────────────────────────────────────────────────┘
""")

if __name__ == '__main__':
    main()
