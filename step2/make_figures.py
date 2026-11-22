"""Publication figures for the ShapeShift write-up (reads step2/results/all_results.csv; fig 3 also needs the
Tahoe signature tables produced by step2/tahoe_signatures.py).

    python step2/make_figures.py            -> figures/fig1..fig4 *.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "step2" / "results" / "all_results.csv"
FIG = ROOT / "figures"
TIERS = [("T1_LINCS_A549", "T1  LINCS\nA549 (= Tahoe line)"), ("T2_EUOS_HepG2", "T2  EU-OPENSCREEN\nHepG2 (≈ C3A)"),
         ("T3_JUMP_U2OS", "T3  JUMP\nU2OS (not in Tahoe)")]
C = {"primary": "#5f6368", "centred": "#1a73e8", "twoplate": "#e8710a", "ceiling": "#c5221f", "one": "#9aa0a6"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def row(r, tier, label, analysis, **kw):
    s = r[(r.tier == tier) & (r.label == label) & (r.analysis == analysis)]
    for k, v in kw.items():
        s = s[s[k] == v]
    if s.empty:
        raise RuntimeError(f"missing result row {tier}/{label}/{analysis}/{kw}")
    return s.iloc[0]


def fig1(r):
    """Agreement per tier: pre-registered primary, plate-centred (all units), two-plate units; 95% bootstrap CIs."""
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    series = [("Pre-registered primary\n(annotated, within-MoA null)", "PRIMARY", "annotated_within_moa", C["primary"]),
              ("Plate-centred, all units\n(secondary)", "PRIMARY_CENTERED", "all_permutation", C["centred"]),
              ("Plate-centred, units on two\nTahoe plates (secondary)", "SENS_multiplate_units", "all_permutation", C["twoplate"])]
    w = 0.25
    for j, (lab, label, an, col) in enumerate(series):
        for i, (t, _) in enumerate(TIERS):
            q = row(r, t, label, an)
            x = i + (j - 1) * w
            ax.errorbar(x, q.r, yerr=[[q.r - q.ci_low], [q.ci_high - q.r]], fmt="o", c=col, capsize=3, ms=5,
                        label=lab if i == 0 else None)
            ax.text(x, q.ci_high + 0.008, f"p={q.p:.2g}", ha="center", fontsize=6.5, color=col)
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xticks(range(3), [b for _, b in TIERS])
    ax.set_ylabel("Mantel r (morphology vs transcriptome\ncompound-similarity matrices)")
    ax.set_title("Agreement is near zero after correction; it rises where Tahoe is measured twice", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_agreement_per_tier.png", dpi=180)


def fig2(r):
    """Observed agreement vs noise ceiling on the same (two-plate) units; disattenuated ratio annotated."""
    fig, ax = plt.subplots(figsize=(6.4, 3.7))
    for i, (t, name) in enumerate(TIERS):
        ceil = row(r, t, "CEILING", "max_attainable_mantel_r_multiplate_units")
        obs2 = row(r, t, "SENS_multiplate_units", "all_permutation")
        obs1 = row(r, t, "SENS_singleplate_units", "all_permutation")
        ax.bar(i, ceil.r, width=0.6, color=C["ceiling"], alpha=0.18, edgecolor=C["ceiling"],
               label="noise ceiling (same two-plate units)" if i == 0 else None)
        ax.errorbar(i - 0.12, obs2.r, yerr=[[obs2.r - obs2.ci_low], [obs2.ci_high - obs2.r]], fmt="o", c=C["twoplate"],
                    capsize=3, label="observed, two-plate units" if i == 0 else None)
        ax.errorbar(i + 0.12, obs1.r, yerr=[[obs1.r - obs1.ci_low], [obs1.ci_high - obs1.r]], fmt="s", c=C["one"],
                    capsize=3, label="observed, one-plate units" if i == 0 else None)
        ax.text(i, ceil.r + 0.015, f"disattenuated\n≈{obs2.r / ceil.r:.2f}", ha="center", fontsize=7)
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xticks(range(3), [b for _, b in TIERS])
    ax.set_ylabel("Mantel r")
    ax.set_ylim(-0.1, 0.62)
    ax.set_title("Observed agreement is 15–35% of what the data could support", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_noise_ceiling.png", dpi=180)


def fig3():
    """Tahoe plate / shared-DMSO effect: cosine of consensus signatures for same-plate vs different-plate pairs,
    before and after plate centring (5 uM, JUMP working-set units)."""
    sys.path.insert(0, str(ROOT / "step2"))
    import compare_tiers as ct
    tah = ct.Tahoe()
    units = pd.read_csv(ROOT / "step2" / "units.csv")
    us = sorted(set(units.loc[units.tier == "T3_JUMP_U2OS", "unit"]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, centred, title in zip(axes, (False, True), ("raw signatures", "plate-centred signatures")):
        s = tah.sig("consensus", 5.0, centered=centred)
        s = s[s.unit.isin(us)].reset_index(drop=True)
        X = s[tah.genes].to_numpy(np.float64)
        S = ct.cos(X)
        P = [set(str(p).split(";")) for p in s.plates]
        i, j = np.triu_indices(len(s), 1)
        same = np.array([bool(P[a] & P[b]) for a, b in zip(i, j)])
        v = S[i, j]
        bins = np.linspace(-0.6, 1.0, 60)
        ax.hist(v[~same], bins=bins, density=True, alpha=0.55, color=C["one"], label=f"different plate (mean {v[~same].mean():.2f})")
        ax.hist(v[same], bins=bins, density=True, alpha=0.55, color=C["ceiling"], label=f"same plate (mean {v[same].mean():.2f})")
        ax.set_title(f"Tahoe 5 µM consensus, {title}", fontsize=9)
        ax.set_xlabel("cosine similarity between compounds")
        ax.legend(fontsize=7, frameon=False)
    axes[0].set_ylabel("density")
    fig.suptitle("Every drug on a Tahoe plate shares the same DMSO control: plate centring removes the offset", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_plate_effect.png", dpi=180)


def fig4(r):
    """T1 dose curve: agreement vs LINCS concentration (nearest Tahoe dose), A549 and consensus signatures."""
    d = r[(r.tier == "T1_LINCS_A549") & (r.label == "DOSE_matched") & (r.analysis == "all_permutation")]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for sig, col in (("A549", C["centred"]), ("consensus", C["twoplate"])):
        s = d[d.tahoe_sig == sig].sort_values("morph_dose")
        ax.errorbar(s.morph_dose, s.r, yerr=[s.r - s.ci_low, s.ci_high - s.r], fmt="o-", c=col, capsize=3,
                    label=f"Tahoe {sig} signature (nearest Tahoe dose)")
    for x in (0.0412, 0.37, 3.33):
        ax.axvspan(x / 1.06, x * 1.06, color="#34a853", alpha=0.15, lw=0)
    ax.set_xscale("log")
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xlabel("LINCS concentration (µM); green = within 1.5× of a Tahoe dose")
    ax.set_ylabel("Mantel r, all units")
    ax.set_title("T1 (A549): agreement is flat across concentration", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig4_dose_curve.png", dpi=180)


def main():
    FIG.mkdir(exist_ok=True)
    r = pd.read_csv(RES)
    fig1(r), fig2(r), fig4(r)
    if "--no-tahoe" not in sys.argv:
        fig3()
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
