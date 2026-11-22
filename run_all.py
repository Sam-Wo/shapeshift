"""One-command rerun of the full ShapeShift pipeline (Steps 1-2).

    python run_all.py                 # everything (streams ~22.5 GB, peak disk <= ~7.5 GB, several hours)
    python run_all.py --from compare  # only re-run the statistics + figures from the cached tables
    python run_all.py --list          # show the steps

Each step downloads exactly what it needs from the public sources (Hugging Face, Cell Painting Gallery S3,
Zenodo, GitHub, PubChem) and raises if a source is unreachable - no tier is skipped silently.
"""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    ("overlap", "overlap.py", "Step 1: match Tahoe-100M drugs to JUMP-CP compounds (InChIKey / tautomer tiers)"),
    ("celllines", "cellcheck/obs_cell_lines.py", "which 50 Tahoe cell lines were screened (range reads)"),
    ("cellmatch", "cellcheck/cellline_match.py", "match Tahoe drugs to LINCS / EU-OPENSCREEN / CPJUMP1"),
    ("units", "step2/build_units.py", "compound units per tier + T3 working set"),
    ("tahoe", "step2/tahoe_signatures.py", "Tahoe pseudobulk signatures (streams ~15 GB of row groups)"),
    ("perplate", "step2/tahoe_signatures.py --stage perplate", "Tahoe per-plate signatures (from cache; noise ceiling)"),
    ("gate", "step2/gate_t2.py", "T2 HepG2/C3A cell-count gate"),
    ("jump", "step2/jump_profiles.py", "T3 JUMP U2OS morphology (416 plates)"),
    ("lincs", "step2/morph_lincs.py", "T1 LINCS A549 morphology (136 plates, 6 doses)"),
    ("euos", "step2/morph_euos.py", "T2 EU-OPENSCREEN HepG2 morphology (Zenodo range reads)"),
    ("cpjump1", "step2/morph_cpjump1.py", "CPJUMP1 A549-vs-U2OS calibration"),
    ("compare", "step2/compare_tiers.py", "all statistics (~25 min)"),
    ("figures", "step2/make_figures.py", "publication figures"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", default=STEPS[0][0], choices=[s[0] for s in STEPS])
    ap.add_argument("--only", choices=[s[0] for s in STEPS])
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for k, f, d in STEPS:
            print(f"{k:10s} {f:32s} {d}")
        return
    names = [s[0] for s in STEPS]
    todo = [STEPS[names.index(a.only)]] if a.only else STEPS[names.index(a.start):]
    for k, f, d in todo:
        print(f"\n=== [{k}] {d}\n", flush=True)
        script, *extra = f.split()
        subprocess.run([sys.executable, "-X", "utf8", str(ROOT / script), *extra], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
