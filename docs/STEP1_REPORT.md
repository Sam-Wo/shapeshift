# ShapeShift – Step 1 report: Tahoe-100M × JUMP-CP compound overlap

*Generated 2026-11-22 by `overlap.py` (metadata only, ~17 MB in `data/`). Per-drug detail: `overlap.csv` (379 rows).*

## Verdict: **GO**

- **Tier 1 alone (InChIKey match): 315 / 379 Tahoe drugs, 292 unique JUMP compounds, 21 MoA classes with ≥3 compounds.** That is about 6× the go bar.
- **With tier 2b** (the same molecule drawn as a different tautomer in JUMP, flagged for review): **357 drugs, 331 unique compounds, 22 MoA classes.**

The verdict does not depend on tier 2b.

The real constraint is **dose and replicates**, not overlap. JUMP compound plates were run at one concentration per source:

- **10 µM** in every source except source_7: 189 matched compounds.
- **0.625 µM** in source_7: 315 matched compounds.

Replicate depth is low (median 5–7 wells per compound per dose). See the caveats below.

## Inputs

| Source | File | Notes |
|---|---|---|
| Tahoe-100M (HF `tahoebio/Tahoe-100M`) | `metadata/drug_metadata.parquet` (40 kB) | 379 drugs: name, SMILES, PubChem CID, targets, `moa-broad`, `moa-fine` |
| | `metadata/sample_metadata.parquet` (65 kB) | Doses. Every drug has 3 doses: **0.05, 0.5, 5 µM** |
| | `metadata/cell_line_metadata.parquet` (19 kB) | Cell-line check only. The 50 lines include no U2OS and no osteosarcoma |
| JUMP (`jump-cellpainting/datasets`) | `compound`, `well`, `plate`, `perturbation_control` | 115,796 compounds, 1.15 M wells |
| JUMP Repurposing Hub join | `repurposing_hub_annotations_processed.tsv` (cellpainting-gallery) | JCP2022 → name, MoA, target. Used for the name bridge, name checks and extra MoA |
| PubChem PUG REST | Properties for 376 of 377 CIDs; parent CIDs and name search for unmatched drugs | Cached in `data/pubchem_cache.json` |

`obs_metadata.parquet` (2.3 GB), expression data, profiles and images were **not** downloaded.

## Method

1. **Standardize** both sides with RDKit: cleanup, then the largest organic fragment (strips salts and solvates), then neutralize, then **remove stereo**, then compute the InChIKey.
   - Stereo is removed because **JUMP is stereo-flat**: 99.999% of its InChIKeys end in `-UHFFFAOYSA-`, and no JUMP SMILES contains `@`. Tahoe's SMILES are also flat (5/379 have `@`).
   - Each Tahoe drug is keyed from its own SMILES and from its PubChem CID structure. For all 377 drugs with both, the two keys agree on connectivity.
   - JUMP is matched on its own key and on the recomputed one.
2. **Tiers**, each assigned in order:
   - **Tier 1, full InChIKey on stereo-flat structures: 315 drugs.**
     - **1a, exact and achiral: 160.** PubChem's isomeric structure has no stereocenters, so the match is truly exact.
     - **1b, stereo unresolved: 155.** The drug is chiral, and JUMP's flat record cannot confirm which stereoisomer was plated.
   - **Tier 2a, InChIKey block 1: 0.** Once both sides are flat and neutralized, block 1 adds nothing beyond tier 1.
   - **Tier 2b, same formula + same heavy-atom graph: 42 drugs, flagged REVIEW.** This tier ignores bond orders and H placement (RDKit `ElementGraph` hash). It exists because **JUMP stores many well-known drugs as unusual tautomers or with shifted double bonds**, which changes even block 1:
     - sunitinib as a hydroxy-indole;
     - dexamethasone and the other corticosteroids as Δ5 dienes;
     - piroxicam and meloxicam as keto-imines;
     - rapamycin, everolimus, paclitaxel and docetaxel.

     JUMP's own control list confirms the approach: positive control `JCP2022_025848` is named "dexamethasone", but only tier 2b matches it. Requiring an identical formula keeps true analogs out; for example, belzutifan is correctly *not* matched to its 2,2-difluoro analog.
   - **Tier 3 (REVIEW).** Three checks run here:
     - a name bridge: the Tahoe name or PubChem title against JUMP-side names (Repurposing Hub and the control list);
     - a PubChem bridge: the parent CID for salts, plus name search;
     - a metal check: matches where standardization stripped a metal are moved to tier 3.

     The bridges found **no new drugs**. One match was demoted: **auranofin**, whose JUMP hit `JCP2022_082888` is the gold-free tetra-acetyl-thioglucose ligand.
3. **Well, plate and source counts** use only `COMPOUND` plates. Wells on TARGET1/2, ORF, CRISPR and POSCON plates are listed separately, so the positive controls don't inflate replicate counts.
4. **Concentrations** are nominal: JUMP publishes no per-well dose. The values come from the [JUMP Hub quirks page](https://github.com/broadinstitute/jump_hub/blob/main/explanations/quirks_details.md), checked against the page itself:
   - COMPOUND plates: 10 µM, except source_7 at 0.625 µM;
   - positive controls and TARGET1/2 plates: 5 µM.

## Summary counts

| | Drugs | Unique JUMP cpds |
|---|---|---|
| Tahoe drugs | 379 | – |
| Tier 1 (1a achiral exact 160 + 1b stereo unresolved 155) | **315** | **292** |
| Tier 1 + 2b | **357** | **331** |
| … of which approved (`human-approved = yes`) | 249 of 263 approved | |
| … at 10 µM (any source other than source_7)¹ | 200 | 189 |
| … at 0.625 µM (source_7)¹ | 341 | 315 |
| … at both doses¹ | 186 | 175 |
| Tier 3 (needs review) | 1 (auranofin) | – |
| Unmatched | 21 | – |

¹ Excludes the two JUMP positive controls, dexamethasone and quinidine. Their roughly 8,700 wells each are 5 µM control wells.

**Replicates on COMPOUND plates:**

- At 10 µM: median 7 wells (IQR 6–11).
- At 0.625 µM (source_7): median 5 wells (IQR 5–10).

**Source spread:** 155 drugs are in source_7 only. The rest are mostly in 5–7 sources. 23 drugs also appear on TARGET1/TARGET2 plates.

**Name check:** for 204 matches the Repurposing Hub name equals the Tahoe name, and for 115 the Hub has no annotation. That includes every tier-2b drug, because the Hub join was also done by InChIKey. All 38 "disagreements" are synonyms, codes or salt/stereo variants, for example capivasertib = AZD5363, belumosudil = KD025, sapanisertib = MLN0128.

## MoA coverage (Tahoe `moa-fine`, unique JUMP compounds)

I used `moa-fine`; `moa-broad` (inhibitor / agonist / unclear) is too coarse for this test.

| MoA class (Tahoe `moa-fine`) | Tahoe total | Tier 1+2b | Tier 1 only | at 10 µM¹ | at both doses¹ |
|---|---|---|---|---|---|
| DNA synthesis/repair inhibitor | 27 | 22 | 19 | 11 | 10 |
| Other TK inhibitor | 15 | 11 | 11 | 6 | 4 |
| Cyclooxygenase inhibitor | 13 | 11 | 9 | 6 | 6 |
| Other MAPK inhibitor | 11 | 11 | 11 | 11 | 11 |
| EGFR/ERBB inhibitor | 12 | 9 | 8 | 4 | 4 |
| Multi-TK inhibitor | 10 | 9 | 8 | 4 | 4 |
| Adrenoceptor agonist | 8 | 8 | 8 | 5 | 5 |
| MTOR inhibitor | 7 | 7 | 4 | 6 | 6 |
| PI3K/AKT inhibitor | 6 | 6 | 6 | 4 | 4 |
| Protein synthesis inhibitor | 6 | 6 | 6 | 5 | 5 |
| Microtubule inhibitor | 7 | 5 | 3 | 1 | 1 |
| JAK/STAT inhibitor | 6 | 5 | 5 | 4 | 4 |
| CDK inhibitor | 5 | 5 | 4 | 0 | 0 |
| Glucocorticoid receptor agonist | 5 | 5 | 1 | 0 | 0 |
| MEK inhibitor | 5 | 4 | 4 | 3 | 3 |
| Androgen receptor antagonist | 4 | 4 | 3 | 2 | 2 |
| DNA methyltransferase inhibitor | 4 | 4 | 4 | 1 | 1 |
| HDAC inhibitor | 4 | 4 | 4 | 3 | 3 |
| Retinoic receptor agonist | 4 | 4 | 4 | 1 | 0 |
| GSK3 inhibitor | 4 | 3 | 3 | 3 | 3 |
| RAF inhibitor | 3 | 3 | 3 | 3 | 2 |
| Sonic inhibitor | 3 | 3 | 3 | 3 | 3 |
| Glucose transporter inhibitor | 4 | 2 | 2 | 0 | 0 |
| Proteasome inhibitor | 3 | 2 | 2 | 2 | 1 |
| RAS inhibitor | 4 | 1 | 1 | 0 | 0 |

**Classes with ≥3 unique compounds:**

- 22 for tier 1+2b, 21 for tier 1 only;
- 16 at 10 µM, 21 at 0.625 µM;
- 15 among compounds measured at both doses.

**Note:** Tahoe labels 199/379 drugs `moa-fine = unclear`, including 188 of the matched drugs. For 126 of those, `overlap.csv` has a Repurposing Hub MoA (`jump_hub_moa`) that could be used to enlarge the classes later.

## Notable unmatched Tahoe drugs (21)

- **Approved oncology agents missing from JUMP:**
  - bortezomib (proteasome; ixazomib *is* matched);
  - oxaliplatin (a platinum complex, unmatchable by structure after standardization; it is not in the Repurposing Hub join either);
  - adagrasib (KRAS G12C);
  - pemigatinib (FGFR; infigratinib and erdafitinib are matched);
  - trimetrexate (DHFR);
  - fostamatinib (SYK prodrug);
  - belzutifan (HIF-2α; JUMP has only a 2,2-difluoro analog).
- **Investigational oncology agents:**
  - the newer RAS(ON) and SOS1 agents RMC-6236 and BI-3406 (the RAS class is therefore essentially absent);
  - 9-ING-41 (GSK3β);
  - simotinib (EGFR);
  - darinaparsin (organoarsenical);
  - "Tubulin inhibitor 6".
- **Other:**
  - talc, sacubitril/valsartan and verteporfin (all three unparseable or not single small molecules);
  - pasireotide (peptide);
  - megestrol, digitoxin, minodronic acid, pentoxifylline.

## Caveats

1. **Dose mismatch is the main design issue.**
   - Tahoe tests 0.05, 0.5 and 5 µM.
   - JUMP's 10 µM is 2× above Tahoe's top dose.
   - Source_7's 0.625 µM is close to Tahoe's middle dose (0.5 µM). That suits a matched-dose comparison, but it is a single lab with about 5 replicate wells per compound.
   - 155 matched drugs exist *only* in source_7.
   - The 186 drugs at both doses give a natural low/high-dose pairing.
2. **The concentrations are nominal.** They come from the JUMP Hub documentation, not per-well records, which JUMP does not publish.
3. **Cell line mismatch.** JUMP is U2OS only; none of Tahoe's 50 lines is U2OS (checked in `cell_line_metadata`). Any morphology↔transcriptome comparison is therefore cross-cell-line at the drug level.
4. **Stereo is not resolved** for the 155 tier-1b drugs or for the chiral tier-2b drugs. Enantiomers and epimers collapse:
   - doxorubicin and epirubicin → `JCP2022_063287`;
   - quinidine and quinine;
   - tretinoin and isotretinoin;
   - levobupivacaine and bupivacaine;
   - (R)-verapamil and verapamil.

   We cannot tell which stereoisomer JUMP plated.
5. **Tier 2b relies on structure alone.** The match requires an identical formula plus an identical heavy-atom graph.
   - Only dexamethasone is independently confirmed, by JUMP's control name. Spot checks of the others look like the same drugs drawn as odd tautomers.
   - The test is weakest for very small molecules. Lactate's JUMP match `CC(=O)C(O)O` is equally methylglyoxal hydrate.
   - A ring double-bond positional isomer could also pass in principle.
   - All 42 rows carry a REVIEW flag in `overlap.csv`.
6. **Tahoe duplicates:** 25 groups of Tahoe entries share one JUMP compound, mostly salt and free-base pairs such as crizotinib ×3 and docetaxel ±trihydrate. Use the unique-compound counts for anything statistical.

## Suggested next step

Proceed to planning. The best-supported first analysis is the **175 unique compounds (186 Tahoe drugs) measured at both JUMP doses**: they cover 15 MoA classes with ≥3 compounds each, including MAPK/MEK/RAF, MTOR, PI3K, EGFR, TK and DNA damage. Pair source_7 (0.625 µM) with Tahoe's 0.5 µM, and the 10 µM sources with Tahoe's 5 µM. Adding the LINCS Cell Painting set (cpg0004, A549, 6 doses) remains an option to fix the dose and replicate limits, but overlap is not the limiting factor.

## Reproduce

```bash
py -3.13 overlap.py   # downloads ~17 MB of metadata into ./data on first run; about 6 min first run, about 30 s after (cached)
```

Outputs: `overlap.csv`, `data/summary.json`, `data/moa_coverage.csv`. The script stops with an error if PubChem is unreachable, rather than silently skipping the cross-check.
