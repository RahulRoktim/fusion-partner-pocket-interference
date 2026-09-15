# Phase 2F — Scientific QC Report

**Explanatory only.** No observation was removed from any endpoint on the basis of this QC. Every
structure that produced a successful detector run remains in every reported number, including the
structures that contradict the hypothesis.

Scope: all **26** rank-1 fusion-associated results (11 P2Rank + 15 fpocket) across **17** distinct
structures.

## Method

Two outcome-independent probes, neither of which requires any externally curated site definition:

1. **FUSION-ONLY condition (QC-only third condition).** The FUSION + LINKER residues from the same
   prepared coordinate file were run through both detectors *in isolation*. If the offending
   cavity reappears when the partner stands alone (≥ 50% of the chimera pocket's fusion residues
   recovered), the cavity is **intrinsic to the partner**; if not, it is chimera-specific.
2. **Deposited-group occupancy.** Any deposited heteroatom group within 5 Å of the pocket's
   residues is reported by chemical ID, including groups on the pre-declared artifact list.

## Result

| Mechanism | n | Partners |
|---|---|---|
| `FUSION_INTRINSIC_CAVITY` | **21** | MBP 13, BRIL 4, T4L 4 |
| `TARGET_FUSION_INTERFACE` | **4** | MBP 2, T4L 1, BRIL 1 |
| `AMBIGUOUS` | 1 | T4L 1 |
| `MAPPING_ARTIFACT` | **0** | — |
| `DETECTOR_ARTIFACT` | **0** | — |
| `LINKER_RELATED` | **0** | — |

**No rank-1 failure was caused by a mapping or parser artifact.** Every offending pocket's residues
resolved cleanly onto SIFTS entity positions; the unmapped fraction was below the 0.5 threshold in
all 26 cases.

## The dominant mechanism is the partner's own native ligand site

In 21 of 26 cases the cavity is intrinsic: it is found on the isolated partner, and in 19 of those
21 it is the isolated partner's **rank-1** pocket.

The strongest single piece of evidence is chemical rather than geometric. In six of the eight MBP
structures, the offending pocket contains **two deposited glucose units (`GLC` × 2)** — i.e.
**maltose, MBP's own natural ligand, sitting inside the pocket the detector ranked first**:

| PDB | Detector | f_fusion | Isolated-partner overlap (rank) | Groups in pocket |
|---|---|---|---|---|
| 3OAI | P2Rank, fpocket | 1.00 | 1.00 (rank 1) | GLC×2 |
| 5GPP | P2Rank, fpocket | 1.00 | 1.00 (rank 1) | GLC×2 / GLC×1 |
| 5H7Q | P2Rank, fpocket | 1.00 | 1.00 (rank 1) | GLC×2, ACT×1 |
| 7UAJ | P2Rank, fpocket | 1.00 | 1.00 (rank 1) | GLC×2 |
| 3N94 | P2Rank | 0.96 | 0.91 (rank 1) | GLC×2 |
| 5AZA | P2Rank, fpocket | 1.00 | 1.00 (rank 1) | GLC×2 |

This is not a subtle scoring artifact. The detectors are correctly finding a real, high-quality,
**occupied** carbohydrate-binding cleft — it simply belongs to the crystallization tag rather than
to the protein the user is studying. The failure is one of *attribution*, not of cavity detection.

## The four interface cases

`4Z35` (BRIL, fpocket), `5JQE` and `5K94` (MBP, fpocket), `6WSK` (T4L, P2Rank) produced rank-1
pockets drawing substantially on both sides of the junction (f_fusion 0.40–0.58). For 5JQE and
5K94 the isolated-partner overlap was 0.61–0.62 but only at isolated ranks 16 and 19 — i.e. the
constituent surface exists on the partner but is *promoted* to rank 1 by the presence of the
target. These are the genuinely chimera-specific artifacts and the most novel objects in the study.

## The one ambiguous case, and a caveat on the probe

`6A73` (T4L, P2Rank) was labelled `AMBIGUOUS`: the pocket is 77% fusion residues, but the
isolated-partner overlap is 0.00. **This is a limitation of the probe, not evidence of
chimera-specificity** — P2Rank returned **zero** pockets for isolated T4L in 6A73 (and for isolated
BRIL in 7ZL9), so there was nothing to overlap with. The same applies to any case where the
isolated run is empty. A refined probe for the full study should treat "isolated run produced no
pockets" as a distinct, uninformative outcome rather than as a zero overlap.

## Detector divergence worth carrying into the full study

The two detectors disagree systematically on BRIL, and the reason is interpretable:

| | P2Rank | fpocket |
|---|---|---|
| BRIL rank-1 failures | **0/8** | **5/8** |
| BRIL first fusion-associated pocket, median rank | 2–3 | 1 |
| BRIL score margin (best fusion / best target) | **0.19** | **1.19** |

P2Rank still *finds* the BRIL cavity — it is present in the top 3 for 4/8 BRIL structures — but its
learned scoring function assigns it roughly a fifth of the score of the GPCR orthosteric site, so it
does not reach rank 1. fpocket's geometric score puts the two on par. The same pattern holds for
T4L (P2Rank margin 0.19) and **reverses for MBP** (P2Rank margin 3.38, fpocket 2.45), where both
detectors rank the maltose cleft above anything on the target.

The practical reading: **rank-1 capture is detector-dependent, but presence in the top 5 is not**
(P2Rank 20/24, fpocket 22/24). A pipeline that inspects more than the single top pocket is exposed
regardless of detector choice.

## Structures that contradict the hypothesis

Reported explicitly, not buried:

- **All 8 BRIL structures under P2Rank** produced a target-dominated rank-1 pocket.
- `9JRT` (BRIL) produced no fusion-associated pocket at all under P2Rank, at any rank.
- `4N6H`, `4JKV`, `6RZ7`(P2Rank) and several T4L cases ranked the target site first.
- Under P2Rank, the paired removal changed the rank-1 class in only 5 of 24 structures; in 13 of
  24 the rank-1 pocket was already target-dominated and stayed so.

## Illustrative examples selected for the next phase

Not optimised as publication figures; chosen to span the mechanisms:

1. **5H7Q** (MBP–MNDA PYD, 1.45 Å) — cleanest case: 94-residue target fused to a 370-residue MBP;
   both detectors rank MBP's maltose-occupied cleft first.
2. **7ZL9** (BRIL–HCAR2, 2.70 Å) — fpocket ranks the BRIL apo-heme cavity first; P2Rank does not.
3. **5K94** (MBP, 2.10 Å) — interface cavity promoted from isolated rank 19 to chimera rank 1.
4. **6WSK** (T4L–CNRIP1, 1.55 Å) — P2Rank interface case, f_fusion 0.57.
5. **9JRT** (BRIL) — negative control from within the cohort: no fusion-associated pocket found.
