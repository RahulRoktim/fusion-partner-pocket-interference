# Phase 2G — Pilot Go/No-Go Report

> **CORRECTION (2026-09-15, Phase 3A).** The statement below that "the `E9` disorder filter
> excluded 1,223 of 1,971 candidates — more than any other criterion" is **incorrect** and should
> not be cited. Most of those flags were spurious: an empty-segment fallback assigned a disorder of
> 1.0 to 988 non-chimeric entities that had already failed criterion E1, and most of the remainder
> were redundant with E7/E8. E9's true marginal effect is **37 entities** (eligible pool 455 -> 492
> without it, +8.1%). See `DEVIATIONS.md` D2 and `results/e9_audit/`. No pilot result changes.


Protocol v1.1. Pilot is a feasibility and signal-detection exercise; no confirmatory inference is
drawn (amendment A4). All p-values below are **descriptive only**.

```
PILOT SAMPLE SIZE:            24  (8 BRIL, 8 T4L, 8 MBP)

ELIGIBLE POOL BEFORE SAMPLING:
    1,971 chimeric candidate entities enumerated from RCSB/SIFTS
      ->   455 eligible after protocol v1.1 §3.2/3.3 filters
             (BRIL 216, T4L 116, MBP 123)
      ->   153 one-per-(partner,target) collapsed units
             (BRIL 38, T4L 46, MBP 69 distinct target accessions)
      ->    24 drawn, seed "20260915-<PARTNER>", outcome-blind

FUSION PARTNERS:              BRIL (P0ABE7) 8 | T4L (P00720) 8 | MBP (P0AEX9) 8

P2RANK SUCCESSFUL RUNS:       48 / 48
FPOCKET SUCCESSFUL RUNS:      48 / 48
                              (96 / 96 total; zero technical failures)

P2RANK RANK-1 FUSION-ASSOCIATED:
    count / 24                11 / 24
    percentage                45.8 %
    95% CI                    27.9 – 64.9 %   (Wilson)

FPOCKET RANK-1 FUSION-ASSOCIATED:
    count / 24                15 / 24
    percentage                62.5 %
    95% CI                    42.7 – 78.8 %   (Wilson)

BRIL:
    P2Rank rate               0 / 8   =   0.0 %   [0.0 – 32.4]
    fpocket rate              5 / 8   =  62.5 %   [30.6 – 86.3]
T4L:
    P2Rank rate               3 / 8   =  37.5 %   [13.7 – 69.4]
    fpocket rate              3 / 8   =  37.5 %   [13.7 – 69.4]
MBP:
    P2Rank rate               8 / 8   = 100.0 %   [67.6 – 100.0]
    fpocket rate              7 / 8   =  87.5 %   [52.9 – 97.8]

TOP-3 RATES:                  P2Rank  17 / 24 = 70.8 %  [50.8 – 85.1]
                              fpocket 20 / 24 = 83.3 %  [64.1 – 93.3]

TOP-5 RATES:                  P2Rank  20 / 24 = 83.3 %  [64.1 – 93.3]
                              fpocket 22 / 24 = 91.7 %  [74.2 – 97.7]

TARGET REFERENCE-SITE COVERAGE:
                              13 / 24 = 54 %  -> BELOW the pre-frozen 60 % threshold.
                              Amendment A2 therefore applies: the cohort-level
                              fusion-capture / ranking endpoint REMAINS primary, and
                              biological target-site recovery is a restricted-subset
                              secondary endpoint. Coverage by partner:
                              BRIL 7/8, T4L 6/8, MBP 0/8.

PAIRED REMOVAL RESULTS:
    Outcome variable          "rank-1 pocket is TARGET_DOMINATED"  [C1 non-degenerate reading]
    P2Rank                    improved 5, worsened 0, unchanged 19   (descriptive p = 0.0625)
    fpocket                   improved 15, worsened 0, unchanged 9   (descriptive p = 0.0001)
    Best target-pocket rank   P2Rank  improved 5 / worsened 0 / unchanged 13 (of 18 comparable)
                              fpocket improved 15 / worsened 0 / unchanged 9 (of 24)
    Reference-site recovery   both detectors: lost 0, gained 1 (of 13 with a site)
    Empty predictions [C2]    P2Rank returned ZERO pockets after removal for 6 structures
                              (3N94, 3OAI, 5GPP, 5H7Q, 5JQE = MBP; 5YQR = T4L; all targets
                              <= 182 residues). Scored as NOT rescued — the conservative
                              direction, working against the hypothesis. fpocket: 0 such cases.
    Degenerate version        On "rank-1 is fusion-associated": 11->0 and 15->0 with the reverse
                              cell structurally fixed at 0, because the REMOVED condition
                              contains no fusion residues. Reported for audit only; this is a
                              restatement of the manipulation, not a result.

DELETION-BOUNDARY SENSITIVITY (8 A, pre-specified, not re-tuned):
    P2Rank                    5 apparent rank-1 rescues -> 4 survive the filter.
                              5K94 discounted (centroid 7.8 A from the cut).
                              1 of 86 REMOVED pockets lies within 8 A of a cut.
    fpocket                   15 apparent rank-1 rescues -> 15 survive (none discounted).
                              10 of 373 REMOVED pockets lie within 8 A of a cut.
    Reading                   The rescue signal is not an artifact of the in-silico excision.

DETECTOR AGREEMENT:
    structures scored by both 24
    agreement on rank-1 call  16 / 24 = 67 %
    Cohen's kappa             0.347  (fair)
    breakdown                 both 9 | P2Rank only 2 | fpocket only 6 | neither 7
    median rank of first      P2Rank 1.5 (present in 22/24) | fpocket 1.0 (present in 24/24)
      fusion-associated pocket
    score margin (best        P2Rank median 0.24 (n=16) | fpocket median 1.31 (n=24)
      fusion-assoc / best     by partner, P2Rank: BRIL 0.19, T4L 0.19, MBP 3.38
      target-dominated)       by partner, fpocket: BRIL 1.19, T4L 0.95, MBP 2.45

STRONGEST PARTNER EFFECT:     MBP — 8/8 (P2Rank) and 7/8 (fpocket) rank-1 capture, and the only
                              partner passing G2 in both detectors. Mechanistically the clearest:
                              in 6 of 8 MBP structures the top-ranked pocket contains two
                              deposited glucose units (GLC x2) — maltose, MBP's own natural
                              ligand — sitting inside the pocket.

MOST IMPORTANT EXAMPLE PDBs:
    5H7Q   MBP / MNDA PYD, 1.45 A   94-residue target fused to 370-residue MBP; both detectors
                                    rank MBP's maltose-occupied cleft first.
    7ZL9   BRIL / HCAR2, 2.70 A     fpocket ranks the BRIL apo-heme cavity first; P2Rank does not.
                                    The cleanest illustration of detector-dependent rank-1 capture.
    5K94   MBP, 2.10 A              interface cavity promoted from isolated-partner rank 19 to
                                    chimera rank 1; also the one boundary-discounted rescue.
    6WSK   T4L / CNRIP1, 1.55 A     P2Rank interface case, f_fusion 0.57.
    9JRT   BRIL, 3.28 A             negative control from inside the cohort: no fusion-associated
                                    pocket at any rank under P2Rank.

MECHANISTIC PATTERN:
    Of 26 rank-1 fusion-associated results (11 P2Rank + 15 fpocket, 17 distinct structures):
      21  FUSION_INTRINSIC_CAVITY   cavity reappears when the partner is run in isolation;
                                    in 19 of 21 it is the isolated partner's own rank-1 pocket
       4  TARGET_FUSION_INTERFACE   4Z35, 5JQE, 5K94, 6WSK
       1  AMBIGUOUS                 6A73 — probe limitation: P2Rank returned no pockets for the
                                    isolated partner, so overlap 0.00 is uninformative, not
                                    evidence of chimera-specificity
       0  MAPPING_ARTIFACT
       0  DETECTOR_ARTIFACT
       0  LINKER_RELATED
    The dominant failure is one of ATTRIBUTION, not of cavity detection: the detectors correctly
    find a real, high-quality, often occupied cavity that belongs to the crystallization tag.

PROTOCOL DEVIATIONS:
    A1–A5  PRE-OUTCOME   requested amendments: TAG deleted from both conditions; reference-site
                         coverage rule resolving §7.2/§7.4; deletion-boundary control at 8 A;
                         go/no-go language (10 % is a decision threshold, not a null); benchmark
                         audit adopted as a second aim with binding language constraints.
    B1–B3  PRE-OUTCOME   fpocket built from pinned source in WSL (no sudo, no Docker daemon);
                         gemmi for parsing; both detectors fed the identical prepared PDB.
    C1     POST-OUTCOME  paired fusion-association endpoint found structurally degenerate;
                         G3 re-read on "rank-1 is TARGET_DOMINATED", thresholds unchanged;
                         both versions reported.
    C2     POST-OUTCOME  empty REMOVED prediction sets scored as NOT rescued (conservative).
    C3     PRE-RUN       5IU7 smoke-test output seen during parser validation (not in cohort).
    C4     PRE-RUN       fpocket pocket-file off-by-one found and fixed by test D3 before any
                         pilot run; would have shifted every fpocket classification by one rank.
    No scientific threshold was altered after seeing outcomes. No observation was excluded on the
    basis of post-hoc QC.

TECHNICAL FAILURES:           None. 96 / 96 detector runs returned exit code 0 with parseable
                              output. The 6 empty P2Rank REMOVED prediction sets are successful
                              runs with genuinely zero predictions, handled under C2 — not
                              failures, and not silently replaced.

G1:                           PASS   (>=6 in one detector and >=4 in the other:
                                      P2Rank 11, fpocket 15)
G2:                           PASS   (MBP >=50 % in both detectors: 100 % and 87.5 %)
G3:                           PASS   (>=5 improved with <=1 reverse, both detectors:
                                      P2Rank 5/0, fpocket 15/0)

OVERALL:                      GO
```

## WHY

Three independent gates pass, and they fail in different ways if the effect were spurious.

1. **The effect is present in both detectors and is not threshold-sensitive.** Rank-1 capture is
   46 % (P2Rank) and 63 % (fpocket). P2Rank's rate is identical at 0.50, 0.70 and 0.90 dominance —
   its pockets are compositionally pure, so the headline number does not depend on where the
   dominance line is drawn. fpocket moves only from 54 % to 63 % across the same range.

2. **The paired removal is causal, directional, and survives the pre-specified artifact control.**
   Removal improved the rank-1 class in 5/24 (P2Rank) and 15/24 (fpocket) with **zero** reversals
   in either detector, and 4/5 and 15/15 of those rescues survive the 8 Å deletion-boundary filter.
   The conservative handling of empty prediction sets (C2) removed 6 potential P2Rank "rescues"
   that would otherwise have inflated this.

3. **The mechanism is specific and chemically verifiable, not a scoring quirk.** 21 of 26 failures
   are cavities intrinsic to the partner, and in six MBP structures the top-ranked pocket literally
   contains maltose. This is exactly the failure mode the study hypothesised, demonstrated from
   deposited coordinates without any externally curated site definition.

Two findings sharpen the full-study design rather than weakening the case:

- **Rank-1 capture is detector-dependent; top-5 exposure is not.** P2Rank ranks the BRIL cavity
  0/8 at rank 1 but finds it in the top 3 for 4/8, scoring it at ~0.19 of the target site. fpocket
  scores the same cavity at parity. Top-5 exposure is 83–92 % in both detectors. A pipeline that
  inspects more than the single best pocket is exposed whichever detector it uses.
- **MBP is the dominant hazard and is structurally distinct.** In MBP fusions the tag is routinely
  *larger* than the target (e.g. 5H7Q: 370 vs 94 residues) and carries an occupied native ligand
  cleft. MBP also has **zero** reference-site coverage — these targets usually have no bound
  biological ligand — which is precisely why the A2 coverage rule mattered.

## What justifies scaling — and what must change first

Evidence justifying the full cohort:

- Effect size is large enough that the planned n = 120–150 is comfortably powered: observed paired
  discordance is 15/24 (fpocket) and 5/24 (P2Rank) with zero reversals, against the protocol's
  power table which required p01 ≈ 0.15–0.25 for adequate power at n = 60–150.
- The annotation pipeline is validated end-to-end (51 mapping + 21 parser checks), runs at
  ~4.5 s/structure for P2Rank and ~1 s for fpocket, and produced zero technical failures.
- The eligible pool (455 entities, 153 distinct partner×target units) supports the target n without
  relaxing any filter.

Three design changes to settle **before** scaling, flagged now and not acted on:

1. ~~**The `E9` disorder filter is doing most of the excluding**~~ **[CORRECTED — see the note at the top of this file; E9's marginal effect is 37 entities, not 1,223.]** The original text read: — 1,223 of 1,971 candidates failed
   the ≤ 20 % segment-disorder rule, more than any other criterion. It is pre-specified and was
   applied as written, but it shapes the cohort heavily and deserves an explicit, documented
   decision (keep, relax, or stratify) rather than inheriting it by default.
2. **The isolated-partner probe needs a third outcome.** "Isolated run produced no pockets" is
   currently indistinguishable from "zero overlap," which is what made 6A73 ambiguous.
3. **MBP will dominate a stratified cohort unless capped.** With 69 eligible MBP targets vs 38
   BRIL and 46 T4L, and MBP carrying the strongest effect, the stratification plan needs a
   pre-specified per-partner cap so the headline number is not an MBP number wearing a cohort's
   clothes.

**The full cohort has not been started. No manuscript text has been written. No third detector has
been added.**
