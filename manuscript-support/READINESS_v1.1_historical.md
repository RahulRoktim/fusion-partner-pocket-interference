# Final Readiness Report

**State:** `v1.1-paper-final-freeze`. All pre-specified analyses are complete. No manuscript prose
has been written.

## What is complete

| Phase | Outcome |
|---|---|
| 0 — Literature / novelty audit | Verdict MODERATE; no equivalent prior study found |
| 1 — Protocol | Frozen v1.0, amended to v1.1 (pre-outcome), v1.2 (confirmatory), v1.3 (final) |
| 2 — Development pilot, n = 24 | GO; permanently exploratory, never pooled |
| 3A — Methods freeze | Degenerate endpoint replaced by target-cavity correspondence; E9 audited |
| 3B — Confirmatory, n = 128 | 512/512 runs; PARTIAL CONFIRMATION of a refined, partner-dependent hypothesis |
| 5 — Dataset provenance audit | Evidence LEVEL 3, narrowly, on the label clause; secondary analysis |
| 6 — Publication package | Tables 1–5, Figures 1–6, S1–S10, claims ledger, outline, freeze v1.0 |
| **J — E9 = 0.35 sensitivity** | **ROBUST.** The last outstanding pre-specified analysis |

## Sensitivity coverage: all six pre-specified axes executed

| Axis | Primary | Alternatives run | Conclusion changed? |
|---|---|---|---|
| Pocket dominance | 0.70 | 0.50, 0.90 | No |
| Target-cavity correspondence | J ≥ 0.40 / 8 Å | J ≥ 0.25/12 Å, J ≥ 0.60/5 Å, overlap-coef 0.50 | No |
| Deletion boundary | all predictions | exclude ≤ 8 Å from a cut | No |
| Structure representation | single designated chain | biological assembly 1 | No |
| HETATM handling | stripped | ions retained | No |
| **E9 segment disorder** | **0.20** | **0.35 (n = 141)** | **No** |

The E9 gap declared in `DEVIATIONS.md` I2 is now closed. Nothing pre-specified remains unexecuted.

## Standing constraints carried into writing

1. The original broad hypothesis is **not** fully confirmed. The confirmed hypothesis is the
   refined, partner-dependent one.
2. Primary manuscript values are the **E9 = 0.20, n = 128** confirmatory results. The E9 = 0.35
   analysis is sensitivity only and never replaces a primary value.
3. Partner strata are reported before any pooled number; every overall exposure estimate appears as
   both raw pooled and partner-standardised.
4. Development and confirmatory cohorts are never pooled for inference.
5. One structure = one inferential observation. Cavity-level statistics are descriptive only.
6. MBP identity and target size are **associated** predictors, not proven causal factors.
7. `DETECTOR_ABSTENTION` and `TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL` are reported as states,
   never assigned an imputed rank.

## Known limitations, all declared

- Three fusion partners; two detectors; no third detector.
- Association, not causation, for partner identity and target size.
- Biological reference-site coverage is below the 60% promotion rule in both cohorts
  (46.9% primary, 44.0% sensitivity), so Aim C stays secondary.
- sc-PDB and LIGYSIS were not measured in the dataset audit — access limitations, measured rather
  than asserted (`DEVIATIONS.md` H3, H4).
- Dataset-audit partner recognition is limited to 12 UniProt accessions.
- The Tier-1 dataset negative may reflect dataset vintage and curation rather than construct
  rarity; stated as a hypothesis, not a measurement.
- Mechanism (Aim D) was not recomputed for the E9 = 0.35 cohort; out of scope for that
  sensitivity.
- Figure 5 panels are Cα traces, not molecular-viewer cartoons; `FIG5_pymol.pml` reproduces them
  for publication rendering.

## Verification at freeze time

| Check | Result |
|---|---|
| Mapping invariants | 51 / 51 pass |
| Detector parsers | 21 / 21 pass |
| Development ↔ confirmatory leakage | 9 / 9 pass, target overlap 0 |
| Variant-classifier equivalence | 46,448 / 46,448 assertions pass |
| Phase 5 prose denominators | 11 / 11 assertions pass |
| Detector runs | primary 512/512, ions 512/512, assembly1 504/504, E9 new 64/64 — zero failures |

## READINESS: **READY FOR MANUSCRIPT**

No blocking issue remains. The next step is drafting prose against `CLAIMS_LEDGER.md` and
`MANUSCRIPT_OUTLINE.md`, which has not been started and is not authorised here.
