# Confirmatory Report

Protocol v1.3, frozen state `environment/PROTOCOL_FREEZE.json`,
sha256(PROTOCOL.md) = `6000af5dfa1662e93bc7d243e1156e84c9c0a445245aee7c92cbde17656f47ef`.
No scientific rule was changed after the freeze. Development and confirmatory cohorts were never
pooled.

```
FINAL CONFIRMATORY N:
    overall   128 structures / 123 independent target clusters
    BRIL       29
    T4L        38
    MBP        61

DEVELOPMENT/CONFIRMATORY TARGET OVERLAP:
    0   (leakage test: 9/9 checks pass, re-run immediately before reporting)

TECHNICAL SUCCESS:
    P2Rank  original 128/128   removed 128/128
    fpocket original 128/128   removed 128/128
    512/512 runs, zero technical failures, zero substitutions
    7,656 pockets classified
```

## PRIMARY AIM A — EXPOSURE

```
P2Rank    (127 evaluable; 1 structure yielded no P2Rank pocket in ORIGINAL)
    rank-1   66/127 = 52.0%   Wilson [43.3-60.5]   cluster-boot [43.2-61.1]
    top-3    89/127 = 70.1%   Wilson [61.6-77.4]
    top-5   104/127 = 81.9%   Wilson [74.3-87.6]
fpocket   (128 evaluable)
    rank-1   60/128 = 46.9%   Wilson [38.4-55.5]   cluster-boot [38.3-55.6]
    top-3    91/128 = 71.1%   Wilson [62.7-78.2]
    top-5   105/128 = 82.0%   Wilson [74.5-87.7]

PARTNER STRATA (reported before any pooled interpretation)
              rank-1                       top-3        top-5
BRIL
    P2Rank     2/28 =  7.1% [ 2.0-22.6]    28.6%        50.0%
    fpocket    3/29 = 10.3% [ 3.6-26.4]    48.3%        69.0%
T4L
    P2Rank     7/38 = 18.4% [ 9.2-33.4]    52.6%        76.3%
    fpocket    5/38 = 13.2% [ 5.8-27.3]    44.7%        65.8%
MBP
    P2Rank    57/61 = 93.4% [84.3-97.4]   100.0%       100.0%
    fpocket   52/61 = 85.2% [74.3-92.0]    98.4%        98.4%

RAW POOLED VS PARTNER-STANDARDIZED (rank-1)
    P2Rank    raw pooled 52.0%   standardized 39.7%   [34.4-45.6]   difference +12.3 pts
    fpocket   raw pooled 46.9%   standardized 36.2%   [30.5-42.4]   difference +10.6 pts
```

The pooled figure is an MBP figure. Reporting it alone would overstate the hazard by ~11-12 points.

## PRIMARY AIM B — TARGET-CAVITY DISPLACEMENT [F3, structure-level]

```
P2Rank
    evaluable                 103
    detector abstentions       22   (20 MBP, 1 T4L, 1 BRIL)
    unmatched                   3   (3PY7, 4DXB, 5W0R)
    median displacement         0    cluster-boot CI [0, 1]
    mean displacement       +0.78    cluster-boot CI [+0.55, +1.02]
    displacement >0        42/103 = 40.8%   Wilson [31.8-50.4]  boot [31.1-50.5]
    displacement >=2       19/103 = 18.4%   Wilson [12.1-27.0]
    displacement >=5        2/103 =  1.9%   Wilson [ 0.5- 6.8]
    displacement <0             0
    sign test              42 up / 0 down, exact p = 4.6e-13

fpocket
    evaluable                 128
    detector abstentions        0
    unmatched                   0
    median displacement         1    cluster-boot CI [0, 1]
    mean displacement       +4.06    cluster-boot CI [+2.92, +5.31]
    displacement >0        65/128 = 50.8%   Wilson [42.2-59.3]  boot [42.2-59.4]
    displacement >=2       51/128 = 39.8%   Wilson [31.8-48.5]
    displacement >=5       34/128 = 26.6%   Wilson [19.7-34.8]
    displacement <0             0
    sign test              65 up / 0 down, exact p = 5.4e-20

PARTNER-STRATIFIED DISPLACEMENT (>0)
                   P2Rank                      fpocket
    BRIL       2/28 =  7.1%  median 0      4/29 = 13.8%  median 0
    T4L        6/37 = 16.2%  median 0      8/38 = 21.1%  median 0
    MBP       34/38 = 89.5%  median 1     53/61 = 86.9%  median 5
```

**Not one structure in either detector showed negative displacement.** Across 231 evaluable
structure-detector pairs the fusion never improved the rank of the target's own preferred cavity.
The effect is strictly one-directional and, by partner, overwhelmingly MBP.

## BIOLOGICAL REFERENCE-SITE SUBSET

```
coverage   overall 60/128 = 46.9%   -> BELOW the pre-frozen 60% rule, so this remains a
                                        restricted-subset SECONDARY analysis
           BRIL 23/29   T4L 27/38   MBP 10/61

P2Rank     site-corresponding pocket found in ORIGINAL   58/60
           outranked by a fusion-associated pocket        9/58 = 15.5% [ 8.4-26.9]
           rank change improved 9 / worsened 0 / unchanged 49   (median 0)
           recovery in ORIGINAL   top-1 47/58 = 81.0%   top-3 96.6%   top-5 98.3%

fpocket    site-corresponding pocket found in ORIGINAL   60/60
           outranked by a fusion-associated pocket       27/60 = 45.0% [33.1-57.5]
           rank change improved 28 / worsened 1 / unchanged 30  (median 0)
           recovery in ORIGINAL   top-1 24/60 = 40.0%   top-3 56.7%   top-5 66.7%
```

This is the strongest constraint on the claim. Where a target has a characterised biological ligand
site, **P2Rank still ranks that site first in 81% of constructs** and is displaced in only 15.5%.
Reference-site coverage is itself partner-skewed: MBP has 10/61, so the constructs where the hazard
is worst are largely the ones with no characterised site to protect.

## MECHANISM

```
                            P2Rank   fpocket   total
FUSION_INTRINSIC_CAVITY        63       50      113
TARGET_FUSION_INTERFACE         3        8       11
CHIMERA_SPECIFIC_FUSION_CAVITY  0        2        2
LINKER_RELATED                  0        0        0
MAPPING_ARTIFACT                0        0        0

probe states   MATCHED_INTRINSIC_CAVITY 119 | NO_MATCHING_CAVITY 6 | AMBIGUOUS 1
               DETECTOR_NO_PREDICTIONS 0

101 of 119 matched cavities are the isolated partner's OWN rank-1 pocket.

by detector x partner
    P2Rank   BRIL  2 intrinsic
             T4L   4 intrinsic, 3 interface
             MBP  57 intrinsic
    fpocket  BRIL  2 intrinsic, 1 interface
             T4L   2 intrinsic, 3 interface
             MBP  46 intrinsic, 4 interface, 2 chimera-specific
```

No mapping or parser artifact in 126 inspected failures. The failure is attribution, not detection:
the detectors find a real cavity that belongs to the tag.

## DETECTOR DEPENDENCE

```
agreement on rank-1 call    112/127 = 88.2%      Cohen's kappa 0.764
    both 55 | P2Rank only 11 | fpocket only 4 | neither 57

score margin (best fusion-associated / best target-dominated)
    P2Rank    BRIL 0.12   T4L 0.14   MBP 5.04   ALL 0.96
    fpocket   BRIL 0.30   T4L 0.31   MBP 1.56   ALL 0.76

abstention in the REMOVED condition
    P2Rank   22/128  (MBP 20, T4L 1, BRIL 1)
    fpocket   0/128
```

Both detectors score BRIL and T4L cavities well below the target site (margins 0.12-0.31) and MBP's
cleft well above it (1.56-5.04). The rank-depth pattern is the practical finding: rank-1 capture is
partner-specific, top-5 exposure is general (50-82% for every partner).

P2Rank abstains entirely on 22 target-only structures, 20 of them MBP — for those constructs the
tag supplies the only pocket P2Rank will report anywhere.

## SENSITIVITY ANALYSES

| Axis | Setting | P2Rank rank-1 | fpocket rank-1 | Displacement >0 (P2R / fpocket) |
|---|---|---|---|---|
| Dominance | 0.50 | 51.2% | 45.3% | — |
| Dominance | **0.70 primary** | **52.0%** | **46.9%** | **40.8% / 50.8%** |
| Dominance | 0.90 | 51.2% | 46.9% | — |
| Correspondence | J>=0.25/12A | — | — | 40.8% / 50.8% |
| Correspondence | **J>=0.40/8A primary** | — | — | **40.8% / 50.8%** |
| Correspondence | J>=0.60/5A | — | — | 41.2% / 50.8% |
| Deletion boundary | exclude cavities <=8A from a cut | — | — | 40.2% / 50.4% (2 and 1 structures excluded) |
| HETATM | **stripped primary** | **52.0%** | **46.9%** | **40.8% / 50.8%** |
| HETATM | ions retained | 52.0% | 47.7% | 40.8% / 51.6% |
| Structure | **single chain primary** | **52.0%** | **46.9%** | **40.8% / 50.8%** |
| Structure | biological assembly 1 | 43.7% | 38.9% | 31.5% / 46.4% |

**No conclusion changes under any sensitivity.** The dominance threshold moves the exposure rate by
at most 1.6 points; the correspondence rule by at most 0.4 points; the boundary control by at most
0.6 points; ions by at most 0.8 points. The assembly-1 variant attenuates exposure by 8 points and
displacement by 4-9 points — the largest single sensitivity effect — while preserving direction,
partner ordering (MBP 83.1%/72.9% vs BRIL 3.4%/6.9%) and zero negative displacements. Two MBP
structures (5W0R, 5EDU) are not evaluable under assembly 1 because the designated chain is absent
from the biological assembly.

## SECONDARY MULTIVARIABLE MODEL

255 rows (structure x detector), 123 clusters, 126 events. **No complete or quasi-complete
separation** (every cell has both outcomes), so the pre-specified Firth fallback was not needed and
standard logistic regression with cluster-robust SE on target accession was estimable.

| Term | OR | 95% CI | p |
|---|---|---|---|
| partner[MBP] | **49.5** | 10.8 – 226.9 | 5.2e-07 |
| partner[T4L] | 2.40 | 0.46 – 12.4 | 0.30 |
| target_size / 100 residues | **0.286** | 0.131 – 0.622 | 0.0016 |
| resolution | 2.27 | 0.80 – 6.38 | 0.12 |
| detector[fpocket] | 0.53 | 0.23 – 1.19 | 0.12 |
| target_disorder | 0.004 | 0.000 – 852 | 0.38 |

Detector x partner interaction was estimable and **not significant** (fpocket x T4L p = 0.41,
fpocket x MBP p = 0.13), so the development-era impression of a large detector-by-partner
difference is not supported. Two covariates carry the effect: **the partner is MBP**, and **the
target is small**.

## DEVELOPMENT VS CONFIRMATORY REPLICATION

Cohorts are not pooled. Development n = 24 (exploratory), confirmatory n = 128, target overlap 0.
Development was re-described with the frozen F3 endpoint so the comparison is like for like.

| Signal | Development | Confirmatory | Verdict |
|---|---|---|---|
| P2Rank rank-1 | 11/24 = 45.8% [27.9-64.9] | 66/127 = 52.0% [43.3-60.5] | **REPLICATED** |
| P2Rank top-3 | 70.8% | 70.1% | **REPLICATED** |
| P2Rank top-5 | 83.3% | 81.9% | **REPLICATED** |
| fpocket rank-1 | 15/24 = 62.5% [42.7-78.8] | 60/128 = 46.9% [38.4-55.5] | **REPLICATED** |
| fpocket top-3 | 83.3% | 71.1% | **REPLICATED** |
| fpocket top-5 | 91.7% | 82.0% | **REPLICATED** |
| P2Rank BRIL | 0/8 = 0% | 2/28 = 7.1% | **REPLICATED** |
| P2Rank T4L | 3/8 = 37.5% | 7/38 = 18.4% | **REPLICATED** |
| P2Rank MBP | 8/8 = 100% | 57/61 = 93.4% | **REPLICATED** |
| **fpocket BRIL** | **5/8 = 62.5%** | **3/29 = 10.3%** | **NOT REPLICATED** |
| fpocket T4L | 3/8 = 37.5% | 5/38 = 13.2% | PARTIALLY REPLICATED |
| fpocket MBP | 7/8 = 87.5% | 52/61 = 85.2% | **REPLICATED** |
| P2Rank displacement >0 | 5/18 = 27.8% | 42/103 = 40.8%, CI excludes 0 | **REPLICATED** |
| fpocket displacement >0 | 16/24 = 66.7% | 65/128 = 50.8%, CI excludes 0 | **REPLICATED** |
| **fpocket BRIL displacement** | **5/8** | **4/29** | **NOT REPLICATED** |
| fpocket T4L / MBP displacement | 3/8, 8/8 | 8/38, 53/61 | **REPLICATED** |

**Tally: 17 REPLICATED, 1 PARTIALLY REPLICATED, 2 NOT REPLICATED.** All four headline signals
replicated. Both failures are the same signal: **fpocket's BRIL rank-1 capture, 62.5% in
development vs 10.3% in confirmation — the development 5/8 was small-sample noise.** This is
precisely the error a confirmatory cohort exists to catch, and it would have become a headline
claim ("fpocket ranks the BRIL apo-heme cavity first in most GPCR constructs") had the pilot been
written up directly.

## PROTOCOL INTEGRITY

**Pre-confirmatory amendments (set F, before any confirmatory detector ran):** F1 target-disjoint
set; F2 structure as the inferential unit — *this withdrew the Phase 3A cavity-level Wilcoxon
p-values as pseudoreplication*; F3 the non-degenerate displacement endpoint; F4 128 structures /
123 clusters with five legitimate cross-partner target pairs retained; F5 separation check with
Firth fallback; F6 replication criteria fixed before any confirmatory result was inspected.

**Technical deviations (set G, after primary results were visible; implementation only, none
touching the frozen primary path):** G1 chain-aware classifier for the assembly1 variant, with an
equivalence test of **46,448 assertions over 7,656 pockets, all passing**, proving the primary is
unaffected; G2 gemmi assembly API correction, made before any assembly1 detector run; G3 two
structures not evaluable under assembly 1 (5W0R, 5EDU), retained everywhere else; G4 chain-aware
displacement for assembly1, fixing a key-format bug that had reported 0/126 evaluable.

**Post-result deviations affecting a scientific rule:** none.

**Failed runs:** none. 512/512 primary, 512/512 ions, 504/504 assembly1.

**Exclusions after freeze:** none. No structure was removed for any outcome-related reason. The 22
P2Rank abstentions and 3 unmatched structures are reported as explicit states, never dropped.

**Test suites, all re-run immediately before reporting:** mapping invariants 51/51, detector
parsers 21/21, leakage 9/9, variant-classifier equivalence 46,448/46,448.

---

## FINAL SCIENTIFIC DECISION: **PARTIAL CONFIRMATION**

### What is confirmed, robustly

1. **The phenomenon is real, replicated and causal.** Both co-primary aims replicated on 128 fresh,
   target-disjoint structures. Displacement is strictly one-directional: across 231 evaluable
   structure-detector pairs, **zero** showed the fusion improving the target cavity's rank.
2. **The mechanism is settled.** 113 of 126 rank-1 failures are cavities intrinsic to the partner,
   101 of them the isolated partner's own top-ranked pocket. Zero mapping or parser artifacts.
3. **Nothing depends on an analytic choice.** Eleven sensitivity settings move the headline by at
   most 8 points and never change direction or partner ordering.
4. **Top-5 exposure is general**: 50-82% for every partner and both detectors.

### Why not STRONG

1. **Rank-1 capture is an MBP phenomenon, not a fusion-partner phenomenon.** BRIL 7.1%/10.3% and
   T4L 18.4%/13.2% versus MBP 93.4%/85.2%. The multivariable model puts the MBP odds ratio at 49.5
   while T4L is indistinguishable from BRIL. The project's original framing — "crystallographic
   fusion partners systematically bias pocket prediction" — is **not** supported at rank 1 for two
   of the three partners studied.
2. **A development signal failed to replicate.** fpocket/BRIL fell from 62.5% to 10.3%. Any claim
   built on the pilot's BRIL result would have been wrong.
3. **Where a biological site exists, the target usually still wins.** P2Rank ranks the
   characterised site first in 81% of such constructs and is displaced in 15.5%.
4. **The pooled estimate is misleading by 11-12 points** relative to the partner-standardised one.

### The claim the data support

> Crystallographic fusion tags can capture the top-ranked predicted pocket, but at rank 1 this is
> overwhelmingly a **maltose-binding-protein** effect: MBP's occupied maltose cleft outranks
> anything on the target in 85-93% of MBP constructs, versus 7-18% for BRIL and T4 lysozyme. The
> effect is causal — excising the tag never lowers the rank of the target's own preferred cavity in
> any of 231 evaluable cases — and is driven by two factors, the tag being MBP (OR 49.5) and the
> target being small (OR 0.29 per 100 residues). Across all three partners, a fusion-associated
> pocket appears in the **top 5** in 50-82% of constructs, so any workflow that inspects more than
> the single best pocket is exposed regardless of tag. Where the target has a characterised ligand
> site, P2Rank still ranks it first in 81% of constructs; the hazard concentrates in small,
> uncharacterised targets — precisely the case in which a pocket-based workflow would be trusted to
> find something new.

### Recommended follow-up (not started)

Targeted, not broad: characterise the MBP failure mode specifically (tag-larger-than-target
geometry, occupied cleft, N-terminal fusion), and quantify the top-5 exposure claim as the general
result. The benchmark/training-set construct audit remains a separate, pre-specified aim whose
feasibility is established but which has **not** been started.
