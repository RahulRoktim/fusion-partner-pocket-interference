# Phase 5 — Dataset Provenance / Training / Benchmark Construct Audit

> **CORRECTION (Phase 6B denominator audit).** Two numbers in the first version of this report
> were wrong and are corrected here: CHEN11's prevalence denominator (251 is the MEMBER count; the
> unique-PDB-entry count is 241), and the Tier-1 total (5,207 omitted FPTRAIN and used CHEN11's
> member count; the correct non-deduplicated sum is 5,419 and the deduplicated union is 5,304).
> The PLINDER partner breakdown is also restated in both units. No conclusion changes — the
> Tier-1 fusion-construct count remains 1. Canonical denominators:
> `results/PUBLICATION_DENOMINATORS.tsv`.

Descriptive audit. It does not touch the frozen confirmatory study; no confirmatory number was
recomputed and no confirmatory rule was consulted or changed.

Language states are kept distinct throughout: `PRESENT_IN_DATASET`, `PRESENT_IN_TRAINING_SET`,
`PRESENT_IN_DEVELOPMENT_SET`, `PRESENT_IN_VALIDATION_SET`, `PRESENT_IN_TEST_SET`,
`FUSION_POCKET_IS_GROUND_TRUTH_POSITIVE`, `FUSION_POCKET_IS_UNLABELLED_DECOY`,
`MEASURABLE_EFFECT_ON_BENCHMARK_SCORE`. The word *contamination* is not used as a conclusion.

A **standalone** MBP / T4L / BRIL structure is not a fusion construct and is counted separately
everywhere.

---

## P2RANK

Membership frozen from the published `.ds` files at `github.com/rdk/p2rank-datasets`
(master, accessed 2026-09-15), **not** from directory listings — the repository README warns the
directories hold more PDB files than the datasets define.

### CHEN11
```
role:                       TRAIN (distributed default P2Rank model)
version/source:             p2rank-datasets/chen11.ds, sha256 91caf227eb213f5a...
n:                          251 members / 241 unique PDB entries / 242 chains
entries resolved by SIFTS:  241/241 unique entries
fusion constructs:          0/241 unique PDB entries = 0.00%  [95% CI 0.00-1.57]
                            (0/251 members = 0.00% [0.00-1.51]; 10 entries contribute 2 chains)
                            chain-level: 0/242 unique benchmarked chains are chimeric
partners:                   none
fusion-labelled ground truth (FUSION_POCKET_IS_GROUND_TRUTH_POSITIVE):   0
unlabelled fusion decoys (FUSION_POCKET_IS_UNLABELLED_DECOY):            0
direct benchmark effect:    none measurable (no chimeric member exists)
also present, NOT constructs: 1 standalone partner protein (1LWG, T4 lysozyme)
```

### JOINED  (ASTEX 78 + B210 176 + BU48 93 + DT198 190)
```
role:                       DEVELOPMENT / VALIDATION (model selection, hyperparameters)
version/source:             p2rank-datasets/joined.ds, sha256 749f23f5da3bac82...
n:                          537 members / 537 entries / 190 chains (DT198 is chain-level)
entries resolved by SIFTS:  534/537
fusion constructs:          1/534 entries = 0.19%  [95% CI 0.03-1.05]   -> 1DUG (DT198)
partners:                   GST-Sj26 (P08515), TERMINAL_N topology
fusion-labelled ground truth:  1   -- see below, this is a clean CASE B
unlabelled fusion decoys:      0
direct benchmark effect:    P2Rank scores a TOP-1 SUCCESS on this member by predicting a pocket
                            that is 100% carrier residues; fpocket does not recover the
                            ground-truth site in any pocket at >=25% overlap
also present, NOT constructs: 13 standalone partner proteins (MBP 3, T4L 2, staph nuclease 4,
                            hen lysozyme 2, adenylate kinase 1, GST 1)
```

**1DUG in detail.** Title: *"Structure of the fibrinogen γ chain integrin binding and factor XIIIa
crosslinking sites obtained through carrier protein driven crystallization."* The benchmarked chain
A is 217 GST-carrier residues plus **10** fibrinogen target residues. The MOAD-2013 ground-truth
ligand is **GSH (glutathione)**, contacting the carrier with 72 contacts and the target with **0**.
The benchmark's positive label therefore *is* the carrier's own native glutathione site, and a
method is credited for finding it. This is the single cleanest instance of
`FUSION_POCKET_IS_GROUND_TRUTH_POSITIVE` in any Tier-1 dataset — and it is n = 1 of 534.

### COACH420
```
role:                       TEST
version/source:             p2rank-datasets/coach420.ds, sha256 4c33519448ae2df7...
n:                          420 members / 420 entries / 420 chains
entries resolved by SIFTS:  418/420
fusion constructs:          0/418 entries = 0.00%  [95% CI 0.00-0.91]
                            chain-level: 0/418 benchmarked chains are chimeric
partners:                   none
fusion-labelled ground truth:  0        unlabelled fusion decoys: 0
direct benchmark effect:    none measurable
other (biological, non-crystallization) chimeras: 2
also present, NOT constructs: 3 standalone partner proteins (148L T4L, 1STH staph nuclease,
                            3HPI MBP)
```

### HOLO4K
```
role:                       TEST
version/source:             p2rank-datasets/holo4k.ds, sha256 442a1abb07ee5ff9...
n:                          4,009 members / 4,009 entries (entry-level only, multi-chain)
entries resolved by SIFTS:  4,004/4,009
fusion constructs:          0/4,004 entries = 0.00%  [95% CI 0.00-0.10]
partners:                   none
fusion-labelled ground truth:  0        unlabelled fusion decoys: 0
direct benchmark effect:    none measurable
other (biological) chimeras: 6   (Gt/Gi chimeric Gα, myosin–dynamin, PXR–SRC1 linker constructs)
also present, NOT constructs: 27 standalone partner proteins (T4L 11, MBP 10, adenylate kinase 3,
                            staph nuclease 2, hen lysozyme 1)
```

**Overall for P2Rank's own data: 1 crystallization-fusion construct across 5,419 resolved entries
summed over its training, development and test sets — 5,304 entries after deduplication, since 115
entries appear in more than one Tier-1 dataset — and it sits in the development set.**

*Denominator note: per-dataset prevalence must use the per-dataset denominator; the union is only
"how many distinct entries were examined". See `results/PUBLICATION_DENOMINATORS.tsv`.*

---

## FPOCKET

### Version provenance (this corrects an earlier error in this project)

Phase 3A stated that fpocket "has no training set." **That was wrong and is withdrawn**
(`DEVIATIONS.md` H1).

```
version audited:            source tag 4.2.3, commit 4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066
scoring implementation:     src/pscoring.c, function score_pocket()
                            a FITTED LINEAR MODEL with hard-coded coefficients:
                              -0.03783394
                              + 0.48461469 * nas_norm
                              + 0.09093926 * as_density
                              + 0.0004155899 * convex_hull_volume
                              - 0.003995233 * surf_pol_vdw14
                              - 0.004072336 * surf_apol_vdw14
in-source documentation:    "The current scoring function has been determined using a logistic
                            regression based on an analysis of pocket descriptors"; adjacent
                            comments record PLS model variants and training-set performance
                            ("Train : 62/86 - 65/89")
descends from / modified:   DESCENDS from the original trained function AND was MODIFIED after
                            2009. git log on src/pscoring.c shows 2017-03-22 (3.0.4 import),
                            2017-03-23, 2017-03-27 "checking why scores are so off for explicit
                            pocket prediction", and 2017-03-27 "correcting mlhd normalized
                            descriptors by observed distribution in the full PDB" — i.e. a
                            recalibration of descriptor normalisation against the whole PDB.
                            4.2.3 provenance is therefore NOT inferable from the 2009 paper alone.
```

### Training / optimisation datasets

```
fpocket paper (Le Guilloux, Schmidtke & Tuffery 2009, BMC Bioinformatics 10:168):
   scoring function fitted by a PLS approach on a training set of 307 PROTEINS, derived from the
   protein test set of An et al. 2005 (PocketFinder evaluation set: 5,616 holo + 11,510 apo),
   filtered to <=50% sequence identity with manual validation. Evaluation was on the PocketPicker
   set (48), Cheng et al. (20) and Astex diverse (82) — validation, NOT training.
   The 307-protein membership list was not located in a downloadable form.

FPTRAIN (p2rank-datasets/fptrain.ds, sha256 af16cd4819b5947c..., 222 PDB entries):
   README verbatim: "FPTRAIN: dataset used by Fpocket for training its pocket scoring function".
   CAUTION: 222 != 307. FPTRAIN is a THIRD-PARTY characterisation by the P2Rank authors; it is
   NOT established to be identical to the set that produced the coefficients in fpocket 4.2.3.
   Reported as such, not as fpocket's own documented training set.
```

### Construct prevalence

```
FPTRAIN:  role TRAIN (fpocket scoring function, per p2rank-datasets)
          n: 222 entries, 222/222 resolved
          fusion constructs: 0/222 = 0.00%  [95% CI 0.00-1.70]
          partners: none
          fusion-labelled ground truth: 0     unlabelled fusion decoys: 0
          other (biological) chimeras: 1
          also present, NOT constructs: 3 standalone partner proteins
                                        (184L T4L, 1A2T staph nuclease, 1JJ0 hen lysozyme)
direct implications: none measurable. No crystallization-fusion construct was found in any
          dataset plausibly connected to fpocket's scoring provenance.
NOT counted as fpocket training data: the fpocket PREDICTIONS distributed in p2rank-datasets
          (README: "Fpocket, used version: v1.0 with default parameters") are the P2Rank authors
          running fpocket for comparison. Those are evaluation artefacts, not fpocket training.
```

---

## LIGYSIS

```
n:                     ~65,000 ligand binding sites across ~25,000 proteins, >100,000 PDBe
                       structures (authors' own figures, LIGYSIS-web README)
construct prevalence:  NOT MEASURED
handling:              NOT DETERMINED
limitation:            bulk membership is unavailable by design. LIGYSIS-web README, verbatim:
                       "We do not offer full LIGYSIS dataset download." Per-protein records are
                       browsable on the public server, but no canonical flat membership list was
                       obtainable. No prevalence figure is reported and none is inferred from the
                       dataset's design philosophy.
role:                  REFERENCE / evaluation substrate (Utgés & Barton 2024). Not a training set
                       for P2Rank or fpocket.
```

---

## sc-PDB

```
n:                     16,034 entries / 4,782 distinct proteins / 6,326 distinct ligands
                       (release v.2017, built on frozen PDB data 2016-11)
construct prevalence:  NOT MEASURED
training relevance:    documented training source for Kalasanty, PUResNet (5,020 structures),
                       DeepSite, DeepSurf, EquiPocket. NOT documented as a training or evaluation
                       set for P2Rank or fpocket, the two detectors in this study.
                       "structure is present in sc-PDB" is NOT equivalent to "this structure was
                       used by model X for training" — no per-model split was obtained, so no
                       per-model training claim is made for any of those predictors.
limitations:           the only distribution route located is a single 4,079,179,798-byte archive
                       (ressources/2016/scPDB.tar.gz). A measured transfer-rate test returned
                       ~11 kB/s, i.e. roughly four days for the full archive, and no lighter
                       membership index was found. Recorded as an objective, measured limitation.
```

---

## PLINDER

```
release:               2024-06 / v2
source:                https://storage.googleapis.com/plinder/2024-06/v2/splits/split.parquet
                       sha256 2959fb4b32f8c5cc3e12fc810b8905... (10,651,188 bytes)
scale:                 409,726 systems across 107,963 PDB entries
split sizes (systems): train 309,140 | removed 98,718 | test 1,036 | val 832
split sizes (entries): train 76,901 | removed 34,181 | test 1,020 | val 573
method:                reverse intersection. 1,405 PDB entries containing a recognised-partner
                       chimeric entity were enumerated from RCSB/SIFTS (12 partners) and
                       intersected with PLINDER membership. Exact for those partners; silent about
                       any partner outside the list.
```

### PDB ENTRY contains a construct vs PLINDER SYSTEM receptor IS the construct

```
systems whose PDB ENTRY contains a recognised construct:   1,894  (633 distinct entries)
   train 1,700 systems / 573 entries | removed 193 / 69 | val 1 / 1 | TEST 0 / 0

systems whose RECEPTOR CHAIN is itself the chimeric chain: 1,556  (512 distinct entries)
   train 1,385 systems / 454 entries | removed 170 / 61 | val 1 / 1 | TEST 0 / 0
   partners (systems / unique PDB entries; multi-partner entries counted under each partner):
             BRIL 662/161 | MBP 496/188 | T4L 211/96 | GFP 108/28 | SUMO-Smt3 31/20
             rubredoxin 28/15 | GST-Sj26 22/9 | thioredoxin 14/3 | ubiquitin 1/1
```

### Where the ground-truth ligand sits (chimeric-receptor systems, train/val/test)

```
train  (n = 1,385)   TARGET 939 | FUSION 342 | INTERFACE 7 | mixed 11
                     not on the chimeric entity 40 | unresolved 46
val    (n = 1)       TARGET 1
test   (n = 0)       -- no chimeric-receptor system in the test split
```

**`FUSION_POCKET_IS_GROUND_TRUTH_POSITIVE` in PLINDER train: 342 SYSTEMS across 158 UNIQUE PDB
ENTRIES** — both units must always be quoted together. Partner breakdown (systems / entries):
**MBP 330/148**, T4L 9/8, GST-Sj26 2/1, BRIL 1/1. The labelled ligand is
overwhelmingly **GLC** (754 ligand instances — glucose units, i.e. maltose in MBP's cleft), with
small numbers of SO4, GOL, GSH, BGC, FLC and others. **186 of the 342 pass PLINDER's own
`system_pass_validation_criteria`.**

**`FUSION_POCKET_IS_UNLABELLED_DECOY` in PLINDER train: 939 systems** where the receptor is a
fusion construct, the labelled site is on the target, and the partner's cavity is present but
unlabelled — e.g. `2rh1__1__1.A__1.I` (carazolol on β2AR, T4L present) and
`3eml__1__1.A__1.B` (ZM241385 on A2A, T4L present).

This is the same mechanism the confirmatory study measured, appearing independently in a training
corpus: MBP's occupied maltose cleft, defined as a protein–ligand ground-truth positive.

---

## DIRECT BENCHMARK EFFECT

```
structures evaluable:            1     (1DUG, JOINED/DT198 — the only recognised-partner chimeric
                                        member in any Tier-1 dataset)
target rank improved after removal:  NOT APPLICABLE. The ground truth for this member lies ON the
                                 carrier, so deleting the carrier deletes the label. The
                                 ORIGINAL-vs-FUSION-REMOVED contrast is undefined here and was
                                 not forced.
Top-1 outcome changed:           P2Rank recovers the ground-truth (carrier) site at RANK 1 with a
                                 pocket composed of 15 carrier residues and 0 target residues —
                                 a scored top-1 success for predicting the crystallization
                                 carrier's own site.
                                 fpocket: ground-truth site not recovered by any of its 15 pockets
                                 at >=25% overlap; its rank-1 pocket is also 100% carrier.
Top-3 / Top-(n+2) outcome:       P2Rank success (rank 1); fpocket failure.
fusion-labelled positives:       Tier 1: 1 (JOINED).      PLINDER train: 342.
unlabelled construct decoys:     Tier 1: 0.               PLINDER train: 939.
```

**Sample size is far too small for a metric claim about P2Rank or fpocket.** One member out of 534
cannot move an aggregate benchmark score measurably, and no aggregate score was recomputed. This is
a description of one instructive case, not a demonstration of performance bias.

---

## STATISTICS

Descriptive throughout. Numerators, denominators and Wilson intervals are given per dataset above.
No p-values are generated from these counts, and datasets of different role and provenance are not
pooled.

| Dataset | Role | n entries | constructs | prevalence [95% CI] |
|---|---|---|---|---|
| CHEN11 | TRAIN (P2Rank) | 241 unique entries (251 members) | 0 | 0.00% [0.00–1.57] |
| JOINED | DEVELOPMENT/VALIDATION | 534 | 1 | 0.19% [0.03–1.05] |
| COACH420 | TEST | 418 | 0 | 0.00% [0.00–0.91] |
| HOLO4K | TEST | 4,004 | 0 | 0.00% [0.00–0.10] |
| FPTRAIN | TRAIN (fpocket scoring) | 222 | 0 | 0.00% [0.00–1.70] |
| PLINDER train | TRAIN | 309,140 systems / 76,901 entries | 1,385 systems / 454 entries | — |
| PLINDER val | VALIDATION | 832 systems / 573 entries | 1 system / 1 entry | — |
| PLINDER test | TEST | 1,036 systems / 1,020 entries | **0 / 0** | — |
| LIGYSIS | REFERENCE | ~25,000 proteins | NOT MEASURED | — |
| sc-PDB | REFERENCE / de-facto TRAIN for other CNNs | 16,034 | NOT MEASURED | — |

---

## OVERALL EVIDENCE LEVEL: **3**

**LEVEL 3 — fusion constructs alter benchmark labels.** Reached on the label clause, and narrowly:

- **Level 1 (present in commonly used datasets):** met.
- **Level 2 (present specifically in training/development/evaluation sets used by pocket
  predictors):** met — 1 member of P2Rank's DEVELOPMENT/VALIDATION set, and 1,385
  chimeric-receptor systems in PLINDER's TRAIN split.
- **Level 3 (alter benchmark labels/ranks or measured evaluation outcomes):** met by two directly
  counted observations — 342 PLINDER-train systems in which a crystallization partner's cavity
  *is* the ground-truth positive, and one JOINED member (1DUG) where the same is true and where
  P2Rank is credited with a top-1 success for predicting it.
- **Level 4 NOT reached.** No model was retrained, no aggregate benchmark metric was recomputed,
  and no measured performance change is attributable to construct presence. With 0–1 constructs in
  P2Rank's and fpocket's own training and test sets, a material effect on their published scores is
  not merely unproven — it is implausible.

## DEFENSIBLE CLAIM

> Crystallographic fusion constructs are **essentially absent from the datasets that trained and
> benchmarked P2Rank and fpocket**: one construct across 5,419 resolved entries (5,304 deduplicated) spanning CHEN11
> (train), JOINED (development), COACH420 and HOLO4K (test), and FPTRAIN (fpocket's scoring
> provenance as characterised by the P2Rank authors). That single member, 1DUG in the development
> set, is nonetheless a clean illustration of the failure mode: its MOAD ground-truth ligand is
> glutathione bound entirely to the GST carrier, and P2Rank is scored as a top-1 success for
> predicting the carrier's site. In a modern, much larger protein–ligand corpus the same pattern
> appears at scale: PLINDER 2024-06/v2 contains 1,385 training systems whose receptor is a fusion
> construct, of which **342 systems / 158 unique PDB entries (330 systems of them MBP) define the crystallization
> partner's own cavity — overwhelmingly MBP's maltose cleft — as a ground-truth positive**, while
> a further 939 carry an unlabelled partner cavity alongside a target-located label. PLINDER's
> **test split contains no such system**.

## CLAIMS NOT SUPPORTED

- That P2Rank's or fpocket's published benchmark scores are affected by construct presence. They
  are not measurably affected, and the data argue against it.
- That any pocket-prediction benchmark is "contaminated". The word is not used as a conclusion.
- That any specific model was trained on fusion constructs via sc-PDB — no per-model split was
  obtained, so no such claim is made for Kalasanty, PUResNet, DeepSite, DeepSurf or EquiPocket.
- That LIGYSIS or sc-PDB do or do not contain constructs. Neither was measured.
- That the 342 PLINDER-train fusion-labelled systems have degraded any model. Nothing was trained.
- That PLINDER's labels are erroneous. Maltose genuinely binds MBP; the observation is about what
  the label *refers to*, not about whether it is chemically correct.

## DATASET-AUDIT LIMITATIONS

1. **sc-PDB and LIGYSIS not measured** — access limitations, both documented and measured rather
   than asserted (H3, H4).
2. **Partner recognition is by a fixed list of 12 UniProt accessions.** Constructs using partners
   outside that list are counted only as `OTHER_CHIMERA` (10 across Tier 1, all appearing to be
   biological rather than crystallization chimeras) and are not resolved further.
3. **2.9% of dataset polymer entities have no UniProt mapping** and cannot be assessed; 9 of 5,313
   entries did not resolve at all (likely obsolete).
4. **Entry-level vs chain-level.** HOLO4K, FPTRAIN and most of JOINED are entry-level, so a
   construct elsewhere in an entry would not imply the benchmarked chain was chimeric. This is
   reported separately wherever the distinction exists; for PLINDER it is the central distinction
   (1,894 entry-level hits vs 1,556 receptor-level).
5. **PLINDER ligand-site resolution** left 46 of 1,385 train systems unresolved and 40 with the
   ligand not on the chimeric entity.
6. **FPTRAIN's identity as fpocket's own training set is not established** (222 vs the paper's 307).
7. **The Tier-1 negative may be an artefact of dataset vintage and curation**, not of construct
   rarity: these sets were assembled from MOAD/PocketFinder-style curated complexes with
   drug-like-ligand filters, and predate or exclude the era in which GPCR–BRIL/T4L constructs
   became common. This is a hypothesis, not a measurement.

## PROTOCOL / TECHNICAL DEVIATIONS

`DEVIATIONS.md` amendment set H: **H1** the fpocket "no training set" correction; **H2** the JOINED
chain-identifier parsing gap (entry-level scan unaffected, not re-run); **H3** sc-PDB not audited,
with a measured download rate; **H4** LIGYSIS bulk membership unavailable by design; **H5** PLINDER
audited by reverse intersection rather than exhaustive scan.

## RECOMMENDATION: **INCLUDE AS SECONDARY ANALYSIS**

Not a major aim. The headline that would have justified one — "the benchmarks that trained and
evaluated these detectors are full of crystallization constructs" — **is false**, and this audit is
what establishes that it is false. Reporting the negative honestly is worth more than the
speculative positive was.

It earns a secondary section for three reasons. First, the Tier-1 negative is a genuine,
quantified result that pre-empts an obvious reviewer question about whether the main effect is a
benchmarking artefact — it is not, because the benchmarks contain essentially no constructs.
Second, 1DUG is a compact, fully worked illustration of `FUSION_POCKET_IS_GROUND_TRUTH_POSITIVE`,
including a scored top-1 success for predicting a carrier cavity. Third, the PLINDER-train finding
— 342 systems, 330 of them MBP, ligand almost always maltose — independently reproduces the main
study's central mechanistic result in a completely different corpus assembled by different people
for a different purpose, which is stronger corroboration than another benchmark of our own.

Not a supplement, because the PLINDER numbers speak directly to the paper's MBP conclusion. Not a
major aim, because nothing here demonstrates a performance consequence.
