# Phase 0 — Literature and Novelty Audit

**Project:** Fusion-Tag Hazard — do crystallographic fusion partners bias automated binding-pocket prediction?
**Audit date:** 2026-09-15
**Status:** Complete. Verdict below. Phase 1 protocol frozen in `PROTOCOL.md`.

---

## 1. Search strategy

Databases queried: PubMed, Consensus (Semantic Scholar / Scopus / PubMed / arXiv), alphaXiv,
bioRxiv, and open web. Query families:

| # | Query family | Outcome |
|---|---|---|
| 1 | fusion protein x crystallization x binding-site prediction x artifact | 0 direct hits |
| 2 | P2Rank / fpocket x fusion protein / chimera / construct | 0 direct hits |
| 3 | MBP / T4L / BRIL x induced or artificial cavity | Construct-engineering papers only |
| 4 | pocket-detection artifacts, pocket-ranking errors | Benchmark papers; no construct-provenance analysis |
| 5 | expression / fusion tags x structural-bioinformatics caveats | Wet-lab caveats only |
| 6 | GPCR construct caveats x docking / virtual screening | Tutorial-level folklore only |
| 7 | benchmark-set curation (LIGYSIS, PLINDER, BioLiP, sc-PDB, HOLO4K, COACH420) | Ligand-artifact filtering only, never protein-artifact filtering |

**Important disambiguation.** In the current literature, "fusion protein binding-site prediction"
almost always means *oncogenic gene fusions* (BCR-ABL, EML4-ALK; e.g. FusionTarget, 2026), not
crystallographic fusion partners. This is a live collision risk for our searching and for the
eventual title/abstract. It also partly explains why the question appears unasked: the obvious
search terms are already occupied by a different field.

---

## 2. Closest prior work, and what each actually tested

### 2.1 Utges & Barton (2024), *J Cheminform* 16:126 — DOI 10.1186/s13321-024-00923-z

**Closest in methodology.** Largest benchmark of ligand-binding-site predictors to date: 13 methods
(fpocket, P2Rank, PRANK, DeepPocket, VN-EGNN, IF-SitePred, GrASP, PUResNet, Surfnet, Ligsite,
PocketFinder and others), 15 variants, 10 metrics, against a new curated reference set (LIGYSIS,
~30,000 proteins).

*What it tested:* recall and precision of predictors against biologically relevant ligand
interfaces; the cost of redundant pocket prediction; the benefit of stronger scoring schemes.
Reports that fpocket emits **8.4 predicted sites per defined reference site**. Proposes top-N+2
recall as a universal metric.

*What it did not test:* **where** the surplus pockets are, or whether the input structure is an
engineered construct. Construct provenance is not a variable anywhere in the study.

*Relevance:* this is the framework to adopt (ranking-aware metrics, top-N+2) and the gap to fill.
It establishes that surplus and mis-ranked pockets are a recognised problem, while leaving their
**origin** unexamined.

### 2.2 Bradford et al. (2021), *Chem Sci* — "Temperature artifacts in protein structures bias ligand-binding predictions"

**Closest in argument shape, and the template for framing.** Shows systematically that a
crystallographic *experimental* artifact (cryogenic vs room-temperature data collection) propagates
into downstream computational ligand-binding predictions — docking and rigorous binding free-energy
calculations — across five protein classes, using the T4 lysozyme L99A cavity as the workhorse.

*Why it matters:* proof that "crystallographic artifact leads to computational failure mode" is a
publishable genre in a strong journal, and a model for our rhetorical structure. It concerns
temperature, not construct composition. Note that it uses T4L's L99A cavity as a *model binding
site* — the field already knows T4L has excellent, ligandable cavities, which is exactly why
T4L-as-fusion-partner should be hazardous.

### 2.3 "Predicting binding sites from unbound versus bound protein structures" (2020), *Sci Rep* — DOI 10.1038/s41598-020-72906-7

Tests structure *state* (apo vs holo; 304 proteins with both) as an independent variable affecting
predictor performance. Same design logic as ours — structure provenance as the exposure — applied
to a different variable. Found no significant apo/holo difference for most methods, with a small
statistically significant advantage for fpocket on holo structures.

### 2.4 Method papers: fpocket and P2Rank

- fpocket: Le Guilloux, Schmidtke & Tuffery (2009), *BMC Bioinformatics* 10:168.
- P2Rank: Krivak & Hoksza (2018), *J Cheminform* — DOI 10.1186/s13321-018-0285-8.

Both define the detectors and their benchmarks (HOLO4K, COACH420). Neither analyses chimeric or
engineered constructs. HOLO4K is noted elsewhere to contain multi-chain assemblies that cause
distribution shift — an adjacent observation, but about oligomeric state, not constructs.

### 2.5 Construct-engineering literature — the source of the hazard, never its consequences

- Chun et al. (2012), *Structure* — "Fusion partner toolchest for the stabilization and
  crystallization of G protein-coupled receptors" (PMID 22681902). Defines the partner set: T4L,
  BRIL (thermostabilised apocytochrome b562RIL, M7W/H102I/R106L), rubredoxin, flavodoxin, xylanase,
  PGS.
- Thal et al. (2014), *Structure* — "Modified T4 lysozyme fusion proteins facilitate GPCR
  crystallogenesis" (PMID 25450769).
- Zou, Weis & Kobilka (2012), *PLoS ONE* — N-terminal T4L-beta2AR fusion.
- Miyagi et al. (2023), *Acta Cryst D* — anti-BRIL Fab crystallization chaperone; explicitly framed
  as enabling "structure-based drug design of membrane-protein drug targets."
- Smyth et al. (2003), *Protein Sci* — crystal structures of fusion proteins with large-affinity tags.
- "Design of an expression system to enhance MBP-mediated crystallization" (2017), *Sci Rep* 7:40991.

All of these describe *how to add* fusion partners. **None examines what the added domain does to
downstream computational analysis of the deposited coordinates.**

### 2.6 Dataset-curation literature — the cleanest statement of the gap

BioLiP, sc-PDB, PLINDER (Durairaj et al., bioRxiv 2024, DOI 10.1101/2024.07.17.603955) and LIGYSIS
all invest heavily in filtering **crystallization artifacts** — but exclusively on the *ligand*
side: buffers, cryoprotectants, additives, ions, single-atom entities. The field has a mature,
explicit concept of "this HETATM is an artifact of the experiment, not biology," and no
corresponding concept for "this **polypeptide** is an artifact of the experiment, not biology."

### 2.7 Folklore evidence (grey literature)

Docking tutorials and protocol write-ups routinely instruct users to delete the BRIL or T4 lysozyme
fusion before centring a docking grid box, and note that the fusion is "a separate stretch of
residues in the same chain, usually with a residue numbering jump." The hazard is therefore known
*qualitatively* to careful practitioners — and has never been quantified, localised, or shown to be
causal.

---

## 3. Is this exact question already answered?

**No.** Nothing found asks, let alone quantifies:

- how often a fusion partner or fusion interface captures the rank-1 predicted pocket;
- whether the effect replicates across independent detectors;
- whether removing fusion residues *causally* restores target-pocket ranking;
- whether standard pocket-prediction benchmark and training sets are themselves contaminated with
  fusion constructs.

## 4. What would still be novel

1. **Quantification.** First frequency estimate, with confidence intervals, of fusion-partner and
   interface capture of top-ranked pockets, stratified by partner and detector.
2. **Causality.** The paired in-silico removal experiment (same coordinates, same detector, fusion
   residues deleted) isolates the fusion partner as cause rather than correlate. No prior work does
   this for construct composition.
3. **Localisation.** Deterministic partition of high-ranking false pockets into fusion body vs
   target-fusion interface vs linker — distinguishing "the partner brought its own pocket" from
   "the chimera created a new one." The interface and linker classes are genuinely novel objects;
   nobody has reported an artificial-interface cavity as a prediction hazard.
4. **Benchmark contamination audit (recommended as a co-primary aim).** If HOLO4K, COACH420,
   sc-PDB, PLINDER or LIGYSIS contain fusion constructs, then predictors have been *trained and
   scored* on artifact geometry. This converts the paper from "a caveat for users" into "a defect
   in the field's evaluation substrate," and is the single biggest available upgrade to the study's
   importance.
5. **A deployable pre-flight screen.** A construct-aware check any pocket-based CADD pipeline can
   run before trusting a PDB entry.

## 5. Novelty verdict

### **MODERATE** — proceed.

**Upgradeable to STRONG** if the benchmark-contamination audit (4.4) is adopted as a co-primary aim.

Reasoning, stated plainly because it determines how the paper must be written:

- The *existence* of the effect is close to a priori certain, and a reviewer will say so. MBP has a
  maltose cleft. T4L has a lysozyme active-site cleft and the well-characterised L99A-type
  hydrophobic cavity. BRIL is **apo**-cytochrome b562 — its heme site is, by construction, a
  vacated cavity. These are three of the most ligandable small domains anyone could have chosen.
  "Pocket finder finds MBP's maltose site" is not a finding.
- Novelty therefore cannot rest on *whether* it happens. It must rest on **how often, at what rank,
  where, whether it is causal, and whether it has already contaminated the field's datasets.**
- Reframed that way, the work is clearly unpublished and defensible. Framed as "we discovered that
  fusion partners have pockets," it is rejectable on sight.

**Not ALREADY DONE, not WEAK.** Proceed to Phase 1 and the pilot.

---

## 6. Verification status of citations

Metadata confirmed via PubMed for: Utges & Barton 2024 (PMID 39529176, DOI
[10.1186/s13321-024-00923-z](https://doi.org/10.1186/s13321-024-00923-z)); Lazou et al. 2026
(PMID 42414572, DOI [10.1038/s42003-026-10596-z](https://doi.org/10.1038/s42003-026-10596-z),
consulted as adjacent context on pocket-prediction reliability). Remaining entries were located via
title/DOI-level search and **still require full-text reading before the manuscript's related-work
section is written** — in particular Bradford et al. 2021 and the PLINDER/LIGYSIS curation
sections, which are load-bearing for the gap claim in 2.6.
