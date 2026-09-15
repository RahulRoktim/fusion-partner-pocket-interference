# Phase 3A Report — Methods Freeze and Confirmatory Cohort Construction

No detector has been run on any confirmatory structure. No manuscript text has been written. No
third detector has been added.

```
DEVELOPMENT SET:
    24 structures, frozen, permanently exploratory
    8 BRIL / 8 T4L / 8 MBP, 24 distinct target accessions
    data_manifest/development_set.json  (status: DEVELOPMENT_EXPLORATORY_PERMANENT)
    These 24 observations may never be reused as confirmatory observations.

POCKET-CORRESPONDENCE METHOD:
  exact algorithm
    1. restrict every pocket to its TARGET residues only
    2. discard pockets with < 3 target residues (not matchable)
    3. similarity  = Jaccard index on target-residue sets
    4. geometry    = Euclidean distance between target-residue heavy-atom centroids
                     (target coordinates are byte-identical across conditions, test I5)
    5. assignment  = scipy linear_sum_assignment, one-to-one, maximising total Jaccard
    6. MATCHED if Jaccard >= 0.40 AND centroid distance <= 8.0 A
    7. states: MATCHED | UNMATCHED_ORIGINAL | NEW_AFTER_REMOVAL
    Detector rank and score are NEVER inputs to correspondence.

  thresholds        primary  J >= 0.40, centroid <= 8.0 A, min 3 target residues
  sensitivity       J >= 0.25 / 12 A ; J >= 0.40 / 8 A (primary) ; J >= 0.60 / 5 A ;
                    overlap-coefficient variant >= 0.50

  validation results
    Threshold selection used ONLY matching-quality metrics, computed and written out before any
    rescue quantity (results/correspondence_threshold_sweep.tsv).
      matched pair count is FLAT across Jaccard 0.20-0.60:
        P2Rank   86 pairs at every threshold from 0.20 to 0.60 (85 at 0.70)
        fpocket  367 -> 364 -> 363 across the same range
      matched pairs:  median Jaccard 1.000, median centroid distance 0.00 A,
                      median 8-9 shared target residues
      mean assignment ambiguity (2nd-best / best Jaccard): 0.068 P2Rank, 0.151 fpocket
    The threshold has almost no discretion to exercise; 0.40 sits mid-plateau.

  examples that matched pockets are the same cavity
    p2rank 3ODU rank 1 -> 1   J=1.000  d=0.00 A
      ORIGINAL target residues (40): 32,33,37,41,45,94,97,98,102,112,113,116,117,120,167,170,...
      REMOVED  target residues (40): identical; only-in-original [] ; only-in-removed []
    p2rank 3ODU ranks 2,3,4 -> 2,3,4 : J=1.000, d=0.00 A, residue sets identical
    fpocket 3N94 rank 1 -> 2   J=1.000  d=0.00 A, residue sets identical
    Interpretation: the cavity is re-detected with an identical residue set; only its RANK moves.

  state counts (primary rule, development set)
    P2Rank    MATCHED  86   UNMATCHED_ORIGINAL 16   NEW_AFTER_REMOVAL  0
    fpocket   MATCHED 364   UNMATCHED_ORIGINAL 52   NEW_AFTER_REMOVAL  9

PILOT MATCHED-TARGET-POCKET RESULT:
  P2Rank
    structures with >=1 matched target cavity          18 / 24
    structures with NO matched cavity                  6  (3N94 3OAI 5GPP 5H7Q 5JQE 5YQR)
    DISPLACEMENT (target's best cavity outranked by a fusion pocket in the construct)
                                                       5 / 18 = 27.8%  [12.5-50.9]
    PAIRED RANK CHANGE over 86 matched cavities
                                                       improved 42, worsened 1, unchanged 43
                                                       median 0, mean +0.81, range [-3, +9]
                                                       Wilcoxon p = 5.3e-08 (DESCRIPTIVE)
      restricted to cavities with a fusion pocket above them (n=43):
                                                       median +1, improved 42, worsened 1
  fpocket
    structures with >=1 matched target cavity          24 / 24
    DISPLACEMENT                                       15 / 24 = 62.5%  [42.7-78.8]
    PAIRED RANK CHANGE over 364 matched cavities
                                                       improved 313, worsened 32, unchanged 19
                                                       median +5, mean +5.60, range [-18, +46]
                                                       Wilcoxon p = 4.0e-47 (DESCRIPTIVE)
      restricted to cavities with a fusion pocket above them (n=327):
                                                       median +5, improved 309, worsened 17
  partner-stratified DISPLACEMENT
                        P2Rank              fpocket
    BRIL                0/8   =   0.0%      5/8  = 62.5%
    T4L                 2/7   =  28.6%      3/8  = 37.5%
    MBP                 3/3   = 100.0%      7/8  = 87.5%
    (P2Rank MBP denominator is 3, not 8: the other 5 MBP structures have no matched
     target cavity because P2Rank predicts nothing on those targets at all.)

  comparison with the degenerate pilot endpoint (audit only)
    P2Rank   degenerate 5 improved / 0 worsened   ->  non-degenerate displacement 5/18
    fpocket  degenerate 15 improved / 0 worsened  ->  non-degenerate displacement 15/24
    The effect survives the replacement of the degenerate outcome.

BIOLOGICAL-REFERENCE-SITE RESULT:
  coverage                     13 / 24 = 54%   (BRIL 7/8, T4L 6/8, MBP 0/8)
                               below the pre-frozen 60% rule -> endpoint C stays SECONDARY
  matched-site behaviour       a matched cavity overlapping the reference site was found in
                               13/13 structures, for BOTH detectors
  biological site outranked by a fusion pocket in the deposited construct
      P2Rank                   1 / 13 =  7.7%  [1.4-33.3]   (6A73 only, rank 3 -> 1)
      fpocket                  6 / 13 = 46.2%  [23.2-70.9]
  paired rank change of the biological site
      P2Rank                   improved 1, worsened 0, unchanged 12
      fpocket                  improved 6, worsened 0, unchanged 7
                               4Z35 4->2, 5ZBQ 13->12, 6A73 9->4, 6RZ4 24->18,
                               6RZ7 10->9, 7YXA 2->1
  THIS IS THE MOST IMPORTANT TEMPERING RESULT IN PHASE 3A. Where a genuine biological ligand site
  exists, P2Rank ranks it first despite the fusion in 12 of 13 cases. The hazard is concentrated in
  constructs WITHOUT a defined biological site — which is exactly the MBP stratum (0/8 coverage).

E9 AUDIT:
  implementation
    fusion_disorder = (|FUS_aligned| - |FUS_aligned AND modelled|) / |FUS_aligned|
    target_disorder = (|TGT_aligned| - |TGT_aligned AND modelled|) / |TGT_aligned|
    EXCLUDE if either > 0.20.  Segments are SIFTS-aligned label_seq positions; `modelled` is all
    entity positions minus RCSB UNOBSERVED_RESIDUE_XYZ for the representative chain.

  excluded counts — AND A CORRECTION TO THE PHASE 2G REPORT
    Phase 2G stated E9 excluded 1,223 of 1,971 candidates and called it the dominant population
    determinant. THAT WAS WRONG.
      * 988 of the flags are spurious: the fallback `... if seg else 1.0` assigns disorder 1.0 to
        an EMPTY segment, which happens for every non-chimeric single-UniProt entity — entities
        that had already failed E1 and were never candidates.
      * among the 983 genuinely chimeric entities, E9 flags 451 (45.9%), but only 37 are excluded
        by E9 ALONE; the rest also fail E7/E8 (too few modelled residues), the same underlying
        phenomenon measured differently.
      * TRUE MARGINAL IMPACT: eligible pool 455 -> 492 without E9 (+37 entities, +8.1%),
        adding 16 new independent targets (BRIL 8, T4L 2, MBP 6).
    by partner among genuinely chimeric entities: BRIL 266/492 flagged (20 sole),
    T4L 23/148 (7 sole), MBP 162/343 (10 sole).

  manual validation
    30 chimeric E9-excluded entities, outcome-blind random sample (seed E9AUDIT-20260915).
    Disorder recomputed directly from deposited coordinates, bypassing the RCSB feature.
      implementation agreement: 30 / 30 (100%)
      location of missing target residues: N-term 18%, INTERNAL 71%, C-term 11%
        -> genuine internal disorder, not unmodelled termini
      author-numbering gaps present in 20/30 entities but CANNOT be mistaken for disorder:
        the pipeline works in label_seq space, never author numbering
      construct residues absent from UniProt are captured as LINKER/TAG, not as disorder
      internal-insertion bookkeeping verified: target/fusion segments disjoint in all cases
    No numbering-gap, terminal-truncation, SIFTS-gap or entity-vs-author error was found.

  bias / selection implications
    E9 failure rate differs by topology (INTERNAL_INSERTION 25.5%, TERMINAL_N 61.7%,
    TERMINAL_C 59.6%) and by method (E9-failing set is 30% cryo-EM vs 4% in the passing set).
    E9 therefore does shift the population toward X-ray internal-insertion constructs. With a
    marginal impact of 37 entities this is a modest but real selection effect, and it is the
    reason completeness is added as a covariate rather than simply accepted.

  DECISION:  REVISE
    1. fix the empty-segment fallback (disorder not-evaluable when a segment is absent) — changes
       no eligibility decision, stops the audit trail misattributing exclusions
    2. KEEP the 0.20 primary threshold on measurement-validity grounds: unmodelled residues inside
       a segment create artificial surface of exactly the kind the A3 deletion-boundary control
       exists to guard against
    3. ADD completeness as a pre-specified 0.35 sensitivity threshold AND as a model covariate
    Made without running any comparison of the fusion effect under alternative E9 settings.

ZERO-PREDICTION AUDIT:
  per-structure findings (P2Rank, FUSION-REMOVED condition)
    pdb    partner  target  res   rc  runtime  input validity                 fpocket REMOVED
    3N94   MBP       94    1.80   0   ~4.3s   766 atoms, 0 breaks, 0 NaN     (pockets found)
    3OAI   MBP      121    2.10   0   ~4.3s   valid                          (pockets found)
    5GPP   MBP       86    2.00   0   4.32s   679 atoms, Rg 12.4 A, valid    1 pocket, score 0.269
    5H7Q   MBP       94    1.45   0   4.27s   766 atoms, Rg 13.5 A, valid    5 pockets, best 0.413
    5JQE   MBP      182    3.16   0   4.31s   1510 atoms, Rg 17.4 A, valid   12 pockets, best 0.330
    5YQR   T4L      108    2.40   0   4.22s   908 atoms, Rg 13.3 A, valid    6 pockets, best 0.401
    All six: return code 0, completed normally, no NaN coordinates, no chain breaks > 4.5 A,
    no residues lacking CA. Secondary structure of the removed target: 5GPP 91% helix,
    5H7Q 96% helix, 5JQE 73% helix, 5YQR 24% helix / 43% strand.

  interpretation
    GENUINE_DETECTOR_ABSTENTION in 6/6. No implementation or input problem.
    The pattern is size-dependent and monotone: the five smallest targets in the whole development
    set (86, 94, 94, 108, 121 residues) all yield zero P2Rank pockets, and P2Rank pocket count
    rises steadily with target size thereafter. fpocket always returns something, but with low
    scores (0.27-0.41).
    CONSEQUENCE FOR THE STUDY, not a nuisance: for five of eight MBP constructs, the fusion partner
    supplies the ONLY pocket P2Rank will report anywhere in the construct.
    No structure is removed on account of a zero prediction.

ISOLATED-PARTNER PROBE:
  revised categories
    MATCHED_INTRINSIC_CAVITY | NO_MATCHING_CAVITY | DETECTOR_NO_PREDICTIONS |
    REFERENCE_STRUCTURE_UNAVAILABLE | AMBIGUOUS
    Correspondence uses the SAME frozen rule (Jaccard >= 0.40, centroid <= 8 A) applied to
    FUSION residues. DETECTOR_NO_PREDICTIONS is never read as biological absence.

  revised state counts over the 26 rank-1 fusion-associated results
    MATCHED_INTRINSIC_CAVITY   24
    DETECTOR_NO_PREDICTIONS     1   (6A73 / P2Rank — isolated run returned nothing)
    NO_MATCHING_CAVITY          1   (6WSK / P2Rank — isolated run DID predict, none corresponds)

  revised mechanistic counts
    FUSION_INTRINSIC_CAVITY                  21   (unchanged)
    TARGET_FUSION_INTERFACE                   4   (unchanged)
    UNDETERMINED_DETECTOR_NO_PREDICTIONS      1   (was AMBIGUOUS)
    CHIMERA_SPECIFIC_FUSION_CAVITY            0
    MAPPING_ARTIFACT / LINKER_RELATED         0
    Of the 24 matched intrinsic cavities, 21 match the isolated partner's OWN rank-1 pocket.
    The revision changed exactly one label, and changed it from a conclusion to an abstention.

CONFIRMATORY FRESH POOL:
    BRIL independent targets        30   (195 entities, median 2.80 A;
                                          18 internal insertion, 10 N-term, 2 C-term)
    T4L  independent targets        38   ( 98 entities, median 2.70 A;
                                          29 internal insertion, 6 N-term, 2 C-term, 1 complex)
    MBP  independent targets        61   (113 entities, median 2.30 A;
                                          57 N-term, 4 complex)
    TOTAL                          129   independent (partner, target) units, 406 entities
    Excluded from the eligible 455: 24 development PDB entries + 25 entities sharing a
    (partner, target) pair with development. 1 confirmatory entity has a target that appears in
    development under a DIFFERENT partner — retained and flagged, not silently dropped.

PROPOSED CONFIRMATORY SAMPLE:
  exact selection strategy
    Use ALL 129 independent (partner, target) units. One structure per unit: highest resolution,
    ties broken by lowest PDB ID. No seed, no draw, no sampling step.
  whether all eligible fresh targets should be used
    YES. Using the entire eligible fresh population removes sampling discretion altogether, gives
    the largest defensible cohort, and makes the analysis un-rerunnable in the honest sense: there
    is no second draw available if the first result is unwelcome. Power is comfortable — the
    development displacement effects were 5/18 and 15/24 with 1 and 32 reversals respectively.

FINAL CONFIRMATORY PRIMARY ENDPOINT:
  Co-primary, both reported unconditionally:
  A  EXPOSURE   proportion of deposited constructs with a fusion-associated pocket at rank 1
                (headline), and in the top 3 and top 5; per detector, per partner.
  B1 DISPLACEMENT   proportion of structures in which the target's best matched cavity was
                outranked by >= 1 fusion-associated pocket in the deposited construct.
                Exact McNemar; risk difference with clustered-bootstrap CI.

FINAL SECONDARY ENDPOINTS:
  B2  paired rank change of matched target cavities (Wilcoxon signed-rank, median + CI)
  C   biological reference-site subset: B1/B2 on the matched cavity overlapping that site
      (promotable to primary only if confirmatory coverage reaches 60%, measured before any
      detector is run)
  D   mechanism: probe states and mechanism classes as revised above
  E   detector dependence: rank-1 and top-k rates, score margin, Cohen's kappa, by partner and
      rank depth
  Pre-specified sensitivities: dominance 0.50/0.70/0.90; correspondence 0.25/12A, 0.40/8A,
  0.60/5A, overlap-coefficient 0.50; deletion boundary 8 A; E9 disorder 0.20/0.35;
  single chain vs biological assembly 1; HETATM stripped vs ions retained.

STATISTICAL MODEL:
  * partner-stratified estimates are the PRIMARY presentation; no headline number hides them
  * any overall estimate reported TWICE: raw pooled AND partner-standardised (equal 1/3 weights).
    On development rates applied to confirmatory stratum sizes these differ by 12.5 points for
    P2Rank (58.3% pooled vs 45.8% standardised) and 4.5 points for fpocket.
  * mixed-effects logistic regression: outcome = fusion-associated rank-1; fixed effects =
    partner, topology, target size, resolution, method, detector, segment disorder;
    random intercept for target UniProt accession
  * paired endpoints: exact McNemar (binary), Wilcoxon signed-rank (continuous)
  * all CIs by clustered bootstrap on target accession; Wilson for simple proportions
  * Holm correction within the secondary family; co-primaries uncorrected against each other
  * effect sizes with CIs everywhere; p-values never reported alone
  * BALANCED SUBSAMPLING REJECTED — it discards valid data, cuts power in the most informative
    stratum, and makes the result depend on an arbitrary draw. Standardisation gives the same
    protection while keeping every observation. Chosen on design grounds, not by effect size.

BENCHMARK-AUDIT FEASIBILITY:
  (full table: results/benchmark_audit_feasibility.md — no contamination claim is made)
  HOLO4K    membership YES (holo4k.ds, 4,009 lines, fetched live)  IDs: entry only
            role TEST/BENCHMARK        SIFTS-testable YES at entry level
            P2Rank test set
  COACH420  membership YES (coach420.ds, 420 lines)                IDs: entry + CHAIN
            role TEST/BENCHMARK        SIFTS-testable YES at chain level (best matched)
            P2Rank test set
  CHEN11    membership YES (chen11.ds, 258 lines)                  IDs: entry + chain
            role TRAINING              SIFTS-testable YES
            documented training set of the distributed default P2Rank model
  LIGYSIS   membership PARTLY (pipeline/web public; bulk export unresolved)
            role REFERENCE             SIFTS-testable LIKELY, conditional on bulk export
            evaluation substrate for Utges & Barton 2024; trains nothing here
  sc-PDB    membership PARTLY (9,283 sites / 3,678 proteins; access terms need checking)
            role REFERENCE, de-facto TRAINING source for several CNN predictors
            SIFTS-testable YES if list obtainable
            documented training data for Kalasanty, PUResNet, DeepSite, DeepSurf, EquiPocket
  PLINDER   membership YES, with explicit TRAIN / VALIDATION / TEST splits
            role all three, separated by the authors
            SIFTS-testable YES, and split labels let each role be reported separately
            primarily co-folding/docking models, not the two detectors here
  Cautions: standalone partner entries are NOT chimeric constructs and our E1 criterion already
  separates them; entry-level and chain-level lists must be reported apart; fpocket is not an ML
  method and has no training set, so no training-contamination claim can attach to it.

NEW PROTOCOL DEVIATIONS:
  D1  PILOT-INFORMED  correspondence method replaces the degenerate paired outcome
  D2  PILOT-INFORMED  E9 REVISED; includes an explicit CORRECTION to the Phase 2G E9 claim
  D3  PILOT-INFORMED  isolated-partner probe states; one label changed (6A73)
  D4  PILOT-INFORMED  zero predictions are data, never failures, never removals
  D5  PILOT-INFORMED  development / confirmatory sets frozen and separated
  D6  PILOT-INFORMED  unequal strata modelled, not equalised; balanced subsampling rejected
  D7  PILOT-INFORMED  confirmatory cohort uses the whole eligible fresh population, no sampling
  Protocol v1.1 archived verbatim as PROTOCOL_v1.1_ARCHIVED.md; v1.0 as PROTOCOL_v1.0_ARCHIVED.md.

RECOMMENDATION:  PROCEED TO FULL CONFIRMATORY ANALYSIS
```

## Why PROCEED, and the one thing that should temper the eventual claim

The methodological objections raised against the pilot have been resolved rather than
side-stepped. The degenerate endpoint is replaced by a correspondence method whose thresholds
demonstrably carry no discretion — matched cavities have a median Jaccard of 1.000 and a median
centroid distance of 0.00 Å, and the matched-pair count is flat across the whole plausible
threshold range. The effect survives the replacement: displacement of the target's own best cavity
is 5/18 (P2Rank) and 15/24 (fpocket), with 42 vs 1 and 313 vs 32 rank improvements across matched
cavities.

The E9 audit changed my mind about E9 and forced a correction to the Phase 2G report. E9 is not
the dominant population determinant I claimed; its marginal effect is 37 entities. It is kept for
measurement-validity reasons and hedged with a sensitivity threshold and a covariate.

**The tempering result is endpoint C.** Where a genuine biological ligand site exists, P2Rank ranks
it first despite the fusion in 12 of 13 development structures; fpocket displaces it in 6 of 13.
The hazard is concentrated in constructs that have no defined biological site — overwhelmingly the
MBP stratum, which has 0/8 reference-site coverage and where P2Rank frequently predicts nothing on
the target at all. The defensible claim is therefore narrower than the raw exposure rate suggests,
and closer to:

> In crystallographic fusion constructs, automated pocket detectors frequently rank a cavity
> belonging to the crystallization tag above anything on the target — most severely when the target
> is small and has no characterised ligand site, which is precisely the situation in which a
> pocket-based CADD workflow would be relied upon to find something new.

That claim is worth testing on 129 fresh independent targets. The confirmatory run is defined,
frozen, and has not been executed.

## Paths

| Artifact | Path |
|---|---|
| Confirmatory protocol | `PROTOCOL.md` Part II (v1.2); archives `PROTOCOL_v1.0_ARCHIVED.md`, `PROTOCOL_v1.1_ARCHIVED.md` |
| Deviations | `DEVIATIONS.md` amendment set D |
| Development set | `data_manifest/development_set.json` / `.tsv` |
| Confirmatory pool | `data_manifest/confirmatory_candidate_pool.json` / `.tsv`, `confirmatory_pool_summary.json` |
| Correspondence | `scripts/08_pocket_correspondence.py`, `results/correspondence_threshold_sweep.tsv`, `results/matched_target_pockets.tsv`, `results/correspondence_summary.json` |
| Paired endpoint | `scripts/09_paired_endpoint.py`, `results/paired_endpoint_summary.txt` / `.json`, `results/paired_endpoint_per_structure.tsv` |
| E9 audit | `scripts/10_e9_audit.py`, `results/e9_audit/e9_distributions.tsv`, `e9_manual_validation.tsv` |
| Isolated-partner probe | `scripts/11_isolated_partner_probe.py`, `results/qc/isolated_probe_revised.tsv` / `.json` |
| Zero-prediction audit | `scripts/12_zero_prediction_audit.py`, `results/zero_prediction_audit.tsv` / `.json` |
| Set freezing | `scripts/13_freeze_sets.py` |
| Benchmark feasibility | `results/benchmark_audit_feasibility.md` |
| Corrected pilot report | `results/PHASE2G_GO_NOGO_REPORT.md` (correction note at top) |
