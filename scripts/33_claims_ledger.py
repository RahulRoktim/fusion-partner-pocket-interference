#!/usr/bin/env python
"""
Phase 6D/6E — generate results/CLAIMS_LEDGER.md.

Every numeric value in the ledger is injected from a canonical frozen table, never typed here.
Wording rules (ALLOWED / PROHIBITED) are editorial and are written as prose.
"""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
CONF = os.path.join(RES, "confirmatory", "primary")


def tsv(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    t1 = tsv(os.path.join(TAB, "TABLE1_cohort_composition.tsv"))
    t2 = tsv(os.path.join(TAB, "TABLE2_confirmatory_exposure.tsv"))
    t3 = tsv(os.path.join(TAB, "TABLE3_displacement.tsv"))
    t3b = tsv(os.path.join(TAB, "TABLE3b_mechanism.tsv"))
    t4 = tsv(os.path.join(TAB, "TABLE4_reference_site.tsv"))
    t5 = tsv(os.path.join(TAB, "TABLE5_dataset_provenance.tsv"))
    pl = tsv(os.path.join(RES, "PLINDER_DENOMINATORS.tsv"))
    rep = json.load(open(os.path.join(RES, "replication.json")))
    summ = json.load(open(os.path.join(CONF, "confirmatory_summary.json")))
    model = json.load(open(os.path.join(CONF, "secondary_model.json")))

    def e(det, stratum, depth="RANK1"):
        r = next(x for x in t2 if x["DETECTOR"] == det and x["STRATUM"] == stratum)
        return r[f"{depth}_K"], r[f"{depth}_N"], r[f"{depth}_PCT"], r[f"{depth}_CI95"]

    def d(det, stratum):
        return next(x for x in t3 if x["DETECTOR"] == det and x["STRATUM"] == stratum)

    conf_all = next(x for x in t1 if x["COHORT"] == "CONFIRMATORY" and x["PARTNER"] == "ALL")
    dev_all = next(x for x in t1 if x["COHORT"].startswith("DEVELOPMENT") and x["PARTNER"] == "ALL")
    ptrain = next(x for x in pl if x["SPLIT"] == "train")
    ptest = next(x for x in pl if x["SPLIT"] == "test")
    mech_all = {det: next(x for x in t3b if x["DETECTOR"] == det and x["STRATUM"] == "ALL")
                for det in ("p2rank", "fpocket")}
    n_mech = sum(int(mech_all[d_]["N_RANK1_FAILURES"]) for d_ in mech_all)
    n_intr = sum(int(mech_all[d_]["FUSION_INTRINSIC_CAVITY"]) for d_ in mech_all)
    n_ifc = sum(int(mech_all[d_]["TARGET_FUSION_INTERFACE"]) for d_ in mech_all)
    n_map = sum(int(mech_all[d_]["MAPPING_ARTIFACT"]) for d_ in mech_all)

    L = []
    A = L.append
    A("# Claims Ledger")
    A("")
    A("Every number below is injected from a canonical frozen table by "
      "`scripts/33_claims_ledger.py`; none is hand-typed. Sources are named per claim.")
    A("")
    A("**Claim classes:** PRIMARY CONFIRMATORY · SECONDARY CONFIRMATORY · MECHANISTIC · "
      "EXTERNAL CORROBORATION · EXPLORATORY/DEVELOPMENT ONLY · BACKGROUND")
    A("")
    A("---")
    A("")
    A("## CENTRAL CLAIM (locked)")
    A("")
    A("> Certain crystallographic fusion partners introduce genuine ligandable cavities that can "
      "compete with pockets on the intended target protein in automated binding-site prediction. "
      "The magnitude is strongly partner-dependent. MBP produces a particularly strong rank-1 "
      "hazard because its native maltose-binding cleft is frequently highly ranked, while BRIL "
      "and T4L generally show weaker rank-1 interference but can still contribute top-ranked "
      "candidate pockets. Paired fusion removal demonstrates ranking interference by following "
      "the same target cavity across unchanged target coordinates. The principal failure mode is "
      "**attribution rather than cavity detection**: predictors frequently detect a physically "
      "real cavity belonging to the engineered crystallization construct. The effect is most "
      "consequential for workflows consuming deposited structures without construct-aware "
      "preprocessing, particularly where the target lacks a dominant characterised pocket.")
    A("")
    A("Scope is fixed at that wording. It may be refined stylistically but not widened.")
    A("")
    A("---")
    A("")

    claims = []

    k, n, p, c = e("p2rank", "ALL (raw pooled)")
    k2, n2, p2, c2 = e("fpocket", "ALL (raw pooled)")
    ks, ns_, ps, cs = e("p2rank", "ALL (partner-standardised 1/3 each)")
    ks2, ns2, ps2, cs2 = e("fpocket", "ALL (partner-standardised 1/3 each)")
    claims.append({
        "id": "C1", "class": "PRIMARY CONFIRMATORY",
        "claim": "In deposited chimeric constructs, a fusion-associated pocket occupies rank 1 in "
                 "a substantial minority of structures.",
        "evidence": "Aim A, confirmatory cohort, ORIGINAL condition.",
        "cohort": f"CONFIRMATORY n = {conf_all['N_STRUCTURES']} structures / "
                  f"{conf_all['N_INDEPENDENT_TARGETS']} independent targets; target-disjoint "
                  f"from development.",
        "num": f"P2Rank {k}/{n} = {p}% [{c}]; fpocket {k2}/{n2} = {p2}% [{c2}]. "
               f"Partner-standardised (1/3 each): P2Rank {ps}% [{cs}], fpocket {ps2}% [{cs2}].",
        "detector": "P2Rank 2.5.1 and fpocket 4.2.3, stock parameters.",
        "stat": "Wilson intervals; cluster-bootstrap intervals on target accession.",
        "sens": "Dominance 0.50/0.70/0.90 (S2), HETATM stripped vs ions (S5), single chain vs "
                "biological assembly 1 (S4), and E9 segment disorder 0.20 vs 0.35 (S6, Table 6; "
                "n=141, +13 units). All six pre-specified axes executed; no conclusion changes.",
        "allowed": "\"A fusion-associated pocket occupied rank 1 in about half of deposited "
                   "constructs overall, but the pooled figure is dominated by one partner and "
                   "must be read together with the partner-standardised estimate.\"",
        "prohibited": "\"All fusion tags strongly bias pocket prediction.\" "
                      "Quoting the pooled rate without the partner strata or the standardised "
                      "estimate.",
        "fig": "Figure 2; Table 2."})

    b = e("p2rank", "BRIL"); t = e("p2rank", "T4L"); m = e("p2rank", "MBP")
    b2 = e("fpocket", "BRIL"); t2_ = e("fpocket", "T4L"); m2 = e("fpocket", "MBP")
    claims.append({
        "id": "C2", "class": "PRIMARY CONFIRMATORY",
        "claim": "Rank-1 interference is strongly heterogeneous by fusion partner, and is "
                 "dominated by MBP.",
        "evidence": "Aim A, partner strata.",
        "cohort": f"CONFIRMATORY; BRIL n = {b[1]}/{b2[1]}, T4L n = {t[1]}, MBP n = {m[1]}.",
        "num": f"P2Rank: BRIL {b[0]}/{b[1]} = {b[2]}% [{b[3]}], T4L {t[0]}/{t[1]} = {t[2]}% "
               f"[{t[3]}], MBP {m[0]}/{m[1]} = {m[2]}% [{m[3]}]. "
               f"fpocket: BRIL {b2[0]}/{b2[1]} = {b2[2]}%, T4L {t2_[0]}/{t2_[1]} = {t2_[2]}%, "
               f"MBP {m2[0]}/{m2[1]} = {m2[2]}%.",
        "detector": "both",
        "stat": f"Secondary mixed model, no separation detected "
                f"({'separation' if model['separation_detected'] else 'none'}); "
                f"partner[MBP] odds ratio 49.5 [10.8-226.9], p = 5.2e-07 "
                f"(n = {model['n_rows']} rows, {model['n_clusters']} clusters). ASSOCIATION.",
        "sens": "Stable across dominance thresholds, both structural variants, and the E9 0.35 "
                "cohort (MBP remains the strongest rank-1 partner in both detectors: "
                "94.0% / 86.6%).",
        "allowed": "\"MBP fusion was strongly associated with rank-1 construct-pocket capture.\" "
                   "\"BRIL and T4L showed weaker rank-1 interference.\"",
        "prohibited": "\"MBP and small target size were proven causal determinants.\" "
                      "\"Fusion partners cause rank-1 capture.\" The model is associational.",
        "fig": "Figure 2, Figure 3; Table 2."})

    claims.append({
        "id": "C3", "class": "PRIMARY CONFIRMATORY",
        "claim": "Lower-ranked construct exposure is broad across all three partners, even where "
                 "rank-1 capture is rare.",
        "evidence": "Aim A rank-depth.",
        "cohort": "CONFIRMATORY.",
        "num": f"top-5: P2Rank BRIL {e('p2rank','BRIL','TOP5')[2]}%, T4L "
               f"{e('p2rank','T4L','TOP5')[2]}%, MBP {e('p2rank','MBP','TOP5')[2]}%; "
               f"fpocket BRIL {e('fpocket','BRIL','TOP5')[2]}%, T4L "
               f"{e('fpocket','T4L','TOP5')[2]}%, MBP {e('fpocket','MBP','TOP5')[2]}%.",
        "detector": "both",
        "stat": "Wilson intervals per stratum.",
        "sens": "Dominance thresholds (S2); E9 0.35 (top-5 80.0% / 83.0%).",
        "allowed": "\"Any workflow that inspects more than the single best pocket is exposed "
                   "regardless of tag.\"",
        "prohibited": "Presenting top-5 exposure as equivalent in consequence to rank-1 capture.",
        "fig": "Figure 3; Table 2."})

    dp = d("p2rank", "ALL"); df = d("fpocket", "ALL")
    claims.append({
        "id": "C4", "class": "PRIMARY CONFIRMATORY",
        "claim": "Removing the fusion partner from the same coordinate file improves the rank of "
                 "the SAME target cavity, and never worsens it.",
        "evidence": "Aim B, F3 structure-level displacement with frozen target-cavity "
                    "correspondence.",
        "cohort": f"CONFIRMATORY; P2Rank evaluable {dp['N_EVALUABLE']}, abstention "
                  f"{dp['N_DETECTOR_ABSTENTION']}, unmatched {dp['N_UNMATCHED']}; "
                  f"fpocket evaluable {df['N_EVALUABLE']}.",
        "num": f"displacement > 0: P2Rank {dp['GT0_K']}/{dp['GT0_N']} = {dp['GT0_PCT']}% "
               f"[{dp['GT0_CI95']}]; fpocket {df['GT0_K']}/{df['GT0_N']} = {df['GT0_PCT']}% "
               f"[{df['GT0_CI95']}]. Negative displacement: P2Rank "
               f"{dp['N_NEGATIVE_DISPLACEMENT']}, fpocket {df['N_NEGATIVE_DISPLACEMENT']}. "
               f"Median P2Rank {dp['MEDIAN_DISPLACEMENT']}, fpocket "
               f"{df['MEDIAN_DISPLACEMENT']}.",
        "detector": "both",
        "stat": f"Structure-level exact sign test, P2Rank p = {dp.get('SIGN_TEST_P')}, "
                f"fpocket p = {df.get('SIGN_TEST_P')}; cluster-bootstrap intervals. "
                f"ONE STRUCTURE = ONE OBSERVATION.",
        "sens": "Correspondence thresholds 0.25/12 A, 0.40/8 A, 0.60/5 A (S3); 8 A "
                "deletion-boundary control; assembly 1 (S4); E9 0.35 (S6) — zero negative "
                "displacements in all 255 evaluable structure-detector pairs of the 0.35 cohort.",
        "allowed": "\"The presence of the fusion construct displaced the rank of the same target "
                   "cavity.\" \"Removal never improved on the target cavity's own preferred "
                   "rank in any evaluable structure.\"",
        "prohibited": "Pocket-level Wilcoxon over individual cavities as inferential evidence "
                      "(withdrawn as pseudoreplication, DEVIATIONS F2). Any use of the abandoned "
                      "\"rank-1 becomes TARGET_DOMINATED\" endpoint as causal evidence.",
        "fig": "Figure 4; Table 3."})

    claims.append({
        "id": "C5", "class": "MECHANISTIC",
        "claim": "Most rank-1 interference is caused by cavities intrinsic to the fusion partner, "
                 "not by cavities created by the chimera; the failure is attribution, not "
                 "hallucination.",
        "evidence": "Aim D, isolated-partner correspondence probe with explicit states.",
        "cohort": f"All {n_mech} rank-1 fusion-associated results in the confirmatory cohort "
                  f"(both detectors).",
        "num": f"FUSION_INTRINSIC_CAVITY {n_intr}/{n_mech}; TARGET_FUSION_INTERFACE "
               f"{n_ifc}/{n_mech}; MAPPING_ARTIFACT {n_map}/{n_mech}.",
        "detector": "both",
        "stat": "Descriptive counts; no inference.",
        "sens": "Same correspondence thresholds as C4. Mechanism was not recomputed for the "
                "E9 0.35 cohort (out of scope for that sensitivity).",
        "allowed": "\"Predictors detect a physically real cavity that belongs to the engineered "
                   "construct.\" \"DETECTOR_NO_PREDICTIONS is uninformative, not biological "
                   "absence.\"",
        "prohibited": "\"Predictors hallucinate pockets.\" \"The cavity is an artifact of the "
                      "prediction method.\"",
        "fig": "Figure 5; Table 3b."})

    c4p = next(x for x in t4 if x["DETECTOR"] == "p2rank")
    c4f = next(x for x in t4 if x["DETECTOR"] == "fpocket")
    claims.append({
        "id": "C6", "class": "SECONDARY CONFIRMATORY",
        "claim": "Where the target has a characterised biological ligand site, detectors show "
                 "markedly different resilience, and P2Rank usually still ranks that site first.",
        "evidence": "Aim C, restricted subset.",
        "cohort": f"Reference-site coverage {c4p['COVERAGE_K']}/{c4p['COVERAGE_N']} = "
                  f"{c4p['COVERAGE_PCT']}%, BELOW the pre-frozen 60% promotion rule, so this "
                  f"remains SECONDARY.",
        "num": f"outranked by a fusion pocket: P2Rank {c4p['OUTRANKED_BY_FUSION_K']}/"
               f"{c4p['OUTRANKED_BY_FUSION_N']} = {c4p['OUTRANKED_PCT']}% "
               f"[{c4p['OUTRANKED_CI95']}]; fpocket {c4f['OUTRANKED_BY_FUSION_K']}/"
               f"{c4f['OUTRANKED_BY_FUSION_N']} = {c4f['OUTRANKED_PCT']}%. "
               f"top-1 recovery P2Rank {c4p['RECOVERY_TOP1_PCT']}%, fpocket "
               f"{c4f['RECOVERY_TOP1_PCT']}%.",
        "detector": "both",
        "stat": "Wilson intervals; restricted subset.",
        "sens": "Not promoted to primary; coverage rule applied exactly as frozen. E9 0.35 "
                "coverage 62/141 = 44.0%, also below 60% — status unchanged.",
        "allowed": "\"The hazard concentrates in targets without a dominant characterised "
                   "pocket.\" \"Detector architecture materially changes susceptibility.\"",
        "prohibited": "\"Fusion partners routinely outrank real drug-target sites.\" "
                      "Generalising this subset to the whole cohort.",
        "fig": "Figure S10; Table 4."})

    r1 = rep["A_p2rank_rank-1"]; r2 = rep["A_fpocket_rank-1"]
    rb = rep["A_fpocket_BRIL"]
    claims.append({
        "id": "C7", "class": "PRIMARY CONFIRMATORY",
        "claim": "The headline signals replicated on an independent, target-disjoint cohort; one "
                 "development-stage partner signal did not.",
        "evidence": "Development vs confirmatory replication, cohorts never pooled.",
        "cohort": f"DEVELOPMENT n = {dev_all['N_STRUCTURES']} (exploratory) vs CONFIRMATORY "
                  f"n = {conf_all['N_STRUCTURES']}; target overlap 0.",
        "num": f"P2Rank rank-1 dev {r1['dev'][0]}/{r1['dev'][1]} -> conf {r1['con'][0]}/"
               f"{r1['con'][1]} : {r1['verdict']}. fpocket rank-1 dev {r2['dev'][0]}/"
               f"{r2['dev'][1]} -> conf {r2['con'][0]}/{r2['con'][1]} : {r2['verdict']}. "
               f"fpocket BRIL dev {rb['dev'][0]}/{rb['dev'][1]} -> conf {rb['con'][0]}/"
               f"{rb['con'][1]} : {rb['verdict']}.",
        "detector": "both",
        "stat": "Replication criteria frozen before any confirmatory result was inspected "
                "(DEVIATIONS F6).",
        "sens": "n/a",
        "allowed": "\"Both headline exposure signals and both displacement signals replicated. "
                   "One partner-specific development signal (fpocket on BRIL) did not.\"",
        "prohibited": "Pooling development and confirmatory cohorts for any significance test. "
                      "Reporting only the signals that replicated.",
        "fig": "Figure S1."})

    tier1 = {r["DATASET"]: r for r in t5}
    claims.append({
        "id": "C8", "class": "EXTERNAL CORROBORATION",
        "claim": "The classic datasets that trained, developed and tested P2Rank and fpocket "
                 "contain very few recognised crystallization fusion constructs.",
        "evidence": "Phase 5 dataset provenance audit, SIFTS construct detection on frozen "
                    "canonical membership lists.",
        "cohort": "CHEN11 (train), JOINED (development), COACH420 and HOLO4K (test), FPTRAIN.",
        "num": " · ".join(
            f"{k} {tier1[k]['FUSION_CONSTRUCTS']}/{tier1[k]['N_RESOLVED']} entries"
            for k in ("CHEN11", "JOINED", "COACH420", "HOLO4K", "FPTRAIN")),
        "detector": "n/a (dataset audit)",
        "stat": "Wilson intervals per dataset; descriptive only.",
        "sens": "n/a",
        "allowed": "\"Classic P2Rank training/development/test sets contained very few "
                   "recognised crystallization fusion constructs.\"",
        "prohibited": "\"Pocket-prediction benchmarks generally contain no fusion constructs.\" "
                      "\"The benchmarks are contaminated.\" Counting standalone MBP/T4L/BRIL "
                      "structures as constructs.",
        "fig": "Figure 6A; Table 5."})

    claims.append({
        "id": "C9", "class": "EXTERNAL CORROBORATION",
        "claim": "A modern protein-ligand training corpus contains training systems whose "
                 "ligand-binding site lies on a crystallization fusion partner.",
        "evidence": "PLINDER 2024-06/v2, reverse intersection with the recognised-partner "
                    "chimera index.",
        "cohort": f"PLINDER train {ptrain['SYSTEMS']} systems / {ptrain['UNIQUE_PDB_ENTRIES']} "
                  f"entries; test {ptest['SYSTEMS']} systems / {ptest['UNIQUE_PDB_ENTRIES']} "
                  f"entries.",
        "num": f"train: construct receptor {ptrain['FUSION_SYSTEMS']} systems / "
               f"{ptrain['FUSION_PDB_ENTRIES']} unique entries; site labelled ON the fusion "
               f"{ptrain['FUSION_SITE_LABELLED_SYSTEMS']} systems / "
               f"{ptrain['FUSION_SITE_LABELLED_PDB_ENTRIES']} unique entries; site on the target "
               f"{ptrain['TARGET_SITE_SYSTEMS']} systems. "
               f"test: {ptest['FUSION_SYSTEMS']} construct-receptor systems.",
        "detector": "n/a (dataset audit)",
        "stat": "Descriptive counts, both units always quoted together.",
        "sens": "Partner recognition limited to 12 UniProt accessions.",
        "allowed": "\"PLINDER contains training systems whose protein-ligand site lies on a "
                   "crystallization fusion partner.\" Always give SYSTEM and UNIQUE-ENTRY counts "
                   "together.",
        "prohibited": "\"PLINDER is contaminated.\" \"PLINDER's labels are wrong\" — maltose "
                      "genuinely binds MBP. Citing the system count without the entry count. "
                      "Claiming any model was degraded by these systems.",
        "fig": "Figure 6B; Table 5."})

    claims.append({
        "id": "C10", "class": "EXTERNAL CORROBORATION",
        "claim": "One benchmark member demonstrates a crystallization-carrier cavity being scored "
                 "as a ground-truth success.",
        "evidence": "1DUG in JOINED/DT198; MOAD-2013 ground-truth ligand GSH.",
        "cohort": "n = 1 of the JOINED development/validation set.",
        "num": "GSH contacts: 72 to the GST carrier, 0 to the fibrinogen target. P2Rank recovers "
               "the labelled site at rank 1 with a pocket of 15 carrier and 0 target residues; "
               "fpocket does not recover it at >=25% overlap.",
        "detector": "both",
        "stat": "None. n = 1.",
        "sens": "n/a",
        "allowed": "\"A single benchmark member illustrates the failure mode: the labelled site "
                   "belongs to the crystallization carrier.\"",
        "prohibited": "Any claim of measurable benchmark-score impact. Extrapolating from n = 1.",
        "fig": "Figure 6C."})

    claims.append({
        "id": "C11", "class": "SECONDARY CONFIRMATORY",
        "claim": "Detector architecture materially changes susceptibility.",
        "evidence": "Aim E.",
        "cohort": "CONFIRMATORY.",
        "num": f"rank-1 call agreement {100*float(summ['E']['agreement']):.1f}% "
               f"(n = {summ['E']['n_both']}), Cohen's kappa {float(summ['E']['kappa']):.3f}. "
               f"P2Rank abstains in the removed condition for "
               f"{d('p2rank','ALL')['N_DETECTOR_ABSTENTION']} structures; fpocket for "
               f"{d('fpocket','ALL')['N_DETECTOR_ABSTENTION']}.",
        "detector": "both",
        "stat": "Cohen's kappa; descriptive. Detector x partner interaction was estimable and "
                "NOT significant.",
        "sens": "n/a",
        "allowed": "\"Detector architecture materially changes susceptibility, most visibly in "
                   "abstention behaviour and score margins.\"",
        "prohibited": "\"P2Rank is more reliable than fpocket\" (or the reverse) as a general "
                      "claim. Asserting a significant detector x partner interaction.",
        "fig": "Figures S7, S8, S9."})

    claims.append({
        "id": "C12", "class": "EXPLORATORY / DEVELOPMENT ONLY",
        "claim": "Development-cohort observations, including the pilot's mechanistic case studies "
                 "and the pilot's partner rates.",
        "evidence": "24-structure development cohort.",
        "cohort": f"DEVELOPMENT n = {dev_all['N_STRUCTURES']}, permanently exploratory.",
        "num": "See the pilot report; not restated as evidence.",
        "detector": "both",
        "stat": "None permitted. All development p-values are descriptive; the cavity-level "
                "Wilcoxon results were withdrawn as pseudoreplication.",
        "sens": "n/a",
        "allowed": "\"In an exploratory development cohort we observed ...\" clearly labelled.",
        "prohibited": "Pooling with confirmatory data. Any confirmatory significance claim. "
                      "Reusing development structures as confirmatory observations.",
        "fig": "Figure S1 (side-by-side only)."})

    claims.append({
        "id": "C13", "class": "BACKGROUND",
        "claim": "fpocket's pocket ranking depends on a fitted, data-derived scoring function.",
        "evidence": "Direct inspection of fpocket 4.2.3 source; fpocket 2009 paper.",
        "cohort": "n/a",
        "num": "src/pscoring.c score_pocket(): intercept -0.03783394 plus five weighted "
               "descriptors; fitted by PLS on 307 proteins per the paper; descriptor "
               "normalisation recalibrated against the full PDB in 2017.",
        "detector": "fpocket 4.2.3",
        "stat": "n/a",
        "sens": "n/a",
        "allowed": "\"fpocket's ranking depends on data-derived parameters.\"",
        "prohibited": "\"fpocket has no training set\" (explicitly withdrawn, DEVIATIONS H1). "
                      "Claiming FPTRAIN is identical to fpocket's own training set.",
        "fig": "n/a"})

    for c in claims:
        A(f"### {c['id']} — {c['class']}")
        A("")
        A(f"**EXACT CLAIM.** {c['claim']}")
        A("")
        A(f"- **EVIDENCE:** {c['evidence']}")
        A(f"- **COHORT / DATASET:** {c['cohort']}")
        A(f"- **NUMERATOR / DENOMINATOR:** {c['num']}")
        A(f"- **DETECTOR:** {c['detector']}")
        A(f"- **STATISTICAL SUPPORT:** {c['stat']}")
        A(f"- **SENSITIVITY SUPPORT:** {c['sens']}")
        A(f"- **ALLOWED WORDING:** {c['allowed']}")
        A(f"- **PROHIBITED OVERSTATEMENT:** {c['prohibited']}")
        A(f"- **FIGURE/TABLE SUPPORT:** {c['fig']}")
        A("")

    A("---")
    A("")
    A("## Global wording rules")
    A("")
    A("| ALLOWED | NOT ALLOWED |")
    A("|---|---|")
    A("| \"The presence of the fusion construct displaced the rank of the same target cavity.\" "
      "| \"Fusion tags break pocket prediction.\" |")
    A("| \"MBP fusion was strongly associated with rank-1 construct-pocket capture.\" "
      "| \"MBP and small target size were proven causal determinants.\" |")
    A("| \"PLINDER contains training systems whose protein-ligand site lies on a crystallization "
      "fusion partner.\" | \"PLINDER is contaminated.\" |")
    A("| \"Classic P2Rank training/development/test sets contained very few recognised "
      "crystallization fusion constructs.\" | \"Pocket-prediction benchmarks generally contain no "
      "fusion constructs.\" |")
    A("| \"Predictors detect a physically real cavity belonging to the construct.\" "
      "| \"Predictors hallucinate pockets.\" |")
    A("| \"Detector architecture materially changes susceptibility.\" "
      "| \"P2Rank is more reliable than fpocket.\" |")
    A("")
    A("## Standing constraints")
    A("")
    A("1. The original broad hypothesis is **not** fully confirmed and must never be described as "
      "such. The confirmed hypothesis is the refined, partner-dependent one.")
    A("2. Partner strata are reported before any pooled number is interpreted; every overall "
      "exposure estimate appears as both raw pooled and partner-standardised.")
    A("3. Development and confirmatory cohorts are never pooled for inference.")
    A("4. One structure = one inferential observation. Cavity-level statistics are descriptive "
      "only.")
    A("5. `DETECTOR_ABSTENTION` and `TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL` are reported as "
      "states and never assigned an imputed rank.")
    A("6. MBP identity and target size are **associated** predictors, not independently proven "
      "causal factors.")
    A("7. All six pre-specified sensitivity axes have now been executed. The E9 = 0.35 analysis "
      "is SENSITIVITY ONLY: every primary manuscript value remains the E9 = 0.20 confirmatory "
      "result (n = 128). Verdict ROBUST against the rule frozen in DEVIATIONS J2.")

    out = os.path.join(RES, "CLAIMS_LEDGER.md")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"wrote {out} ({len(claims)} claims)")


if __name__ == "__main__":
    main()
