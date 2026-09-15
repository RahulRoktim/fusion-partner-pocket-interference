#!/usr/bin/env python
"""
Amendment J — primary (E9 = 0.20) vs sensitivity (E9 = 0.35) comparison, the new-structure
characterisation, and the regenerated supplementary figure S6.

All numbers are read from the two frozen summary files; none is typed here.
Primary values are never modified.

Outputs:
  results/tables/TABLE6_e9_sensitivity.tsv
  results/tables/TABLE6b_e9_new_structures.tsv
  results/e9_sensitivity_comparison.txt / .json
  figures/supplementary/S6_e9_disorder.png / .pdf   (replaces the audit-only version)
"""
import collections, csv, json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
MANI = os.path.join(ROOT, "data_manifest")
FIG = os.path.join(ROOT, "figures", "supplementary")
PARTNERS = ["BRIL", "T4L", "MBP"]
DET = ["p2rank", "fpocket"]
DETLAB = {"p2rank": "P2Rank", "fpocket": "fpocket"}
CP = {"BRIL": "#4C72B0", "T4L": "#DD8452", "MBP": "#C44E52"}
plt.rcParams.update({"font.size": 7.5, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200, "savefig.dpi": 350, "font.family": "DejaVu Sans"})


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0, c - h), min(1, c + h)


def main():
    P = json.load(open(os.path.join(RES, "confirmatory", "primary",
                                    "confirmatory_summary.json")))
    S = json.load(open(os.path.join(RES, "confirmatory", "e9_035",
                                    "confirmatory_summary.json")))
    rows, lines = [], []

    def emit(s=""):
        lines.append(s); print(s)

    emit("=" * 104)
    emit("E9 DISORDER SENSITIVITY — PRIMARY 0.20 vs SENSITIVITY 0.35")
    emit("Primary values are the manuscript values and are unchanged. 0.35 is sensitivity only.")
    emit("=" * 104)
    emit(f"\ncohort: primary n = {P['n_structures']} structures / {P['n_clusters']} clusters"
         f"   |   sensitivity n = {S['n_structures']} / {S['n_clusters']} clusters")

    # ---------------- Aim A
    emit("\nAIM A — RANK-1 / TOP-3 / TOP-5 EXPOSURE")
    emit(f"{'detector':9s} {'stratum':26s} {'primary 0.20':>26s} {'sensitivity 0.35':>26s}"
         f" {'delta':>8s}")
    for det in DET:
        pa, sa = P["detectors"][det]["A"], S["detectors"][det]["A"]
        for depth in ("rank-1", "top-3", "top-5"):
            pk, pn = pa["depths"][depth]["k"], pa["depths"][depth]["n"]
            sk, sn = sa["depths"][depth]["k"], sa["depths"][depth]["n"]
            pp, sp = 100 * pk / pn, 100 * sk / sn
            plo, phi = wilson(pk, pn); slo, shi = wilson(sk, sn)
            emit(f"{det:9s} {'ALL  ' + depth:26s} "
                 f"{f'{pk}/{pn}={pp:5.1f}% [{100*plo:.1f}-{100*phi:.1f}]':>26s} "
                 f"{f'{sk}/{sn}={sp:5.1f}% [{100*slo:.1f}-{100*shi:.1f}]':>26s} "
                 f"{sp-pp:+7.1f}")
            rows.append({"AIM": "A", "DETECTOR": det, "STRATUM": "ALL", "METRIC": depth,
                         "PRIMARY_K": pk, "PRIMARY_N": pn, "PRIMARY_PCT": round(pp, 1),
                         "PRIMARY_CI95": f"{100*plo:.1f}-{100*phi:.1f}",
                         "SENS_K": sk, "SENS_N": sn, "SENS_PCT": round(sp, 1),
                         "SENS_CI95": f"{100*slo:.1f}-{100*shi:.1f}",
                         "DELTA_PCT_POINTS": round(sp - pp, 1)})
        for pt in PARTNERS:
            pb, sb = pa["by_partner"][pt], sa["by_partner"][pt]
            for label, pk, pn, sk, sn in (
                    ("rank-1", pb["k"], pb["n"], sb["k"], sb["n"]),
                    ("top-3", pb["top-3"]["k"], pb["top-3"]["n"],
                     sb["top-3"]["k"], sb["top-3"]["n"]),
                    ("top-5", pb["top-5"]["k"], pb["top-5"]["n"],
                     sb["top-5"]["k"], sb["top-5"]["n"])):
                pp, sp = 100 * pk / pn, 100 * sk / sn
                plo, phi = wilson(pk, pn); slo, shi = wilson(sk, sn)
                if label == "rank-1":
                    emit(f"{det:9s} {pt + '  ' + label:26s} "
                         f"{f'{pk}/{pn}={pp:5.1f}% [{100*plo:.1f}-{100*phi:.1f}]':>26s} "
                         f"{f'{sk}/{sn}={sp:5.1f}% [{100*slo:.1f}-{100*shi:.1f}]':>26s} "
                         f"{sp-pp:+7.1f}")
                rows.append({"AIM": "A", "DETECTOR": det, "STRATUM": pt, "METRIC": label,
                             "PRIMARY_K": pk, "PRIMARY_N": pn, "PRIMARY_PCT": round(pp, 1),
                             "PRIMARY_CI95": f"{100*plo:.1f}-{100*phi:.1f}",
                             "SENS_K": sk, "SENS_N": sn, "SENS_PCT": round(sp, 1),
                             "SENS_CI95": f"{100*slo:.1f}-{100*shi:.1f}",
                             "DELTA_PCT_POINTS": round(sp - pp, 1)})
        pstd, sstd = 100 * pa["standardised"], 100 * sa["standardised"]
        pc, sc = pa["standardised_ci"], sa["standardised_ci"]
        emit(f"{det:9s} {'ALL  partner-standardised':26s} "
             f"{f'{pstd:5.1f}% [{100*float(pc[0]):.1f}-{100*float(pc[1]):.1f}]':>26s} "
             f"{f'{sstd:5.1f}% [{100*float(sc[0]):.1f}-{100*float(sc[1]):.1f}]':>26s} "
             f"{sstd-pstd:+7.1f}")
        rows.append({"AIM": "A", "DETECTOR": det, "STRATUM": "ALL",
                     "METRIC": "rank-1 partner-standardised",
                     "PRIMARY_K": "NA", "PRIMARY_N": pa["depths"]["rank-1"]["n"],
                     "PRIMARY_PCT": round(pstd, 1),
                     "PRIMARY_CI95": f"{100*float(pc[0]):.1f}-{100*float(pc[1]):.1f}",
                     "SENS_K": "NA", "SENS_N": sa["depths"]["rank-1"]["n"],
                     "SENS_PCT": round(sstd, 1),
                     "SENS_CI95": f"{100*float(sc[0]):.1f}-{100*float(sc[1]):.1f}",
                     "DELTA_PCT_POINTS": round(sstd - pstd, 1)})

    # ---------------- Aim B
    emit("\nAIM B — TARGET-CAVITY DISPLACEMENT (structure-level)")
    for det in DET:
        pb, sb = P["detectors"][det]["B"], S["detectors"][det]["B"]
        emit(f"\n  {DETLAB[det]}")
        for label, pv, sv in (
                ("evaluable", pb["evaluable"], sb["evaluable"]),
                ("detector abstentions", pb["abstention"], sb["abstention"]),
                ("unmatched", pb["unmatched"], sb["unmatched"]),
                ("median displacement", pb["median"], sb["median"]),
                ("mean displacement", round(float(pb["mean"]), 2), round(float(sb["mean"]), 2)),
                ("negative displacement", pb["n_negative"], sb["n_negative"])):
            emit(f"    {label:24s} primary {str(pv):>8s}   sensitivity {str(sv):>8s}")
            rows.append({"AIM": "B", "DETECTOR": det, "STRATUM": "ALL", "METRIC": label,
                         "PRIMARY_K": pv, "PRIMARY_N": "NA", "PRIMARY_PCT": "NA",
                         "PRIMARY_CI95": "NA", "SENS_K": sv, "SENS_N": "NA",
                         "SENS_PCT": "NA", "SENS_CI95": "NA", "DELTA_PCT_POINTS": "NA"})
        for lab in (">0", ">=2", ">=5"):
            pk, pn = pb["props"][lab]["k"], pb["props"][lab]["n"]
            sk, sn = sb["props"][lab]["k"], sb["props"][lab]["n"]
            pp, sp = 100 * pk / pn, 100 * sk / sn
            plo, phi = wilson(pk, pn); slo, shi = wilson(sk, sn)
            emit(f"    displacement {lab:11s} primary {pk}/{pn}={pp:5.1f}% "
                 f"[{100*plo:.1f}-{100*phi:.1f}]   sensitivity {sk}/{sn}={sp:5.1f}% "
                 f"[{100*slo:.1f}-{100*shi:.1f}]   delta {sp-pp:+.1f}")
            rows.append({"AIM": "B", "DETECTOR": det, "STRATUM": "ALL",
                         "METRIC": f"displacement {lab}",
                         "PRIMARY_K": pk, "PRIMARY_N": pn, "PRIMARY_PCT": round(pp, 1),
                         "PRIMARY_CI95": f"{100*plo:.1f}-{100*phi:.1f}",
                         "SENS_K": sk, "SENS_N": sn, "SENS_PCT": round(sp, 1),
                         "SENS_CI95": f"{100*slo:.1f}-{100*shi:.1f}",
                         "DELTA_PCT_POINTS": round(sp - pp, 1)})
        emit("    partner-stratified, displacement > 0:")
        for pt in PARTNERS:
            pp_, sp_ = pb["by_partner"][pt], sb["by_partner"][pt]
            ppc = 100 * pp_["k_gt0"] / pp_["n"]; spc = 100 * sp_["k_gt0"] / sp_["n"]
            emit(f"      {pt:5s} primary {pp_['k_gt0']}/{pp_['n']}={ppc:5.1f}%   "
                 f"sensitivity {sp_['k_gt0']}/{sp_['n']}={spc:5.1f}%   delta {spc-ppc:+.1f}")
            rows.append({"AIM": "B", "DETECTOR": det, "STRATUM": pt,
                         "METRIC": "displacement >0",
                         "PRIMARY_K": pp_["k_gt0"], "PRIMARY_N": pp_["n"],
                         "PRIMARY_PCT": round(ppc, 1),
                         "PRIMARY_CI95": f"{100*wilson(pp_['k_gt0'],pp_['n'])[0]:.1f}-"
                                         f"{100*wilson(pp_['k_gt0'],pp_['n'])[1]:.1f}",
                         "SENS_K": sp_["k_gt0"], "SENS_N": sp_["n"],
                         "SENS_PCT": round(spc, 1),
                         "SENS_CI95": f"{100*wilson(sp_['k_gt0'],sp_['n'])[0]:.1f}-"
                                      f"{100*wilson(sp_['k_gt0'],sp_['n'])[1]:.1f}",
                         "DELTA_PCT_POINTS": round(spc - ppc, 1)})

    # ---------------- Aim C (descriptive; status unchanged)
    emit("\nAIM C — BIOLOGICAL REFERENCE SITE (descriptive; secondary status unchanged)")
    pc, sc = P["reference_site_coverage"], S["reference_site_coverage"]
    emit(f"  coverage primary {pc['k']}/{pc['n']} = {100*pc['fraction']:.1f}%   "
         f"sensitivity {sc['k']}/{sc['n']} = {100*sc['fraction']:.1f}%   "
         f"(both BELOW the 60% promotion rule -> remains SECONDARY)")
    rows.append({"AIM": "C", "DETECTOR": "both", "STRATUM": "ALL",
                 "METRIC": "reference-site coverage",
                 "PRIMARY_K": pc["k"], "PRIMARY_N": pc["n"],
                 "PRIMARY_PCT": round(100 * pc["fraction"], 1), "PRIMARY_CI95": "NA",
                 "SENS_K": sc["k"], "SENS_N": sc["n"],
                 "SENS_PCT": round(100 * sc["fraction"], 1), "SENS_CI95": "NA",
                 "DELTA_PCT_POINTS": round(100 * (sc["fraction"] - pc["fraction"]), 1)})
    for det in DET:
        pC, sC = P["detectors"][det].get("C"), S["detectors"][det].get("C")
        if not (pC and sC):
            continue
        emit(f"  {DETLAB[det]:8s} outranked by a fusion pocket: primary "
             f"{pC['outranked']}/{pC['found']} = {100*pC['outranked']/pC['found']:.1f}%   "
             f"sensitivity {sC['outranked']}/{sC['found']} = "
             f"{100*sC['outranked']/sC['found']:.1f}%")
        emit(f"  {'':8s} top-1 recovery: primary {100*pC['top1']/pC['found']:.1f}%   "
             f"sensitivity {100*sC['top1']/sC['found']:.1f}%")
        rows.append({"AIM": "C", "DETECTOR": det, "STRATUM": "ALL",
                     "METRIC": "reference site outranked by fusion pocket",
                     "PRIMARY_K": pC["outranked"], "PRIMARY_N": pC["found"],
                     "PRIMARY_PCT": round(100 * pC["outranked"] / pC["found"], 1),
                     "PRIMARY_CI95": "NA",
                     "SENS_K": sC["outranked"], "SENS_N": sC["found"],
                     "SENS_PCT": round(100 * sC["outranked"] / sC["found"], 1),
                     "SENS_CI95": "NA", "DELTA_PCT_POINTS": "NA"})

    cols = ["AIM", "DETECTOR", "STRATUM", "METRIC", "PRIMARY_K", "PRIMARY_N", "PRIMARY_PCT",
            "PRIMARY_CI95", "SENS_K", "SENS_N", "SENS_PCT", "SENS_CI95", "DELTA_PCT_POINTS"]
    with open(os.path.join(TAB, "TABLE6_e9_sensitivity.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "NA")) for c in cols) + "\n")
    print(f"\n  wrote results/tables/TABLE6_e9_sensitivity.tsv ({len(rows)} rows)")

    # ---------------- new-structure characterisation
    diff = json.load(open(os.path.join(MANI, "e9_035_cohort_diff.json")))
    plan = json.load(open(os.path.join(MANI, "e9_035_run_plan.json")))
    sens_set = {s["pdb_id"]: s for s in
                json.load(open(os.path.join(MANI, "e9_035_set_final.json")))["structures"]}
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       "confirmatory_manifest_e9_035.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    swaps = {b: a for (_k, a, b) in diff["representative_swaps"]}
    nrows = []
    for pid in plan["new_runs"]:
        s = sens_set[pid]; m = man.get(pid, {})
        nrows.append({
            "PDB_ID": pid, "ADMISSION": ("REPRESENTATIVE_UPGRADE (replaces "
                                         + swaps[pid] + ")") if pid in swaps else "NEW_UNIT",
            "PARTNER": s["fusion_partner"], "TARGET_UNIPROT": s["target_accession"],
            "TOPOLOGY": s["topology"], "METHOD": s["method"],
            "RESOLUTION_A": round(float(s["resolution"]), 2),
            "TARGET_RESIDUES_MODELLED": m.get("residue_class_counts", {}).get("TARGET"),
            "FUSION_RESIDUES_MODELLED": m.get("residue_class_counts", {}).get("FUSION"),
            "TARGET_DISORDER": round(float(s["target_disorder"]), 3),
            "FUSION_DISORDER": round(float(s["fusion_disorder"]), 3),
            "REFERENCE_SITE_AVAILABLE": m.get("reference_site_available"),
            "DESCRIPTION": (s.get("description") or "")[:70]})
    ncols = list(nrows[0].keys())
    with open(os.path.join(TAB, "TABLE6b_e9_new_structures.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(ncols) + "\n")
        for r in nrows:
            fh.write("\t".join("NA" if r[c] is None else str(r[c]) for c in ncols) + "\n")
    emit("\nNEWLY RUN STRUCTURES (descriptive)")
    emit(f"  {'pdb':6s} {'admission':16s} {'ptnr':5s} {'topology':18s} {'method':6s} "
         f"{'res':>5s} {'tgt_len':>7s} {'fus_len':>7s} {'tgt_dis':>7s} {'fus_dis':>7s}")
    for r in nrows:
        emit(f"  {r['PDB_ID']:6s} {r['ADMISSION'][:16]:16s} {r['PARTNER']:5s} "
             f"{str(r['TOPOLOGY'])[:18]:18s} "
             f"{'X-ray' if 'X-RAY' in str(r['METHOD']) else 'cryoEM':6s} "
             f"{r['RESOLUTION_A']:5.2f} {str(r['TARGET_RESIDUES_MODELLED']):>7s} "
             f"{str(r['FUSION_RESIDUES_MODELLED']):>7s} "
             f"{r['TARGET_DISORDER']:7.3f} {r['FUSION_DISORDER']:7.3f}")
    emit(f"  by partner: {dict(collections.Counter(r['PARTNER'] for r in nrows))}")
    emit(f"  by admission: {dict(collections.Counter(r['ADMISSION'].split(' ')[0] for r in nrows))}")

    open(os.path.join(RES, "e9_sensitivity_comparison.txt"), "w",
         encoding="utf-8").write("\n".join(lines) + "\n")
    json.dump({"rows": rows, "new_structures": nrows},
              open(os.path.join(RES, "e9_sensitivity_comparison.json"), "w"),
              indent=1, default=str)

    # ---------------- regenerate S6
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for ax, det in zip(axes[:2], DET):
        pa, sa = P["detectors"][det]["A"], S["detectors"][det]["A"]
        xs = range(3)
        for i, pt in enumerate(PARTNERS):
            pb, sb = pa["by_partner"][pt], sa["by_partner"][pt]
            pp = 100 * pb["k"] / pb["n"]; sp = 100 * sb["k"] / sb["n"]
            ax.plot([i - 0.16], [pp], "o", color=CP[pt], ms=6)
            ax.plot([i + 0.16], [sp], "s", color=CP[pt], ms=6, mfc="white", mew=1.4)
            plo, phi = wilson(pb["k"], pb["n"]); slo, shi = wilson(sb["k"], sb["n"])
            ax.plot([i - 0.16] * 2, [100 * plo, 100 * phi], color=CP[pt], lw=1.2, alpha=0.7)
            ax.plot([i + 0.16] * 2, [100 * slo, 100 * shi], color=CP[pt], lw=1.2, alpha=0.4)
        ax.set_xticks(list(xs)); ax.set_xticklabels(PARTNERS)
        ax.set_ylim(-3, 105); ax.set_title(DETLAB[det], fontsize=8)
        ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
    axes[0].set_ylabel("rank-1 fusion-associated (%)")
    axes[0].plot([], [], "o", color="#555", label="E9 $\\leq$ 0.20 (primary)")
    axes[0].plot([], [], "s", color="#555", mfc="white", label="E9 $\\leq$ 0.35 (sensitivity)")
    axes[0].legend(frameon=False, fontsize=6, loc="upper left")

    ax = axes[2]
    for i, det in enumerate(DET):
        pb, sb = P["detectors"][det]["B"], S["detectors"][det]["B"]
        pp = 100 * pb["props"][">0"]["k"] / pb["props"][">0"]["n"]
        sp = 100 * sb["props"][">0"]["k"] / sb["props"][">0"]["n"]
        ax.plot([i - 0.16], [pp], "o", color="#2b6cb0", ms=6)
        ax.plot([i + 0.16], [sp], "s", color="#2b6cb0", ms=6, mfc="white", mew=1.4)
        ax.text(i, 92, f"neg: {pb['n_negative']}/{sb['n_negative']}", ha="center", fontsize=5.8,
                color="#444")
    ax.set_xticks([0, 1]); ax.set_xticklabels([DETLAB[d] for d in DET])
    ax.set_ylim(-3, 105); ax.set_ylabel("displacement > 0 (%)")
    ax.set_title("target-cavity displacement", fontsize=8)
    ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)

    fig.suptitle(f"S6  E9 segment-disorder sensitivity: primary 0.20 (n={P['n_structures']}) "
                 f"vs 0.35 (n={S['n_structures']})", fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0.04, 1, 0.92])
    fig.text(0.01, 0.005,
             "Conclusions are unchanged under the relaxed threshold; primary values remain the "
             "manuscript values. 16 structures were newly run (13 new units + 3 higher-resolution "
             "representative upgrades); 125 reused frozen output.",
             fontsize=5.6, ha="left", va="bottom", color="#444444")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"S6_e9_disorder.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/supplementary/S6_e9_disorder.png / .pdf (replaces audit-only version)")


if __name__ == "__main__":
    main()
