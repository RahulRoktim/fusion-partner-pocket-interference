#!/usr/bin/env python
"""
Phase 6F — supplementary figures S1..S10, all generated from frozen result files.

S6 note: protocol III.6 pre-specified an E9 disorder sensitivity at 0.35 for the CONFIRMATORY
analysis. That re-run was NOT executed (it would require a new cohort and new detector runs).
S6 therefore shows the E9 AUDIT evidence only, and is labelled as such. This gap is declared in
the claims ledger and the readiness assessment rather than papered over.
"""
import collections, csv, json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
CONF = os.path.join(RES, "confirmatory")
MANI = os.path.join(ROOT, "data_manifest")
FIG = os.path.join(ROOT, "figures", "supplementary")
os.makedirs(FIG, exist_ok=True)
PARTNERS = ["BRIL", "T4L", "MBP"]
DET = ["p2rank", "fpocket"]
DETLAB = {"p2rank": "P2Rank", "fpocket": "fpocket"}
CP = {"BRIL": "#4C72B0", "T4L": "#DD8452", "MBP": "#C44E52"}
plt.rcParams.update({"font.size": 7.5, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 200, "savefig.dpi": 350, "font.family": "DejaVu Sans"})


def tsv(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, name, note=None):
    if note:
        fig.text(0.01, 0.005, note, fontsize=5.6, ha="left", va="bottom", color="#444444")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/supplementary/{name}.png")


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0, c - h), min(1, c + h)


def s1_replication():
    rep = json.load(open(os.path.join(RES, "replication.json")))
    keys = [k for k in rep if k.startswith("A_") or k.startswith("B_")]
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    y = 0; yt, yl = [], []
    for k in sorted(keys):
        r = rep[k]
        dk, dn = r["dev"]; ck, cn = r["con"]
        if dn == 0 or cn == 0:
            continue
        dl, dh = wilson(dk, dn); cl, ch = wilson(ck, cn)
        ax.plot([100*dl, 100*dh], [y+0.16]*2, color="#BBBBBB", lw=1.4)
        ax.plot([100*dk/dn], [y+0.16], "o", color="#BBBBBB", ms=4)
        ax.plot([100*cl, 100*ch], [y-0.16]*2, color="#2b6cb0", lw=1.4)
        ax.plot([100*ck/cn], [y-0.16], "o", color="#2b6cb0", ms=4)
        v = r["verdict"]
        col = {"REPLICATED": "#2e7d32", "PARTIALLY REPLICATED": "#ef6c00",
               "NOT REPLICATED": "#c62828"}.get(v, "#555")
        ax.text(103, y, v.replace(" REPLICATED", ""), fontsize=5.8, color=col, va="center")
        yt.append(y); yl.append(k.replace("A_", "exposure ").replace("B_", "displacement "))
        y -= 1
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=6)
    ax.set_xlim(-2, 128); ax.set_xlabel("rate (%)")
    ax.plot([], [], "o", color="#BBBBBB", label="development (n=24, exploratory)")
    ax.plot([], [], "o", color="#2b6cb0", label="confirmatory (n=128)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    ax.set_title("S1  Development vs confirmatory replication (cohorts never pooled)", fontsize=8)
    save(fig, "S1_replication")


def s2_s3_sensitivity():
    sens = json.load(open(os.path.join(CONF, "primary",
                                       "confirmatory_summary.json")))["sensitivity"]
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, det in zip(axes, DET):
        d = sens["dominance"][det]
        xs = ["0.50", "0.70", "0.90"]
        ys = [100 * d[t]["k"] / d[t]["n"] for t in xs]
        cis = [wilson(d[t]["k"], d[t]["n"]) for t in xs]
        ax.plot(range(3), ys, "-o", color="#2b6cb0", ms=5)
        for i, (lo, hi) in enumerate(cis):
            ax.plot([i, i], [100*lo, 100*hi], color="#2b6cb0", lw=1.1, alpha=0.6)
        ax.set_xticks(range(3)); ax.set_xticklabels([f"{x}\n{'(primary)' if x=='0.70' else ''}"
                                                     for x in xs])
        ax.set_ylim(0, 75); ax.set_title(DETLAB[det], fontsize=8)
        ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
    axes[0].set_ylabel("rank-1 fusion-associated (%)")
    fig.suptitle("S2  Pocket-dominance threshold sensitivity", fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92]); save(fig, "S2_dominance_sensitivity")

    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, det in zip(axes, DET):
        d = sens["correspondence"][det]
        labs = list(d)
        ys = [100 * d[l]["k_gt0"] / d[l]["n"] for l in labs]
        cis = [wilson(d[l]["k_gt0"], d[l]["n"]) for l in labs]
        ax.plot(range(len(labs)), ys, "-o", color="#6a51a3", ms=5)
        for i, (lo, hi) in enumerate(cis):
            ax.plot([i, i], [100*lo, 100*hi], color="#6a51a3", lw=1.1, alpha=0.6)
        ax.set_xticks(range(len(labs)))
        ax.set_xticklabels([l + ("\n(primary)" if "0.40" in l else "") for l in labs], fontsize=6)
        ax.set_ylim(0, 75); ax.set_title(DETLAB[det], fontsize=8)
        ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
    axes[0].set_ylabel("displacement > 0 (%)")
    fig.suptitle("S3  Target-cavity correspondence threshold sensitivity",
                 fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92]); save(fig, "S3_correspondence_sensitivity")


def s4_s5_variants():
    out = {}
    for v in ("primary", "ions", "assembly1"):
        p = os.path.join(CONF, v, "confirmatory_summary.json")
        if os.path.exists(p):
            out[v] = json.load(open(p))
    a1 = os.path.join(CONF, "assembly1", "displacement_chainaware.json")
    a1d = json.load(open(a1))["summary"] if os.path.exists(a1) else {}
    for name, title, variants in (("S4_assembly1_sensitivity",
                                   "S4  Biological assembly 1 vs single chain",
                                   ["primary", "assembly1"]),
                                  ("S5_hetatm_ions_sensitivity",
                                   "S5  HETATM handling: stripped vs ions retained",
                                   ["primary", "ions"])):
        fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9))
        for ax, det in zip(axes, DET):
            xs, ys, cis, labs = [], [], [], []
            for i, v in enumerate(variants):
                if v not in out:
                    continue
                A = out[v]["detectors"][det]["A"]["depths"]["rank-1"]
                xs.append(i); ys.append(100 * A["k"] / A["n"])
                cis.append(wilson(A["k"], A["n"]))
                labs.append(f"{v}\nn={A['n']}")
            ax.bar(xs, ys, color=["#2b6cb0", "#93a8c4"], width=0.5)
            for i, (lo, hi) in enumerate(cis):
                ax.plot([xs[i], xs[i]], [100*lo, 100*hi], color="#222", lw=1.0)
            ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=6.5)
            ax.set_ylim(0, 75); ax.set_title(DETLAB[det], fontsize=8)
            ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
        axes[0].set_ylabel("rank-1 fusion-associated (%)")
        note = None
        if "assembly1" in variants and a1d:
            note = ("assembly-1 displacement (chain-aware): "
                    + " | ".join(f"{d} {a1d[d]['k_gt0']}/{a1d[d]['n']} >0, median "
                                 f"{a1d[d]['median']}" for d in a1d if d in ("p2rank", "fpocket"))
                    + ".  2 MBP structures (5W0R, 5EDU) are not evaluable under assembly 1.")
        fig.suptitle(title, fontsize=8.5, x=0.02, ha="left")
        fig.tight_layout(rect=[0, 0, 1, 0.92]); save(fig, name, note)


def s6_e9():
    rows = tsv(os.path.join(RES, "e9_audit", "e9_distributions.tsv"))
    chim = [r for r in rows if r["target_disorder"] not in ("", "None")
            and r["fusion_disorder"] not in ("", "None")
            and r["target_accession"] not in ("", "None")]
    td = sorted(float(r["target_disorder"]) for r in chim)
    fd = sorted(float(r["fusion_disorder"]) for r in chim)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, vals, lab in ((axes[0], td, "target segment"), (axes[1], fd, "fusion segment")):
        ax.hist(vals, bins=40, color="#6a8caf", edgecolor="white", lw=0.3)
        ax.axvline(0.20, color="#c62828", lw=1.2, label="0.20 (primary)")
        ax.axvline(0.35, color="#ef6c00", lw=1.2, ls="--", label="0.35 (pre-specified\nsensitivity, NOT RUN)")
        ax.set_xlabel(f"{lab} disorder fraction"); ax.set_ylabel("chimeric entities")
        ax.legend(frameon=False, fontsize=5.8)
    fig.suptitle("S6  E9 segment-disorder audit (genuinely chimeric entities only)",
                 fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    save(fig, "S6_e9_disorder",
         "E9's marginal effect on the eligible pool is +37 entities (455 -> 492, +8.1%). "
         "The pre-specified confirmatory re-run at 0.35 was NOT executed; declared as a gap.")


def s7_abstention():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       "confirmatory_manifest_primary.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(CONF, "primary", "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    for det, mk, col in (("p2rank", "o", "#2b6cb0"), ("fpocket", "s", "#999999")):
        xs, ys = [], []
        for s in sorted(man):
            n = man[s]["residue_class_counts"].get("TARGET", 0)
            xs.append(n); ys.append(len(idx.get((det, s, "REMOVED"), [])))
        ax.scatter(xs, ys, s=11, marker=mk, color=col, alpha=0.7, label=DETLAB[det],
                   edgecolors="none")
    ax.set_xlabel("modelled target residues (fusion-removed structure)")
    ax.set_ylabel("pockets predicted")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(color="#f2f2f2", lw=0.6); ax.set_axisbelow(True)
    ax.set_title("S7  P2Rank abstains on small targets; fpocket never abstains", fontsize=8)
    save(fig, "S7_abstention_by_target_size")


def s8_margins():
    P = json.load(open(os.path.join(CONF, "primary", "pockets_classified.json")))["pockets"]
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       "confirmatory_manifest_primary.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9), sharey=True)
    for ax, det in zip(axes, DET):
        data, labs = [], []
        for pt in PARTNERS:
            v = []
            for s in sorted(man):
                if man[s]["fusion_partner"] != pt:
                    continue
                pk = idx.get((det, s, "ORIGINAL"), [])
                bf = max((p["score"] for p in pk if p["fusion_associated"]), default=None)
                bt = max((p["score"] for p in pk if p["pocket_class"] == "TARGET_DOMINATED"),
                         default=None)
                if bf is not None and bt:
                    v.append(bf / bt)
            data.append(v); labs.append(f"{pt}\nn={len(v)}")
        bp = ax.boxplot(data, tick_labels=labs, showfliers=False, widths=0.55, patch_artist=True)
        for patch, pt in zip(bp["boxes"], PARTNERS):
            patch.set_facecolor(CP[pt]); patch.set_alpha(0.55); patch.set_edgecolor("#333")
        for m in bp["medians"]:
            m.set_color("#111"); m.set_linewidth(1.2)
        ax.axhline(1.0, color="#c62828", lw=1.0, ls="--")
        ax.set_yscale("log"); ax.set_title(DETLAB[det], fontsize=8)
        ax.grid(axis="y", color="#f2f2f2", lw=0.6); ax.set_axisbelow(True)
    axes[0].set_ylabel("best fusion-associated score /\nbest target-dominated score")
    fig.suptitle("S8  Score margin by partner (dashed line = parity)",
                 fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.91]); save(fig, "S8_score_margins")


def s9_agreement():
    summ = json.load(open(os.path.join(CONF, "primary", "confirmatory_summary.json")))
    E = summ["E"]
    P = json.load(open(os.path.join(CONF, "primary", "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])
    structs = sorted({r["pdb_id"] for r in P})
    both = [s for s in structs if idx.get(("p2rank", s, "ORIGINAL"))
            and idx.get(("fpocket", s, "ORIGINAL"))]
    a = [idx[("p2rank", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    b = [idx[("fpocket", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    m = [[sum(1 for x, y in zip(a, b) if x and y), sum(1 for x, y in zip(a, b) if x and not y)],
         [sum(1 for x, y in zip(a, b) if not x and y),
          sum(1 for x, y in zip(a, b) if not x and not y)]]
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    im = ax.imshow(m, cmap="Blues", vmin=0)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, m[i][j], ha="center", va="center", fontsize=11,
                    color="white" if m[i][j] > max(max(m)) * 0.55 else "#111")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["fpocket\nfusion-assoc", "fpocket\ntarget"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["P2Rank\nfusion-assoc", "P2Rank\ntarget"])
    ax.set_title(f"S9  Rank-1 call agreement\nn={len(both)}, "
                 f"agreement {100*float(E['agreement']):.1f}%, "
                 f"$\\kappa$ = {float(E['kappa']):.3f}", fontsize=8)
    save(fig, "S9_detector_agreement")


def s10_refsite():
    t4 = tsv(os.path.join(TAB, "TABLE4_reference_site.tsv"))
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    x = range(3)
    for i, det in enumerate(DET):
        r = next(z for z in t4 if z["DETECTOR"] == det)
        ys = [float(r["RECOVERY_TOP1_PCT"]), float(r["RECOVERY_TOP3_PCT"]),
              float(r["RECOVERY_TOP5_PCT"])]
        ax.plot(x, ys, "-o", color=["#2b6cb0", "#999999"][i], ms=5,
                label=f"{DETLAB[det]}  (outranked by a fusion pocket in "
                      f"{r['OUTRANKED_BY_FUSION_K']}/{r['OUTRANKED_BY_FUSION_N']})")
    ax.set_xticks(list(x)); ax.set_xticklabels(["top 1", "top 3", "top 5"])
    ax.set_ylim(0, 105); ax.set_ylabel("biological reference site recovered (%)")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    ax.grid(axis="y", color="#f2f2f2", lw=0.6); ax.set_axisbelow(True)
    cov = next(z for z in t4 if z["DETECTOR"] == "p2rank")
    ax.set_title(f"S10  Biological reference-site subset "
                 f"(coverage {cov['COVERAGE_K']}/{cov['COVERAGE_N']} = "
                 f"{cov['COVERAGE_PCT']}%, below the 60% promotion rule)", fontsize=7.6)
    save(fig, "S10_reference_site")


if __name__ == "__main__":
    print("generating supplementary figures ...")
    s1_replication(); s2_s3_sensitivity(); s4_s5_variants(); s6_e9()
    s7_abstention(); s8_margins(); s9_agreement(); s10_refsite()
