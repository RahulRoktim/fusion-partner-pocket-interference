#!/usr/bin/env python
"""
Phase 6F — main figures 1, 2, 3, 4, 6. Every number is read from a canonical frozen result file;
none is typed into this script.

Outputs: figures/FIG1..FIG6 (.png and .pdf)
"""
import collections, csv, json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
CONF = os.path.join(RES, "confirmatory", "primary")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

PARTNERS = ["BRIL", "T4L", "MBP"]
DET = ["p2rank", "fpocket"]
DETLAB = {"p2rank": "P2Rank 2.5.1", "fpocket": "fpocket 4.2.3"}
CP = {"BRIL": "#4C72B0", "T4L": "#DD8452", "MBP": "#C44E52"}
CD = {"p2rank": "#2b2b2b", "fpocket": "#7f7f7f"}

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8, "figure.dpi": 200, "savefig.dpi": 400,
                     "font.family": "DejaVu Sans"})


def tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


def ci(s):
    lo, hi = s.split("-")[0], s.split("-")[1].split(" ")[0]
    return float(lo), float(hi)


# ---------------------------------------------------------------- FIGURE 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 46); ax.axis("off")
    boxes = [
        (2, 30, 17, 11, "Deposited PDB entry\nchimeric construct", "#EAEFF7"),
        (23, 30, 17, 11, "SIFTS residue-level\ntarget / fusion / linker\nannotation", "#EAEFF7"),
        (44, 36, 17, 8.5, "ORIGINAL\ntarget + fusion + linker", "#FBEEE6"),
        (44, 24, 17, 8.5, "FUSION-REMOVED\ntarget only\n(target coords unchanged)", "#E8F2EA"),
        (65, 30, 14, 11, "P2Rank 2.5.1\nfpocket 4.2.3\nstock parameters", "#F2F2F2"),
        (83, 30, 15, 11, "pocket classes\nTARGET / FUSION\nINTERFACE / LINKER", "#F2F2F2"),
        (44, 8, 36, 10, "SAME-TARGET-CAVITY CORRESPONDENCE\nJaccard on target residues "
                        "$\\geq$ 0.40 and centroid $\\leq$ 8 $\\AA$\none-to-one; rank and score "
                        "never used for matching", "#EDE7F3"),
        (83, 8, 15, 10, "RANK DISPLACEMENT\noriginal rank\n$-$ removed rank", "#EDE7F3"),
        (2, 8, 36, 10, "DEVELOPMENT n = 24 (exploratory)\n$\\rightarrow$ frozen protocol "
                       "$\\rightarrow$\nCONFIRMATORY n = 128, target-disjoint", "#FDF6DC"),
    ]
    for x, y, w, h, t, c in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.45",
                                    fc=c, ec="#555555", lw=0.7))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=6.6)

    def arrow(x1, y1, x2, y2, style="-|>"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                     mutation_scale=9, lw=0.8, color="#444444",
                                     connectionstyle="arc3,rad=0"))
    arrow(19, 35.5, 23, 35.5)
    arrow(40, 35.5, 44, 40)
    arrow(40, 35.5, 44, 28)
    arrow(61, 40, 65, 36.5)
    arrow(61, 28, 65, 34.5)
    arrow(79, 35.5, 83, 35.5)
    arrow(90, 30, 90, 18)
    arrow(62, 24, 62, 18)
    arrow(80, 13, 83, 13)
    ax.text(50, 44.5, "Figure 1  Study design", ha="center", fontsize=9, fontweight="bold")
    save(fig, "FIG1_study_design")


# ---------------------------------------------------------------- FIGURE 2
def fig2():
    t2 = tsv(os.path.join(TAB, "TABLE2_confirmatory_exposure.tsv"))
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ypos, labels = [], []
    y = 0
    for pt in PARTNERS:
        for det in DET:
            r = next(x for x in t2 if x["DETECTOR"] == det and x["STRATUM"] == pt)
            k, n = int(r["RANK1_K"]), int(r["RANK1_N"])
            p = 100 * k / n
            lo, hi = ci(r["RANK1_CI95"])
            ax.plot([lo, hi], [y, y], color=CP[pt], lw=1.6,
                    alpha=0.55 if det == "fpocket" else 1.0, solid_capstyle="butt")
            ax.plot([p], [y], "o" if det == "p2rank" else "s", color=CP[pt], ms=5,
                    mfc=CP[pt] if det == "p2rank" else "white", mew=1.2)
            ax.text(hi + 1.5, y, f"{k}/{n}", va="center", fontsize=6.3, color="#555555")
            ypos.append(y); labels.append(f"{pt}  {DETLAB[det].split()[0]}")
            y -= 1
        y -= 0.6
    # standardised overall
    y -= 0.3
    for det in DET:
        r = next(x for x in t2 if x["DETECTOR"] == det
                 and x["STRATUM"].startswith("ALL (partner-standardised"))
        p = float(r["RANK1_PCT"])
        lo, hi = ci(r["RANK1_CI95"])
        ax.plot([lo, hi], [y, y], color="#333333", lw=1.6,
                alpha=0.55 if det == "fpocket" else 1.0)
        ax.plot([p], [y], "o" if det == "p2rank" else "s", color="#333333", ms=5,
                mfc="#333333" if det == "p2rank" else "white", mew=1.2)
        ypos.append(y); labels.append(f"standardised  {DETLAB[det].split()[0]}")
        y -= 1
    ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("structures with a fusion-associated pocket at rank 1 (%)")
    ax.set_xlim(-2, 108)
    ax.axvline(0, color="#cccccc", lw=0.6, zorder=0)
    ax.grid(axis="x", color="#eeeeee", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("Figure 2  Confirmatory rank-1 exposure (n = 128, target-disjoint)\n"
                 "point estimate with 95% Wilson interval; circles P2Rank, squares fpocket",
                 fontsize=8, loc="left")
    save(fig, "FIG2_confirmatory_exposure")


# ---------------------------------------------------------------- FIGURE 3
def fig3():
    t2 = tsv(os.path.join(TAB, "TABLE2_confirmatory_exposure.tsv"))
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), sharey=True)
    depths = [("RANK1", "rank 1"), ("TOP3", "top 3"), ("TOP5", "top 5")]
    for ax, det in zip(axes, DET):
        for pt in PARTNERS:
            r = next(x for x in t2 if x["DETECTOR"] == det and x["STRATUM"] == pt)
            ys = [float(r[f"{d}_PCT"]) for d, _ in depths]
            ax.plot(range(3), ys, "-o", color=CP[pt], ms=4.5, lw=1.4, label=pt)
            for i, (d, _) in enumerate(depths):
                lo, hi = ci(r[f"{d}_CI95"])
                ax.plot([i, i], [lo, hi], color=CP[pt], lw=1.0, alpha=0.55)
        ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in depths])
        ax.set_title(DETLAB[det], fontsize=8)
        ax.set_ylim(-3, 105)
        ax.grid(axis="y", color="#eeeeee", lw=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("structures with a fusion-associated\npocket within rank depth (%)")
    axes[0].legend(frameon=False, fontsize=7, loc="center left")
    fig.suptitle("Figure 3  Rank-depth behaviour: rank-1 interference is MBP-dominated, "
                 "lower-ranked exposure is broad", fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save(fig, "FIG3_rank_depth")


# ---------------------------------------------------------------- FIGURE 4
def fig4():
    rows = tsv(os.path.join(CONF, "displacement.tsv"))
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.2, 1.0], hspace=0.55, wspace=0.28)
    for j, det in enumerate(DET):
        sub = [r for r in rows if r["detector"] == det]
        ev = [r for r in sub if r["state"] == "EVALUABLE"]
        ax = fig.add_subplot(gs[0, j])
        for r in ev:
            o, m = int(r["original_rank"]), int(r["removed_rank"])
            ax.plot([0, 1], [o, m], "-", color=CP[r["partner"]], lw=0.7,
                    alpha=0.75 if o != m else 0.25)
            ax.plot([0], [o], "o", color=CP[r["partner"]], ms=2.4, alpha=0.8)
        ax.set_xlim(-0.25, 1.25); ax.set_xticks([0, 1])
        ax.set_xticklabels(["original\n(fusion present)", "fusion-removed"], fontsize=7)
        ax.set_ylim(max(int(r['original_rank']) for r in ev) + 1.0, 0.4)
        ax.set_ylabel("rank of the same target cavity" if j == 0 else "")
        ax.set_title(f"{DETLAB[det]}   n evaluable = {len(ev)}", fontsize=8)
        ax.grid(axis="y", color="#f0f0f0", lw=0.6); ax.set_axisbelow(True)
        if j == 0:
            for pt in PARTNERS:
                ax.plot([], [], "-", color=CP[pt], lw=1.4, label=pt)
            ax.legend(frameon=False, fontsize=7, loc="lower left")

        axb = fig.add_subplot(gs[1, j])
        states = collections.Counter(r["state"] for r in sub)
        order = ["EVALUABLE", "DETECTOR_ABSTENTION", "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL"]
        lab = ["evaluable", "detector\nabstention", "target cavity not\nrecovered in original"]
        vals = [states.get(o, 0) for o in order]
        b = axb.barh(range(3), vals, color=["#4C72B0", "#BBBBBB", "#888888"], height=0.6)
        axb.set_yticks(range(3)); axb.set_yticklabels(lab, fontsize=6.4)
        axb.invert_yaxis()
        for i, v in enumerate(vals):
            axb.text(v + 2, i, str(v), va="center", fontsize=6.5)
        axb.set_xlim(0, 140)
        axb.set_xlabel("structures", fontsize=7)
        axb.grid(axis="x", color="#f0f0f0", lw=0.6); axb.set_axisbelow(True)
    fig.suptitle("Figure 4  Target-cavity rank displacement, one structure = one observation.\n"
                 "States that carry no numerical rank are shown explicitly, never imputed.",
                 fontsize=8.5, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save(fig, "FIG4_displacement")


# ---------------------------------------------------------------- FIGURE 6
def fig6():
    t5 = tsv(os.path.join(TAB, "TABLE5_dataset_provenance.tsv"))
    pl = tsv(os.path.join(RES, "PLINDER_DENOMINATORS.tsv"))
    bench = json.load(open(os.path.join(ROOT, "dataset_audit", "benchmark_impact.json")))
    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.62, wspace=0.3)

    # Panel A
    axA = fig.add_subplot(gs[0, :])
    names = ["CHEN11", "JOINED", "COACH420", "HOLO4K", "FPTRAIN"]
    xs, hs, ls = [], [], []
    for i, nm in enumerate(names):
        r = next(x for x in t5 if x["DATASET"] == nm)
        n = int(r["N_RESOLVED"]); k = int(r["FUSION_CONSTRUCTS"])
        xs.append(i); hs.append(100 * k / n)
        ls.append(f"{nm}\n{r['ROLE'].split('(')[0].strip()}\nn={n} entries")
    axA.bar(xs, hs, color="#4C72B0", width=0.55)
    for i, (h, nm) in enumerate(zip(hs, names)):
        r = next(x for x in t5 if x["DATASET"] == nm)
        axA.text(i, h + 0.012, f"{r['FUSION_CONSTRUCTS']}/{r['N_RESOLVED']}",
                 ha="center", fontsize=6.6)
    axA.set_xticks(xs); axA.set_xticklabels(ls, fontsize=6.2)
    axA.set_ylabel("crystallization fusion\nconstructs (% of entries)", fontsize=7)
    axA.set_ylim(0, 0.30)
    axA.set_title("A   Datasets that trained, developed and tested P2Rank and fpocket",
                  fontsize=8, loc="left")
    axA.grid(axis="y", color="#f0f0f0", lw=0.6); axA.set_axisbelow(True)

    # Panel B
    axB = fig.add_subplot(gs[1, 0])
    splits = ["train", "val", "test"]
    fs = [int(next(x for x in pl if x["SPLIT"] == s)["FUSION_SYSTEMS"]) for s in splits]
    fl = [int(next(x for x in pl if x["SPLIT"] == s)["FUSION_SITE_LABELLED_SYSTEMS"])
          for s in splits]
    fse = [int(next(x for x in pl if x["SPLIT"] == s)["FUSION_PDB_ENTRIES"]) for s in splits]
    fle = [int(next(x for x in pl if x["SPLIT"] == s)["FUSION_SITE_LABELLED_PDB_ENTRIES"])
           for s in splits]
    x = range(3)
    axB.bar([i - 0.19 for i in x], fs, width=0.36, color="#BBBBBB",
            label="construct receptor (systems)")
    axB.bar([i + 0.19 for i in x], fl, width=0.36, color="#C44E52",
            label="site labelled ON the fusion (systems)")
    for i in x:
        axB.text(i - 0.19, fs[i] + 25, f"{fs[i]}\n({fse[i]} entries)", ha="center", fontsize=5.8)
        axB.text(i + 0.19, fl[i] + 25, f"{fl[i]}\n({fle[i]} entries)", ha="center", fontsize=5.8)
    axB.set_xticks(list(x)); axB.set_xticklabels(
        [f"{s}\n{int(next(y for y in pl if y['SPLIT']==s)['SYSTEMS']):,} systems"
         for s in splits], fontsize=6.4)
    axB.set_ylabel("PLINDER systems", fontsize=7)
    axB.set_ylim(0, 1750)
    axB.legend(frameon=False, fontsize=5.8, loc="upper right")
    axB.set_title("B   PLINDER 2024-06/v2, by split\n(system and unique-entry counts both shown)",
                  fontsize=8, loc="left")
    axB.grid(axis="y", color="#f0f0f0", lw=0.6); axB.set_axisbelow(True)

    # Panel C
    axC = fig.add_subplot(gs[1, 1])
    c = bench["cases"][0]
    comp = c["gt_site_composition"]
    vals = [comp.get("FUSION", 0), comp.get("TARGET", 0)]
    axC.barh([1, 0], vals, color=["#C44E52", "#4C72B0"], height=0.5)
    axC.set_yticks([1, 0])
    axC.set_yticklabels(["contacts to GST\ncarrier", "contacts to fibrinogen\ntarget"],
                        fontsize=6.5)
    for i, v in enumerate(vals):
        axC.text(v + 1.5, 1 - i, str(v), va="center", fontsize=7)
    axC.set_xlim(0, 88)
    axC.set_xlabel("ground-truth ligand (GSH) atom contacts", fontsize=7)
    p2 = c["detectors"]["p2rank"]
    axC.set_title(f"C   1DUG, the one construct in a P2Rank-related set\n"
                  f"(JOINED / DT198, development). P2Rank recovers the\n"
                  f"labelled site at rank {p2['gt_rank']} with a 100% carrier pocket.",
                  fontsize=7.4, loc="left")
    axC.grid(axis="x", color="#f0f0f0", lw=0.6); axC.set_axisbelow(True)

    fig.suptitle("Figure 6  Dataset provenance audit (secondary analysis)",
                 fontsize=9, x=0.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "FIG6_dataset_provenance")


if __name__ == "__main__":
    print("generating main figures ...")
    fig1(); fig2(); fig3(); fig4(); fig6()
