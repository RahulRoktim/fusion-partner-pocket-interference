#!/usr/bin/env python
"""
Phase 6F/6H — Figure 5, structural case studies, plus the validation table the figure must pass.

Selection is deterministic and stated, not curated by eye:
  for each mechanistic class below, take the CONFIRMATORY structure with the best (lowest)
  resolution among structures whose rank-1 pocket is fusion-associated in BOTH detectors where
  possible, else in P2Rank, else in fpocket.
    C1  MBP,  FUSION_INTRINSIC_CAVITY
    C2  T4L,  FUSION_INTRINSIC_CAVITY
    C3  BRIL, FUSION_INTRINSIC_CAVITY
    C4  any,  TARGET_FUSION_INTERFACE
These are illustrative, NOT statistically representative, and the figure caption says so.

Geometry is read from the deposited coordinates and the frozen prepared inputs. Nothing is
reconstructed by hand. A PyMOL script reproducing each panel at publication quality is emitted
alongside, for later rendering in a molecular viewer.

Outputs:
  figures/FIG5_structural_cases.png/.pdf
  figures/FIG5_case_validation.tsv
  figures/FIG5_pymol.pml
"""
import collections, csv, json, math, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
CONF = os.path.join(RES, "confirmatory", "primary")
MANI = os.path.join(ROOT, "data_manifest")
FIG = os.path.join(ROOT, "figures")
CACHE = os.path.join(ROOT, ".structure_cache")
os.makedirs(FIG, exist_ok=True)
CP = {"TARGET": "#4C72B0", "FUSION": "#C44E52", "LINKER": "#8C8C8C"}
LIG_DIST = 5.0


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       "confirmatory_manifest_primary.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    mech = json.load(open(os.path.join(CONF, "mechanism.json")))["rows"]
    P = json.load(open(os.path.join(CONF, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])
    corr = [l.rstrip("\n").split("\t") for l in
            open(os.path.join(CONF, "displacement.tsv"), encoding="utf-8")]
    dh = corr[0]
    disp = [dict(zip(dh, r)) for r in corr[1:]]

    # ---------------- deterministic selection
    def pick(partner, mechanism):
        cands = collections.defaultdict(set)
        for m in mech:
            if m["mechanism"] != mechanism:
                continue
            if partner and m["partner"] != partner:
                continue
            cands[m["pdb_id"]].add(m["detector"])
        if not cands:
            return None
        both = [p for p, d in cands.items() if len(d) == 2]
        pool = both or list(cands)
        return sorted(pool, key=lambda p: (man[p]["resolution"], p))[0]

    cases = []
    for label, partner, mechanism in (
            ("C1", "MBP", "FUSION_INTRINSIC_CAVITY"),
            ("C2", "T4L", "FUSION_INTRINSIC_CAVITY"),
            ("C3", "BRIL", "FUSION_INTRINSIC_CAVITY"),
            ("C4", None, "TARGET_FUSION_INTERFACE")):
        p = pick(partner, mechanism)
        if p:
            cases.append({"label": label, "pdb_id": p, "class": mechanism,
                          "partner": man[p]["fusion_partner"]})
    print("selected cases:", [(c["label"], c["pdb_id"], c["partner"], c["class"]) for c in cases])

    # ---------------- gather geometry + validation facts
    for c in cases:
        s = c["pdb_id"]
        rec = man[s]
        c["chain"] = rec["auth_chain"]
        c["resolution"] = round(float(rec["resolution"]), 2)
        c["method"] = rec["method"]
        c["target_accession"] = rec["target_accession"]
        c["fusion_accession"] = rec["fusion_accession"]
        c["topology"] = rec["topology"]
        c["description"] = rec["description"]
        c["target_ranges_entity"] = rec["target_ranges_entity"]
        c["fusion_ranges_entity"] = rec["fusion_ranges_entity"]
        c["n_target"] = rec["residue_class_counts"].get("TARGET", 0)
        c["n_fusion"] = rec["residue_class_counts"].get("FUSION", 0)
        cmap = rec["auth_class_map"]

        # CA coordinates by class, from the frozen prepared ORIGINAL input
        st = gemmi.read_structure(os.path.join(ROOT, rec["original_pdb"]))
        pts = {"TARGET": [], "FUSION": [], "LINKER": []}
        for ch in st[0]:
            for res in ch:
                ca = res.find_atom("CA", "*")
                if ca is None:
                    continue
                k = f"{res.seqid.num}{res.seqid.icode.strip()}"
                cl = cmap.get(k)
                if cl in pts:
                    pts[cl].append((ca.pos.x, ca.pos.y, ca.pos.z))
        c["ca"] = pts

        # the offending rank-1 pocket, per detector
        c["detectors"] = {}
        for det in ("p2rank", "fpocket"):
            pk = idx.get((det, s, "ORIGINAL"), [])
            if not pk:
                continue
            top = pk[0]
            mm = next((m for m in mech if m["pdb_id"] == s and m["detector"] == det), None)
            c["detectors"][det] = {
                "rank1_class": top["pocket_class"], "rank1_score": top["score"],
                "f_fusion": top["f_fusion"], "f_target": top["f_target"],
                "n_residues": top["n_residues"],
                "centroid": (top["centroid_x"], top["centroid_y"], top["centroid_z"]),
                "residues": top["residues"].split(";") if top["residues"] else [],
                "mechanism": mm["mechanism"] if mm else None,
                "isolated_rank": mm["best_isolated_rank"] if mm else None,
                "probe_state": mm["probe_state"] if mm else None}
            d = next((x for x in disp if x["pdb_id"] == s and x["detector"] == det), None)
            if d and d["state"] == "EVALUABLE":
                c["detectors"][det]["target_cavity_original_rank"] = int(d["original_rank"])
                c["detectors"][det]["target_cavity_removed_rank"] = int(d["removed_rank"])
                c["detectors"][det]["displacement"] = int(d["displacement"])
            else:
                c["detectors"][det]["target_cavity_state"] = d["state"] if d else "NA"

        # deposited groups inside the offending pocket, from the DEPOSITED file
        dep = gemmi.read_structure(os.path.join(CACHE, f"{s.lower()}.cif"))
        dep.setup_entities()
        model = dep[0]
        want = set(c["detectors"].get("p2rank", {}).get("residues")
                   or c["detectors"].get("fpocket", {}).get("residues") or [])
        anchor = []
        for ch in model:
            if ch.name != c["chain"]:
                continue
            for res in ch:
                if f"{res.seqid.num}{res.seqid.icode.strip()}" in want:
                    anchor.append(res)
        hits = collections.Counter()
        ligpts = []
        for ch in model:
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                    continue
                if res.name in ("HOH", "DOD"):
                    continue
                close = False
                for a in res:
                    for t in anchor:
                        if any(a.pos.dist(b.pos) <= LIG_DIST for b in t):
                            close = True
                            break
                    if close:
                        break
                if close:
                    hits[res.name] += 1
                    ligpts += [(a.pos.x, a.pos.y, a.pos.z) for a in res]
        c["deposited_groups_in_pocket"] = dict(hits)
        c["ligand_points"] = ligpts

    # ---------------- validation table (6H)
    vcols = ["CASE", "PDB_ID", "CHAIN", "PARTNER", "MECHANISM_CLASS", "RESOLUTION_A", "METHOD",
             "TARGET_UNIPROT", "FUSION_UNIPROT", "TOPOLOGY",
             "TARGET_RANGES_ENTITY", "FUSION_RANGES_ENTITY",
             "N_TARGET_RESIDUES_MODELLED", "N_FUSION_RESIDUES_MODELLED",
             "DETECTOR", "RANK1_CLASS", "RANK1_SCORE", "RANK1_F_FUSION", "RANK1_N_RESIDUES",
             "PROBE_STATE", "ISOLATED_PARTNER_RANK",
             "TARGET_CAVITY_ORIGINAL_RANK", "TARGET_CAVITY_REMOVED_RANK", "DISPLACEMENT",
             "DEPOSITED_GROUPS_IN_POCKET", "ENTRY_TITLE"]
    vrows = []
    for c in cases:
        title = json.load(open(os.path.join(MANI, "confirmatory_manifest_primary.json")))
        for det, d in c["detectors"].items():
            vrows.append({
                "CASE": c["label"], "PDB_ID": c["pdb_id"], "CHAIN": c["chain"],
                "PARTNER": c["partner"], "MECHANISM_CLASS": d.get("mechanism"),
                "RESOLUTION_A": c["resolution"], "METHOD": c["method"],
                "TARGET_UNIPROT": c["target_accession"], "FUSION_UNIPROT": c["fusion_accession"],
                "TOPOLOGY": c["topology"],
                "TARGET_RANGES_ENTITY": str(c["target_ranges_entity"]),
                "FUSION_RANGES_ENTITY": str(c["fusion_ranges_entity"]),
                "N_TARGET_RESIDUES_MODELLED": c["n_target"],
                "N_FUSION_RESIDUES_MODELLED": c["n_fusion"],
                "DETECTOR": det, "RANK1_CLASS": d["rank1_class"],
                "RANK1_SCORE": round(d["rank1_score"], 3),
                "RANK1_F_FUSION": d["f_fusion"], "RANK1_N_RESIDUES": d["n_residues"],
                "PROBE_STATE": d.get("probe_state"),
                "ISOLATED_PARTNER_RANK": d.get("isolated_rank"),
                "TARGET_CAVITY_ORIGINAL_RANK": d.get("target_cavity_original_rank",
                                                     d.get("target_cavity_state")),
                "TARGET_CAVITY_REMOVED_RANK": d.get("target_cavity_removed_rank", "NA"),
                "DISPLACEMENT": d.get("displacement", "NA"),
                "DEPOSITED_GROUPS_IN_POCKET": ";".join(
                    f"{k}x{v}" for k, v in c["deposited_groups_in_pocket"].items()) or "none",
                "ENTRY_TITLE": (c["description"] or "")[:90]})
    with open(os.path.join(FIG, "FIG5_case_validation.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(vcols) + "\n")
        for r in vrows:
            fh.write("\t".join("NA" if r.get(c) is None else str(r.get(c)) for c in vcols) + "\n")
    print(f"  wrote figures/FIG5_case_validation.tsv ({len(vrows)} rows)")

    # ---------------- figure
    fig = plt.figure(figsize=(7.2, 6.4))
    for i, c in enumerate(cases):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        for cl in ("TARGET", "FUSION", "LINKER"):
            p = c["ca"][cl]
            if not p:
                continue
            xs, ys, zs = zip(*p)
            ax.plot(xs, ys, zs, "-", color=CP[cl], lw=0.7, alpha=0.9)
        det = "p2rank" if "p2rank" in c["detectors"] else list(c["detectors"])[0]
        d = c["detectors"][det]
        cx, cy, cz = d["centroid"]
        if not any(math.isnan(v) for v in (cx, cy, cz)):
            ax.scatter([cx], [cy], [cz], s=110, marker="*", color="#111111",
                       depthshade=False, zorder=10)
        if c["ligand_points"]:
            lx, ly, lz = zip(*c["ligand_points"])
            ax.scatter(lx, ly, lz, s=7, color="#E1A100", depthshade=False, zorder=9)
        ax.set_axis_off()
        allp = [q for cl in ("TARGET", "FUSION", "LINKER") for q in c["ca"][cl]]
        if allp:
            xs, ys, zs = zip(*allp)
            cxm, cym, czm = (min(xs)+max(xs))/2, (min(ys)+max(ys))/2, (min(zs)+max(zs))/2
            r = max(max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)) / 2 * 1.02
            ax.set_xlim(cxm-r, cxm+r); ax.set_ylim(cym-r, cym+r); ax.set_zlim(czm-r, czm+r)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        lig = ";".join(f"{k}" for k in c["deposited_groups_in_pocket"]) or "none"
        ax.set_title(
            f"{c['label']}  {c['pdb_id']}  {c['partner']}  {c['resolution']:.2f} $\\AA$\n"
            f"{d['rank1_class']}, $f_{{fusion}}$={d['f_fusion']:.2f}, "
            f"{det} rank 1\n"
            f"target cavity {d.get('target_cavity_original_rank','-')} "
            f"$\\rightarrow$ {d.get('target_cavity_removed_rank','-')} after removal | "
            f"groups: {lig}", fontsize=6.3, pad=-2)
    fig.text(0.02, 0.975,
             "Figure 5  Structural case studies: a real cavity on the wrong biological component",
             fontsize=9, fontweight="bold", ha="left")
    fig.text(0.02, 0.945,
             "blue = target C$\\alpha$ trace · red = fusion partner · grey = linker · "
             "star = rank-1 pocket centroid · gold = deposited groups within 5 $\\AA$.  "
             "Cases selected by a stated deterministic rule; illustrative, not statistically "
             "representative.", fontsize=6.2, ha="left")
    fig.subplots_adjust(left=0.0, right=1.0, top=0.90, bottom=0.0, wspace=0.0, hspace=0.12)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"FIG5_structural_cases.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("  wrote figures/FIG5_structural_cases.png / .pdf")

    # ---------------- PyMOL script for publication-quality rendering
    lines = ["# Figure 5 — publication rendering, generated by scripts/31_figure5_structural.py",
             "# Run:  pymol -cq FIG5_pymol.pml     (structures are fetched from the PDB)",
             "set ray_opaque_background, 0", "bg_color white", "set cartoon_transparency, 0.1", ""]
    for c in cases:
        s, ch = c["pdb_id"], c["chain"]
        tsel = " or ".join(f"resi {a}-{b}" for a, b in c["target_ranges_entity"]) or "none"
        fsel = " or ".join(f"resi {a}-{b}" for a, b in c["fusion_ranges_entity"]) or "none"
        det = "p2rank" if "p2rank" in c["detectors"] else list(c["detectors"])[0]
        pres = "+".join(x for x in c["detectors"][det]["residues"] if x.isdigit()) or "0"
        lines += [
            f"# ---- {c['label']}  {s} chain {ch}  {c['partner']}  {c['class']}",
            f"fetch {s}, async=0", f"hide everything, {s}",
            f"create {c['label']}_chain, {s} and chain {ch} and polymer",
            f"show cartoon, {c['label']}_chain",
            f"color skyblue, {c['label']}_chain",
            f"# NOTE: selections below use ENTITY (label_seq) ranges from SIFTS; map to auth",
            f"#       numbering with the manifest before use if they differ for {s}.",
            f"select {c['label']}_pocket, {c['label']}_chain and resi {pres}",
            f"color firebrick, {c['label']}_pocket",
            f"show sticks, {c['label']}_pocket",
            f"select {c['label']}_lig, {s} and not polymer and not solvent "
            f"within 5 of {c['label']}_pocket",
            f"show spheres, {c['label']}_lig", f"color orange, {c['label']}_lig", ""]
    lines += ["orient", "set ray_trace_mode, 1", "# ray 2400, 1800", "# png FIG5_pymol.png, dpi=400"]
    open(os.path.join(FIG, "FIG5_pymol.pml"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("  wrote figures/FIG5_pymol.pml")

    json.dump([{k: v for k, v in c.items() if k not in ("ca", "ligand_points")} for c in cases],
              open(os.path.join(FIG, "FIG5_cases.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
