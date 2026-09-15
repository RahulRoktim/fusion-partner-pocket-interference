#!/usr/bin/env python
"""
Phase 2F — scientific QC of every rank-1 FUSION_DOMINATED / INTERFACE / LINKER result.

Explanatory only. NO observation is removed from any endpoint on the basis of this QC.

Mechanism is determined with two outcome-independent probes:

  1. FUSION-ONLY condition (QC-only third condition). The fusion+linker residues from the same
     coordinate file are run through both detectors in isolation. If the pocket found on the
     fusion partner inside the chimera is also found when the partner stands alone, the cavity
     is INTRINSIC to the partner. If it is not, the cavity is CHIMERA-SPECIFIC.

  2. Deposited-ligand occupancy. Any deposited heteroatom group sitting inside the offending
     pocket is reported by chemical ID, including groups on the pre-declared artifact list.
     A partner's own native ligand found inside its own predicted pocket is direct evidence
     that the detector recovered the partner's native site.

Mechanism labels:
  FUSION_INTRINSIC_CAVITY, TARGET_FUSION_INTERFACE, LINKER_RELATED,
  MAPPING_ARTIFACT, DETECTOR_ARTIFACT, AMBIGUOUS
"""
import collections, importlib.util, json, os, shutil, subprocess, sys
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")
PREP = os.path.join(ROOT, "prepared")
QCDIR = os.path.join(RESULTS, "qc")
WORK = os.path.join(os.environ.get("TEMP", "/tmp"), "fusiontag_qc")
os.makedirs(QCDIR, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

_s = importlib.util.spec_from_file_location("runner", os.path.join(ROOT, "scripts",
                                                                  "common_runner.py"))
runner = importlib.util.module_from_spec(_s); _s.loader.exec_module(runner)
_s2 = importlib.util.spec_from_file_location("cls", os.path.join(ROOT, "scripts",
                                                                 "05_classify_pockets.py"))
cls = importlib.util.module_from_spec(_s2); _s2.loader.exec_module(cls)

POCKET_LIGAND_DIST = 5.0
OVERLAP_FRACTION = 0.50      # fraction of chimera-pocket fusion residues recovered in isolation


def write_fusion_only(pdb_id, chain, manrec):
    """Write the FUSION+LINKER residues of the prepared ORIGINAL file, alone."""
    src = os.path.join(ROOT, manrec["original_pdb"])
    st = gemmi.read_structure(src)
    cmap = manrec["auth_class_map"]
    out = gemmi.Structure(); out.spacegroup_hm = "P 1"
    m = gemmi.Model("1"); c = gemmi.Chain("A")
    for ch in st[0]:
        for res in ch:
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if cmap.get(key) in ("FUSION", "LINKER"):
                c.add_residue(res)
    m.add_chain(c); out.add_model(m); out.setup_entities()
    path = os.path.join(WORK, f"{pdb_id}_{chain}_FUSIONONLY.pdb")
    out.write_pdb(path)
    return path, len(c)


def run_both(inp):
    base = os.path.splitext(os.path.basename(inp))[0]
    p2dir = os.path.join(os.path.dirname(inp), base + "_p2rank")
    subprocess.run(runner.p2rank_cmd(inp, p2dir), capture_output=True, text=True,
                   env=runner.p2rank_env(), timeout=1800)
    argv, env = runner.fpocket_cmd(inp)
    subprocess.run(argv, capture_output=True, text=True, env=env, timeout=1800)
    fpdir = os.path.join(os.path.dirname(inp), base + "_out")
    return (cls.parse_p2rank(p2dir, base + ".pdb"),
            cls.parse_fpocket(fpdir, base) if os.path.isdir(fpdir) else [])


def ligands_in_pocket(cif_path, chain_id, pocket_res_keys, cmap):
    """Deposited heteroatom groups within POCKET_LIGAND_DIST of the pocket's residues."""
    st = gemmi.read_structure(cif_path); st.setup_entities()
    model = st[0]
    want = set(pocket_res_keys)
    targets = []
    for ch in model:
        if ch.name != chain_id:
            continue
        for res in ch:
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if key in want:
                targets.append(res)
    if not targets:
        return []
    hits = collections.Counter()
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            if res.name in ("HOH", "DOD"):
                continue
            close = 0
            for a in res:
                for t in targets:
                    if any(a.pos.dist(b.pos) <= POCKET_LIGAND_DIST for b in t):
                        close += 1
                        break
                if close:
                    break
            if close:
                hits[res.name] += 1
    return sorted(hits.items(), key=lambda x: -x[1])


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(RESULTS, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    # every rank-1 failure, both detectors
    failures = []
    for det in ("p2rank", "fpocket"):
        for s in sorted(man):
            pk = idx.get((det, s, "ORIGINAL"))
            if pk and pk[0]["fusion_associated"]:
                failures.append((det, s, pk[0]))
    print(f"rank-1 fusion-associated results to inspect: {len(failures)} "
          f"({len({s for _, s, _ in failures})} distinct structures)\n")

    # FUSION-ONLY runs, once per structure that appears in any failure
    need = sorted({s for _, s, _ in failures})
    iso = {}
    for i, s in enumerate(need, 1):
        rec = man[s]
        path, n = write_fusion_only(s, rec["auth_chain"], rec)
        p2, fp = run_both(path)
        iso[s] = {"n_residues": n, "p2rank": p2, "fpocket": fp}
        print(f"  [{i}/{len(need)}] FUSION-ONLY {s}: {n} residues -> "
              f"P2Rank {len(p2)} pockets, fpocket {len(fp)} pockets")

    rows = []
    for det, s, pocket in failures:
        rec = man[s]
        cmap = rec["auth_class_map"]
        pres = set(pocket["residues"].split(";")) if pocket["residues"] else set()
        fus_res = {r for r in pres if cmap.get(r) == "FUSION"}

        # --- probe 1: is this cavity present when the partner stands alone?
        best_ov, best_rank = 0.0, None
        for p in iso[s][det]:
            ov = len(fus_res & set(p["residues"])) / len(fus_res) if fus_res else 0.0
            if ov > best_ov:
                best_ov, best_rank = ov, p["rank"]
        intrinsic = best_ov >= OVERLAP_FRACTION

        # --- probe 2: deposited ligands sitting in the pocket
        cif = os.path.join(ROOT, ".structure_cache", f"{s.lower()}.cif")
        ligs = ligands_in_pocket(cif, rec["auth_chain"], pres, cmap) if os.path.exists(cif) else []

        # --- mechanism
        unmapped_frac = 1 - (pocket["n_classifiable"] / max(pocket["n_residues"], 1))
        if unmapped_frac > 0.5:
            mech = "MAPPING_ARTIFACT"
        elif pocket["pocket_class"] == "LINKER":
            mech = "LINKER_RELATED"
        elif pocket["pocket_class"] == "INTERFACE":
            mech = "TARGET_FUSION_INTERFACE"
        elif pocket["pocket_class"] == "FUSION_DOMINATED" and intrinsic:
            mech = "FUSION_INTRINSIC_CAVITY"
        elif pocket["pocket_class"] == "FUSION_DOMINATED" and not intrinsic:
            mech = "AMBIGUOUS"
        else:
            mech = "AMBIGUOUS"

        rows.append({
            "detector": det, "pdb_id": s, "partner": rec["fusion_partner"],
            "target": rec["target_accession"], "topology": rec["topology"],
            "rank1_class": pocket["pocket_class"], "score": pocket["score"],
            "f_target": pocket["f_target"], "f_fusion": pocket["f_fusion"],
            "f_linker": pocket["f_linker"],
            "n_residues": pocket["n_residues"], "n_classifiable": pocket["n_classifiable"],
            "unmapped_fraction": round(unmapped_frac, 3),
            "isolated_partner_overlap": round(best_ov, 3),
            "isolated_partner_rank": best_rank,
            "cavity_intrinsic_to_partner": intrinsic,
            "deposited_groups_in_pocket": ";".join(f"{k}x{v}" for k, v in ligs[:6]),
            "mechanism": mech,
        })

    json.dump({"protocol_version": "1.1", "qc_only": True,
               "overlap_fraction_criterion": OVERLAP_FRACTION,
               "rows": rows}, open(os.path.join(QCDIR, "qc_mechanism.json"), "w"), indent=1)
    cols = list(rows[0].keys())
    with open(os.path.join(QCDIR, "qc_mechanism.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("\n" + "=" * 110)
    print("MECHANISM OF EVERY RANK-1 FUSION-ASSOCIATED RESULT")
    print("=" * 110)
    print(f"{'det':8s} {'pdb':6s} {'ptnr':5s} {'class':17s} {'fFus':>5s} "
          f"{'isoOv':>6s} {'isoRk':>5s} {'mechanism':24s} groups in pocket")
    for r in rows:
        print(f"{r['detector']:8s} {r['pdb_id']:6s} {r['partner']:5s} {r['rank1_class']:17s} "
              f"{r['f_fusion']:5.2f} {r['isolated_partner_overlap']:6.2f} "
              f"{str(r['isolated_partner_rank']):>5s} {r['mechanism']:24s} "
              f"{r['deposited_groups_in_pocket'][:40]}")
    print("\nmechanism counts:", dict(collections.Counter(r["mechanism"] for r in rows)))
    print("by partner:", dict(collections.Counter(
        (r["partner"], r["mechanism"]) for r in rows)))


if __name__ == "__main__":
    main()
