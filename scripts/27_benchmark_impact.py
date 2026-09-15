#!/usr/bin/env python
"""
Phase 5 STEP 9 — direct benchmark impact for chimeric benchmark members.

Runs the frozen detector versions on the EXACT benchmark input files (as distributed in
p2rank-datasets), not on reconstructions, and asks where the benchmark's ground-truth site ranks
and what that pocket is made of.

This is a dataset-audit analysis. It does not touch and cannot alter the confirmatory study.

NOTE ON SAMPLE SIZE: Tier-1 datasets contain exactly ONE recognised-partner chimeric member
(1DUG in JOINED/DT198). n = 1 supports description only, never a metric claim.
"""
import collections, importlib.util, json, os, subprocess, sys
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "dataset_audit")
MEM = os.path.join(AUD, "membership")
WORK = os.path.join(os.environ.get("TEMP", "/tmp"), "fusiontag_benchmark_impact")
os.makedirs(WORK, exist_ok=True)
CONTACT = 4.0


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


runner = _load("runner", "common_runner.py")
cls = _load("cls", "05_classify_pockets.py")

CASES = [{
    "dataset": "JOINED", "subset": "dt198", "role": "DEVELOPMENT/VALIDATION",
    "pdb_id": "1DUG", "chain": "A",
    "benchmark_file": os.path.join(MEM, "1dug_benchmark.pdb"),
    "ground_truth_ligands": ["GSH"],
    "fusion_accession": "P08515", "fusion_partner": "GST-Sj26",
    "target_accession": "P02679",
}]


def residue_classes(pdb_id, chain, fus_acc, tgt_acc):
    """auth residue key -> TARGET / FUSION / UNMAPPED, from SIFTS via the audit index."""
    ents = {r["pdb_id"]: r for r in
            json.load(open(os.path.join(AUD, "entry_constructs.json")))["entries"]}
    rec = ents[pdb_id]
    ent = next(e for e in rec["entities"]
               if str(e.get("construct_class", "")).startswith("RECOGNISED"))
    tgt, fus = set(), set()
    for a, b in ent["target_ranges"] or []:
        tgt.update(range(a, b + 1))
    for a, b in ent["fusion_ranges"] or []:
        fus.update(range(a, b + 1))
    st = gemmi.read_structure(os.path.join(ROOT, ".structure_cache", f"{pdb_id.lower()}.cif"))
    st.setup_entities()
    cmap = {}
    for ch in st[0]:
        if ch.name != chain:
            continue
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid() or res.label_seq is None:
                continue
            ls = int(res.label_seq)
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            cmap[key] = "TARGET" if ls in tgt else ("FUSION" if ls in fus else "UNMAPPED")
    return cmap, ent


def gt_site_residues(path, lig_codes, cmap):
    """Residues of the benchmark input within CONTACT of the ground-truth ligand."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    model = st[0]
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    site, comp = set(), collections.Counter()
    for ch in model:
        for res in ch:
            if res.name not in lig_codes:
                continue
            for atom in res:
                for m in ns.find_atoms(atom.pos, '\0', radius=CONTACT):
                    cra = m.to_cra(model)
                    ri = gemmi.find_tabulated_residue(cra.residue.name)
                    if ri is None or not ri.is_amino_acid():
                        continue
                    if atom.pos.dist(cra.atom.pos) > CONTACT:
                        continue
                    k = f"{cra.residue.seqid.num}{cra.residue.seqid.icode.strip()}"
                    site.add(k)
                    comp[cmap.get(k, "UNMAPPED")] += 1
    return site, comp


def main():
    print("=" * 96)
    print("STEP 9 — DIRECT BENCHMARK IMPACT (dataset audit; n is tiny, description only)")
    print("=" * 96)
    out = []
    for c in CASES:
        print(f"\n{c['pdb_id']} chain {c['chain']} — {c['dataset']}/{c['subset']} "
              f"({c['role']}), ground-truth ligand {c['ground_truth_ligands']}")
        cmap, ent = residue_classes(c["pdb_id"], c["chain"],
                                    c["fusion_accession"], c["target_accession"])
        nT = sum(1 for v in cmap.values() if v == "TARGET")
        nF = sum(1 for v in cmap.values() if v == "FUSION")
        print(f"  construct: {ent['description']}")
        print(f"  chain composition: TARGET {nT} residues, FUSION({c['fusion_partner']}) "
              f"{nF} residues, topology {ent.get('topology')}")

        site, comp = gt_site_residues(c["benchmark_file"], set(c["ground_truth_ligands"]), cmap)
        print(f"  ground-truth site: {len(site)} residues, composition {dict(comp)}")
        gt_on = ("FUSION" if comp.get("FUSION", 0) > comp.get("TARGET", 0) else "TARGET")
        print(f"  -> the benchmark's ground-truth positive lies on the {gt_on}")

        base = f"{c['pdb_id']}_{c['chain']}_BENCHMARK"
        inp = os.path.join(WORK, base + ".pdb")
        import shutil
        shutil.copyfile(c["benchmark_file"], inp)
        p2dir = os.path.join(WORK, base + "_p2rank")
        subprocess.run(runner.p2rank_cmd(inp, p2dir), capture_output=True, text=True,
                       env=runner.p2rank_env(), timeout=1800)
        argv, env = runner.fpocket_cmd(inp)
        subprocess.run(argv, capture_output=True, text=True, env=env, timeout=1800)
        fpdir = os.path.join(WORK, base + "_out")
        preds = {"p2rank": cls.parse_p2rank(p2dir, base + ".pdb"),
                 "fpocket": cls.parse_fpocket(fpdir, base) if os.path.isdir(fpdir) else []}

        res = {"case": c["pdb_id"], "dataset": c["dataset"], "role": c["role"],
               "ground_truth_on": gt_on, "gt_site_size": len(site),
               "gt_site_composition": dict(comp), "detectors": {}}
        for det, pk in preds.items():
            hit_rank, hit_frac, hit_comp = None, None, None
            for p in pk:
                inter = site & set(p["residues"])
                frac = len(inter) / len(site) if site else 0
                if frac >= 0.25 and len(inter) >= 3:
                    hit_rank = p["rank"]
                    hit_frac = round(frac, 3)
                    kn = [cmap.get(r, "UNMAPPED") for r in p["residues"]]
                    hit_comp = {"TARGET": kn.count("TARGET"), "FUSION": kn.count("FUSION"),
                                "UNMAPPED": kn.count("UNMAPPED")}
                    break
            top = pk[0] if pk else None
            topk = ([cls.classify(
                0, 0, 0, 0, 0.7)] if False else None)
            topcomp = None
            if top:
                kn = [cmap.get(r, "UNMAPPED") for r in top["residues"]]
                topcomp = {"TARGET": kn.count("TARGET"), "FUSION": kn.count("FUSION"),
                           "UNMAPPED": kn.count("UNMAPPED")}
            print(f"  {det:8s}: {len(pk)} pockets | ground-truth site recovered at rank "
                  f"{hit_rank} (overlap {hit_frac}) | rank-1 pocket composition {topcomp}")
            res["detectors"][det] = {
                "n_pockets": len(pk), "gt_rank": hit_rank, "gt_overlap": hit_frac,
                "gt_pocket_composition": hit_comp, "rank1_composition": topcomp,
                "top1_success": hit_rank == 1,
                "top3_success": hit_rank is not None and hit_rank <= 3,
                "topn_plus2_success": hit_rank is not None and hit_rank <= 3}
        out.append(res)

    json.dump({"note": "Tier-1 datasets contain exactly one recognised-partner chimeric member; "
                       "n = 1 supports description only.", "cases": out},
              open(os.path.join(AUD, "benchmark_impact.json"), "w"), indent=1)

    print("\n" + "=" * 96)
    print("INTERPRETATION")
    print("=" * 96)
    for r in out:
        for det, d in r["detectors"].items():
            if r["ground_truth_on"] == "FUSION" and d["top1_success"]:
                print(f"  {r['case']} / {det}: the method scores a TOP-1 SUCCESS on this "
                      f"benchmark member by predicting a cavity that belongs to the "
                      f"crystallization carrier, not to the biological target.")
            elif r["ground_truth_on"] == "FUSION":
                print(f"  {r['case']} / {det}: ground-truth (carrier) site recovered at rank "
                      f"{d['gt_rank']}; top-1 credit {d['top1_success']}.")


if __name__ == "__main__":
    main()
