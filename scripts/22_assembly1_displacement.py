#!/usr/bin/env python
"""
[G4] Chain-aware F3 displacement for the assembly1 variant ONLY.

Bug being fixed: scripts/16 computes correspondence with the frozen
`corr.load_target_coords`, which keys residues as "<resnum><icode>". For the assembly1 variant the
pocket residue strings and the class map are chain-qualified ("A:167"), so the two key spaces never
matched, every centroid came back None, and all 126 structures were misreported as
TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL. That is an artifact of key formatting, not a result.

Fix: a chain-aware coordinate loader, monkeypatched into the frozen correspondence module for this
variant only. The correspondence ALGORITHM, its thresholds and the F3 endpoint definition are
unchanged; only the key format changes. The primary and ions variants are untouched.
"""
import collections, importlib.util, json, math, os
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = "assembly1"
RES = os.path.join(ROOT, "results", "confirmatory", VARIANT)


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


corr = _load("corr", "08_pocket_correspondence.py")
ana = _load("ana", "16_confirmatory_analysis.py")


def load_target_coords_chainaware(manrec):
    """Chain-qualified 'CHAIN:resnum' -> heavy-atom coordinates, TARGET residues only."""
    st = gemmi.read_structure(os.path.join(ROOT, manrec["original_pdb"]))
    cmap = manrec["auth_class_map"]
    out = {}
    for ch in st[0]:
        for res in ch:
            key = f"{ch.name}:{res.seqid.num}{res.seqid.icode.strip()}"
            if cmap.get(key) == "TARGET":
                out[key] = [(a.pos.x, a.pos.y, a.pos.z) for a in res]
    return out


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       f"confirmatory_manifest_{VARIANT}.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(RES, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    structures = sorted(man)
    coords = {s: load_target_coords_chainaware(man[s]) for s in structures}
    nz = sum(1 for s in structures if coords[s])
    print(f"chain-aware target coordinates built for {nz}/{len(structures)} structures")

    out, rows = {}, []
    for det in ("p2rank", "fpocket"):
        recs = []
        for s in structures:
            cmap = man[s]["auth_class_map"]
            orig = idx.get((det, s, "ORIGINAL"), [])
            rem = idx.get((det, s, "REMOVED"), [])
            if not rem:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": man[s]["target_accession"],
                             "state": "DETECTOR_ABSTENTION", "displacement": None})
                continue
            top = rem[0]
            if len(corr.target_set(top, cmap)) < corr.MIN_TARGET_RESIDUES:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": man[s]["target_accession"],
                             "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                             "displacement": None})
                continue
            cand, O, R = corr.build_pairs(orig, rem, coords[s], cmap)
            pairs = corr.assign(cand, O, R, corr.PRIMARY_JACCARD, corr.PRIMARY_MAX_CENTROID_DIST)
            jrem = next((j for j, r in enumerate(R) if r["rank"] == top["rank"]), None)
            hit = next(((i, j, m) for i, j, m in pairs if j == jrem), None)
            if hit is None:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": man[s]["target_accession"],
                             "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                             "displacement": None})
                continue
            i, j, m = hit
            recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                         "target": man[s]["target_accession"], "state": "EVALUABLE",
                         "removed_rank": top["rank"], "original_rank": O[i]["rank"],
                         "displacement": O[i]["rank"] - top["rank"],
                         "jaccard": round(m["jaccard"], 3),
                         "centroid_dist": round(m["dist"], 2)})
        for r in recs:
            r["detector"] = det
        rows += recs
        ev = [r for r in recs if r["state"] == "EVALUABLE"]
        d = [r["displacement"] for r in ev]
        k = sum(1 for x in d if x > 0)
        lo, hi = ana.cluster_bootstrap([(r["target"], 1 if r["displacement"] > 0 else 0)
                                        for r in ev], ana.mean)
        print(f"\n{det.upper()}")
        print(f"  evaluable                            {len(ev)}")
        print(f"  DETECTOR_ABSTENTION                  "
              f"{sum(1 for r in recs if r['state']=='DETECTOR_ABSTENTION')}")
        print(f"  TARGET_CAVITY_NOT_RECOVERED_IN_ORIG  "
              f"{sum(1 for r in recs if r['state']=='TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL')}")
        if ev:
            print(f"  median displacement                  {ana.median(d):g}")
            print(f"  mean displacement                    {ana.mean(d):+.2f}")
            print(f"  displacement >0                      {k}/{len(ev)} = "
                  f"{100*k/len(ev):.1f}%  boot[{100*lo:.1f}-{100*hi:.1f}]")
            print(f"  displacement <0                      {sum(1 for x in d if x < 0)}")
            for pt in ("BRIL", "T4L", "MBP"):
                sub = [r for r in ev if r["partner"] == pt]
                if sub:
                    dd = [r["displacement"] for r in sub]
                    print(f"    {pt:5s} n={len(sub):3d} median {ana.median(dd):g} "
                          f">0 in {sum(1 for x in dd if x>0)}/{len(sub)}")
            out[det] = {"evaluable": len(ev), "median": ana.median(d), "mean": ana.mean(d),
                        "k_gt0": k, "n": len(ev), "boot": (lo, hi),
                        "n_negative": sum(1 for x in d if x < 0)}
    json.dump({"variant": VARIANT, "chain_aware": True, "summary": out, "rows": rows},
              open(os.path.join(RES, "displacement_chainaware.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
