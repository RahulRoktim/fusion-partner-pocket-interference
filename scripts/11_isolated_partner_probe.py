#!/usr/bin/env python
"""
Phase 3A.3 — revised isolated-partner mechanism probe.

The pilot probe conflated "the isolated partner genuinely lacks this cavity" with "the detector
returned nothing". This version makes the states explicit and applies the same deterministic
correspondence machinery used for target cavities, now on FUSION residues.

States:
  MATCHED_INTRINSIC_CAVITY        a corresponding cavity exists on the isolated partner
  NO_MATCHING_CAVITY              the isolated run DID predict pockets, none corresponds
  DETECTOR_NO_PREDICTIONS         the isolated run predicted nothing -> uninformative, NOT absence
  REFERENCE_STRUCTURE_UNAVAILABLE no isolated reference structure could be built
  AMBIGUOUS                       too few fusion residues to match, or borderline correspondence

Correspondence uses the SAME frozen rule as target matching (Jaccard >= 0.40 on fusion residues
AND fusion-centroid distance <= 8 A). Fusion coordinates are identical between the chimera and the
isolated extraction, so centroid distance is directly comparable.
"""
import collections, importlib.util, json, math, os, subprocess
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")
QCDIR = os.path.join(RESULTS, "qc")
WORK = os.path.join(os.environ.get("TEMP", "/tmp"), "fusiontag_qc")

JACCARD_MIN = 0.40
MAX_CENTROID_DIST = 8.0
MIN_FUSION_RESIDUES = 3


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m


runner = _load("runner", "common_runner.py")
cls = _load("cls", "05_classify_pockets.py")


def fusion_coords(manrec):
    st = gemmi.read_structure(os.path.join(ROOT, manrec["original_pdb"]))
    cmap = manrec["auth_class_map"]
    out = {}
    for ch in st[0]:
        for res in ch:
            k = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if cmap.get(k) in ("FUSION", "LINKER"):
                out[k] = [(a.pos.x, a.pos.y, a.pos.z) for a in res]
    return out


def centroid(keys, coords):
    pts = [p for k in keys for p in coords.get(k, [])]
    if not pts:
        return None
    n = len(pts)
    return (sum(p[0] for p in pts)/n, sum(p[1] for p in pts)/n, sum(p[2] for p in pts)/n)


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

    # all rank-1 fusion-associated results, both detectors
    failures = []
    for det in ("p2rank", "fpocket"):
        for s in sorted(man):
            pk = idx.get((det, s, "ORIGINAL"))
            if pk and pk[0]["fusion_associated"]:
                failures.append((det, s, pk[0]))

    rows = []
    for det, s, pocket in failures:
        rec = man[s]
        cmap = rec["auth_class_map"]
        fcoords = fusion_coords(rec)
        pres = set(pocket["residues"].split(";")) if pocket["residues"] else set()
        fus_res = {r for r in pres if cmap.get(r) in ("FUSION", "LINKER")}

        base = f"{s}_{rec['auth_chain']}_FUSIONONLY"
        iso_path = os.path.join(WORK, base + ".pdb")
        if not os.path.exists(iso_path):
            state, best_ov, best_rank, n_iso, dist = "REFERENCE_STRUCTURE_UNAVAILABLE", None, None, None, None
        else:
            if det == "p2rank":
                iso = cls.parse_p2rank(os.path.join(WORK, base + "_p2rank"), base + ".pdb")
            else:
                d = os.path.join(WORK, base + "_out")
                iso = cls.parse_fpocket(d, base) if os.path.isdir(d) else []
            n_iso = len(iso)
            if len(fus_res) < MIN_FUSION_RESIDUES:
                state, best_ov, best_rank, dist = "AMBIGUOUS", None, None, None
            elif n_iso == 0:
                state, best_ov, best_rank, dist = "DETECTOR_NO_PREDICTIONS", None, None, None
            else:
                c_chim = centroid(fus_res, fcoords)
                best_ov, best_rank, dist = 0.0, None, None
                for p in iso:
                    ip = set(p["residues"])
                    j = len(fus_res & ip) / len(fus_res | ip) if (fus_res | ip) else 0.0
                    c_iso = centroid(ip & set(fcoords), fcoords)
                    d = math.dist(c_chim, c_iso) if (c_chim and c_iso) else float("inf")
                    if j > best_ov:
                        best_ov, best_rank, dist = j, p["rank"], d
                if best_ov >= JACCARD_MIN and dist is not None and dist <= MAX_CENTROID_DIST:
                    state = "MATCHED_INTRINSIC_CAVITY"
                elif best_ov >= JACCARD_MIN or (dist is not None and dist <= MAX_CENTROID_DIST):
                    state = "AMBIGUOUS"
                else:
                    state = "NO_MATCHING_CAVITY"

        # mechanism, now conditioned on a state that cannot silently mean "absent"
        unmapped = 1 - (pocket["n_classifiable"] / max(pocket["n_residues"], 1))
        if unmapped > 0.5:
            mech = "MAPPING_ARTIFACT"
        elif pocket["pocket_class"] == "LINKER":
            mech = "LINKER_RELATED"
        elif pocket["pocket_class"] == "INTERFACE":
            mech = "TARGET_FUSION_INTERFACE"
        elif state == "MATCHED_INTRINSIC_CAVITY":
            mech = "FUSION_INTRINSIC_CAVITY"
        elif state == "NO_MATCHING_CAVITY":
            mech = "CHIMERA_SPECIFIC_FUSION_CAVITY"
        else:
            mech = "UNDETERMINED_" + state

        rows.append({"detector": det, "pdb_id": s, "partner": rec["fusion_partner"],
                     "rank1_class": pocket["pocket_class"],
                     "f_fusion": pocket["f_fusion"],
                     "n_fusion_residues_in_pocket": len(fus_res),
                     "isolated_n_pockets": n_iso,
                     "best_jaccard": None if best_ov is None else round(best_ov, 3),
                     "best_isolated_rank": best_rank,
                     "centroid_dist": None if dist in (None, float("inf")) else round(dist, 2),
                     "probe_state": state, "mechanism": mech})

    json.dump({"jaccard_min": JACCARD_MIN, "max_centroid_dist": MAX_CENTROID_DIST,
               "rows": rows}, open(os.path.join(QCDIR, "isolated_probe_revised.json"), "w"),
              indent=1)
    cols = list(rows[0].keys())
    with open(os.path.join(QCDIR, "isolated_probe_revised.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("=" * 108)
    print("REVISED ISOLATED-PARTNER PROBE")
    print("=" * 108)
    print(f"{'det':8s} {'pdb':6s} {'ptnr':5s} {'class':17s} {'isoN':>5s} {'J':>6s} "
          f"{'isoRk':>5s} {'dist':>6s} {'state':30s} mechanism")
    for r in rows:
        print(f"{r['detector']:8s} {r['pdb_id']:6s} {r['partner']:5s} {r['rank1_class']:17s} "
              f"{str(r['isolated_n_pockets']):>5s} {str(r['best_jaccard']):>6s} "
              f"{str(r['best_isolated_rank']):>5s} {str(r['centroid_dist']):>6s} "
              f"{r['probe_state']:30s} {r['mechanism']}")
    print("\nprobe states:", dict(collections.Counter(r["probe_state"] for r in rows)))
    print("mechanisms  :", dict(collections.Counter(r["mechanism"] for r in rows)))
    print("\ncomparison with the pilot QC labelling:")
    old = json.load(open(os.path.join(QCDIR, "qc_mechanism.json")))["rows"]
    oldc = collections.Counter(r["mechanism"] for r in old)
    newc = collections.Counter(r["mechanism"] for r in rows)
    for k in sorted(set(oldc) | set(newc)):
        print(f"  {k:36s} pilot={oldc.get(k,0):3d}  revised={newc.get(k,0):3d}")


if __name__ == "__main__":
    main()
