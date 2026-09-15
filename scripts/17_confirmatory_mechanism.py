#!/usr/bin/env python
"""
Confirmatory SECONDARY AIM D — mechanism, using the frozen isolated-partner correspondence
method and the frozen category list. No new categories may be created.

States:  MATCHED_INTRINSIC_CAVITY | NO_MATCHING_CAVITY | DETECTOR_NO_PREDICTIONS |
         REFERENCE_STRUCTURE_UNAVAILABLE | AMBIGUOUS
Classes: FUSION_INTRINSIC_CAVITY | CHIMERA_SPECIFIC_FUSION_CAVITY | TARGET_FUSION_INTERFACE |
         LINKER_RELATED | MAPPING_ARTIFACT | UNDETERMINED_<state>
"""
import collections, importlib.util, json, math, os, shutil, subprocess, sys
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = os.environ.get("FUSIONTAG_VARIANT", "primary")
RES = os.path.join(ROOT, "results", "confirmatory", VARIANT)
WORK = os.path.join(r"C:\Users\user\AppData\Local\Temp\claude"
                    r"\C--AI-PROJECTS-Softwares-Fusion-Tag-Hazard"
                    r"\5c4a9bb4-b9dd-40f9-b5cc-b8b58139430b\scratchpad",
                    "conf_mech_" + VARIANT)
os.makedirs(WORK, exist_ok=True)

JACCARD_MIN = 0.40
MAX_CENTROID_DIST = 8.0
MIN_FUSION_RESIDUES = 3


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


runner = _load("runner", "common_runner.py")
cls = _load("cls", "05_classify_pockets.py")


def fusion_coords(rec):
    st = gemmi.read_structure(os.path.join(ROOT, rec["original_pdb"]))
    cmap = rec["auth_class_map"]
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


def write_fusion_only(rec):
    st = gemmi.read_structure(os.path.join(ROOT, rec["original_pdb"]))
    cmap = rec["auth_class_map"]
    out = gemmi.Structure(); out.spacegroup_hm = "P 1"
    m = gemmi.Model("1"); c = gemmi.Chain("A")
    for ch in st[0]:
        for res in ch:
            k = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if cmap.get(k) in ("FUSION", "LINKER"):
                c.add_residue(res)
    m.add_chain(c); out.add_model(m); out.setup_entities()
    p = os.path.join(WORK, f"{rec['pdb_id']}_{rec['auth_chain']}_FUSIONONLY.pdb")
    out.write_pdb(p)
    return p, len(c)


def run_both(inp):
    base = os.path.splitext(os.path.basename(inp))[0]
    p2 = os.path.join(os.path.dirname(inp), base + "_p2rank")
    if not os.path.exists(os.path.join(p2, base + ".pdb_predictions.csv")):
        subprocess.run(runner.p2rank_cmd(inp, p2), capture_output=True, text=True,
                       env=runner.p2rank_env(), timeout=3600)
    fpd = os.path.join(os.path.dirname(inp), base + "_out")
    if not os.path.isdir(fpd):
        argv, env = runner.fpocket_cmd(inp)
        subprocess.run(argv, capture_output=True, text=True, env=env, timeout=3600)
    return (cls.parse_p2rank(p2, base + ".pdb"),
            cls.parse_fpocket(fpd, base) if os.path.isdir(fpd) else [])


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

    failures = []
    for det in ("p2rank", "fpocket"):
        for s in sorted(man):
            pk = idx.get((det, s, "ORIGINAL"))
            if pk and pk[0]["fusion_associated"]:
                failures.append((det, s, pk[0]))
    need = sorted({s for _, s, _ in failures})
    print(f"rank-1 fusion-associated results: {len(failures)} over {len(need)} structures")

    iso = {}
    for i, s in enumerate(need, 1):
        try:
            path, n = write_fusion_only(man[s])
            p2, fp = run_both(path)
            iso[s] = {"p2rank": p2, "fpocket": fp, "n_residues": n}
        except Exception as exc:
            iso[s] = {"error": f"{type(exc).__name__}: {exc}"}
            print(f"  ISO FAIL {s}: {exc}")
        if i % 15 == 0:
            print(f"  isolated-partner runs {i}/{len(need)}")

    rows = []
    for det, s, pocket in failures:
        rec = man[s]
        cmap = rec["auth_class_map"]
        fc = fusion_coords(rec)
        pres = set(pocket["residues"].split(";")) if pocket["residues"] else set()
        fus = {r for r in pres if cmap.get(r) in ("FUSION", "LINKER")}
        info = iso.get(s, {})
        if "error" in info:
            state, bj, br, dist, niso = "REFERENCE_STRUCTURE_UNAVAILABLE", None, None, None, None
        else:
            pk = info.get(det, [])
            niso = len(pk)
            if len(fus) < MIN_FUSION_RESIDUES:
                state, bj, br, dist = "AMBIGUOUS", None, None, None
            elif niso == 0:
                state, bj, br, dist = "DETECTOR_NO_PREDICTIONS", None, None, None
            else:
                cch = centroid(fus, fc)
                bj, br, dist = 0.0, None, None
                for p in pk:
                    ip = set(p["residues"])
                    j = len(fus & ip) / len(fus | ip) if (fus | ip) else 0.0
                    ci = centroid(ip & set(fc), fc)
                    d = math.dist(cch, ci) if (cch and ci) else float("inf")
                    if j > bj:
                        bj, br, dist = j, p["rank"], d
                if bj >= JACCARD_MIN and dist is not None and dist <= MAX_CENTROID_DIST:
                    state = "MATCHED_INTRINSIC_CAVITY"
                elif bj >= JACCARD_MIN or (dist is not None and dist <= MAX_CENTROID_DIST):
                    state = "AMBIGUOUS"
                else:
                    state = "NO_MATCHING_CAVITY"
        unm = 1 - (pocket["n_classifiable"] / max(pocket["n_residues"], 1))
        if unm > 0.5:
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
                     "target": rec["target_accession"], "rank1_class": pocket["pocket_class"],
                     "f_fusion": pocket["f_fusion"],
                     "n_fusion_residues_in_pocket": len(fus),
                     "isolated_n_pockets": niso,
                     "best_jaccard": None if bj is None else round(bj, 3),
                     "best_isolated_rank": br,
                     "centroid_dist": None if dist in (None, float("inf")) else round(dist, 2),
                     "probe_state": state, "mechanism": mech})

    json.dump({"variant": VARIANT, "jaccard_min": JACCARD_MIN,
               "max_centroid_dist": MAX_CENTROID_DIST, "rows": rows},
              open(os.path.join(RES, "mechanism.json"), "w"), indent=1)
    cols = list(rows[0].keys())
    with open(os.path.join(RES, "mechanism.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    print("\n" + "=" * 92)
    print("SECONDARY AIM D — MECHANISM (confirmatory)")
    print("=" * 92)
    print("\nprobe states, per detector:")
    for det in ("p2rank", "fpocket"):
        c = collections.Counter(r["probe_state"] for r in rows if r["detector"] == det)
        print(f"  {det:8s}: {dict(c)}")
    print("\nmechanism, per detector:")
    for det in ("p2rank", "fpocket"):
        c = collections.Counter(r["mechanism"] for r in rows if r["detector"] == det)
        print(f"  {det:8s}: {dict(c)}")
    print("\nmechanism, per detector x partner:")
    for det in ("p2rank", "fpocket"):
        for pt in ("BRIL", "T4L", "MBP"):
            c = collections.Counter(r["mechanism"] for r in rows
                                    if r["detector"] == det and r["partner"] == pt)
            print(f"  {det:8s} {pt:5s}: {dict(c)}")
    print("\nisolated-partner rank of the matched intrinsic cavity:")
    mi = [r for r in rows if r["probe_state"] == "MATCHED_INTRINSIC_CAVITY"]
    r1 = sum(1 for r in mi if r["best_isolated_rank"] == 1)
    print(f"  {r1}/{len(mi)} matched cavities are the isolated partner's OWN rank-1 pocket")


if __name__ == "__main__":
    main()
