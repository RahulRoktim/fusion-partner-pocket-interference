#!/usr/bin/env python
"""
Chain-aware classification, used ONLY for the assembly1 sensitivity variant.

Why this exists. The frozen classifier keys residues by "<resnum><icode>" with no chain, which is
unambiguous for a single-chain input but collides in a multi-chain assembly: residue 167 of chain B
would otherwise inherit the class of residue 167 of chain A. This module parses the chain
identifier from each detector's output and looks up "<chain>:<resnum><icode>".

The frozen single-chain path (scripts/05 and scripts/15) is NOT modified. Equivalence on
single-chain inputs is asserted by tests/test_variant_classifier_equivalence.py.

Writes results/confirmatory/assembly1/pockets_classified.json|tsv, replacing the placeholder that
the frozen driver produced for this variant.
"""
import collections, csv, importlib.util, json, math, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = "assembly1"
RES = os.path.join(ROOT, "results", "confirmatory", VARIANT)


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


cls = _load("cls", "05_classify_pockets.py")


def parse_p2rank_chainaware(outdir, input_basename):
    path = os.path.join(outdir, input_basename + "_predictions.csv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            row = {k.strip(): (v.strip() if isinstance(v, str) else v)
                   for k, v in row.items() if k}
            res = []
            for tok in (row.get("residue_ids") or "").split():
                m = re.match(r"^([^_]*)_(-?\d+)([A-Za-z]?)$", tok)
                if m:
                    res.append(f"{m.group(1)}:{m.group(2)}{m.group(3)}")
            out.append({"rank": int(float(row["rank"])), "score": float(row["score"]),
                        "extra_score": float(row.get("probability") or "nan"),
                        "centroid": (float(row["center_x"]), float(row["center_y"]),
                                     float(row["center_z"])),
                        "residues": res})
    out.sort(key=lambda p: p["rank"])
    return out


def parse_fpocket_chainaware(outdir, base):
    info = os.path.join(outdir, base + "_info.txt")
    pockdir = os.path.join(outdir, "pockets")
    if not os.path.exists(info):
        return []
    scores, cur = {}, None
    for line in open(info, encoding="utf-8", errors="replace"):
        m = re.match(r"^Pocket\s+(\d+)\s*:", line.strip())
        if m:
            cur = int(m.group(1)); scores[cur] = {}
            continue
        if cur is not None and ":" in line:
            k, v = line.split(":", 1)
            try:
                scores[cur][k.strip()] = float(v.strip())
            except ValueError:
                pass
    out = []
    for pid in sorted(scores):
        atm = os.path.join(pockdir, f"pocket{pid}_atm.pdb")
        vert = os.path.join(pockdir, f"pocket{pid}_vert.pqr")
        res, seen = [], set()
        if os.path.exists(atm):
            for line in open(atm, encoding="utf-8", errors="replace"):
                if line.startswith(("ATOM", "HETATM")):
                    key = f"{line[21].strip()}:{line[22:26].strip()}{line[26].strip()}"
                    if key not in seen:
                        seen.add(key); res.append(key)
        cx = cy = cz = 0.0; n = 0
        if os.path.exists(vert):
            for line in open(vert, encoding="utf-8", errors="replace"):
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        cx += float(line[30:38]); cy += float(line[38:46])
                        cz += float(line[46:54]); n += 1
                    except ValueError:
                        pass
        s = scores[pid]
        out.append({"rank": pid, "score": s.get("Score", float("nan")),
                    "extra_score": s.get("Druggability Score", float("nan")),
                    "centroid": (cx/n, cy/n, cz/n) if n else (float("nan"),)*3,
                    "residues": res})
    out.sort(key=lambda p: p["rank"])
    return out


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       f"confirmatory_manifest_{VARIANT}.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    runs = json.load(open(os.path.join(RES, "detector_runs.json")))["runs"]
    rows = []
    for run in runs:
        if not run["success"] or run["pdb_id"] not in man:
            continue
        meta = man[run["pdb_id"]]
        od = os.path.join(ROOT, run["archived_output"])
        base = run["key"]
        pockets = (parse_p2rank_chainaware(od, base + ".pdb") if run["detector"] == "p2rank"
                   else parse_fpocket_chainaware(od, base))
        cmap = meta["auth_class_map"]
        refset = {f"A:{r}" for r in meta["reference_site_residues"]}
        bnds = meta["deletion_boundaries"]
        for p in pockets:
            classes = [cmap.get(r, "UNMAPPED") for r in p["residues"]]
            known = [c for c in classes if c in ("TARGET", "FUSION", "LINKER")]
            n = len(known)
            ft = known.count("TARGET")/n if n else 0.0
            ff = known.count("FUSION")/n if n else 0.0
            fl = known.count("LINKER")/n if n else 0.0
            dmin = None
            if run["condition"] == "REMOVED" and bnds and not math.isnan(p["centroid"][0]):
                dmin = min(math.dist(p["centroid"], (b["x"], b["y"], b["z"])) for b in bnds)
            inter = refset & set(p["residues"])
            rf = len(inter)/len(refset) if refset else None
            row = {"pdb_id": run["pdb_id"], "chain": run["chain"],
                   "fusion_partner": run["partner"], "target_accession": meta["target_accession"],
                   "topology": meta["topology"], "resolution": meta["resolution"],
                   "detector": run["detector"], "condition": run["condition"],
                   "rank": p["rank"], "score": p["score"], "extra_score": p.get("extra_score"),
                   "centroid_x": p["centroid"][0], "centroid_y": p["centroid"][1],
                   "centroid_z": p["centroid"][2],
                   "n_residues": len(p["residues"]), "n_classifiable": n,
                   "n_other_chain_residues": classes.count("OTHER_CHAIN"),
                   "f_target": round(ft, 4), "f_fusion": round(ff, 4), "f_linker": round(fl, 4),
                   "dist_to_deletion_boundary": None if dmin is None else round(dmin, 2),
                   "ref_site_overlap_n": len(inter),
                   "ref_site_overlap_frac": None if rf is None else round(rf, 4),
                   "recovers_reference_site": bool(refset) and len(inter) >= 3 and rf >= 0.25,
                   "residues": ";".join(p["residues"])}
            for t in cls.SENSITIVITY_THRESHOLDS:
                row[f"class_at_{t:.2f}"] = cls.classify(ft, ff, fl, n, t)
            row["pocket_class"] = row["class_at_0.70"]
            row["fusion_associated"] = row["pocket_class"] in ("FUSION_DOMINATED", "INTERFACE",
                                                               "LINKER")
            rows.append(row)
    rows.sort(key=lambda r: (r["detector"], r["pdb_id"], r["condition"], r["rank"]))
    json.dump({"protocol_version": "1.3", "variant": VARIANT, "chain_aware": True,
               "primary_dominance": cls.PRIMARY_DOMINANCE, "pockets": rows},
              open(os.path.join(RES, "pockets_classified.json"), "w"), indent=1)
    cols = [c for c in rows[0].keys() if c != "residues"] + ["residues"]
    with open(os.path.join(RES, "pockets_classified.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"chain-aware classification: {len(rows)} pockets over "
          f"{len({r['pdb_id'] for r in rows})} structures")
    oc = [r for r in rows if r["n_other_chain_residues"] > 0]
    print(f"  pockets touching another chain of the assembly: {len(oc)}/{len(rows)}")


if __name__ == "__main__":
    main()
