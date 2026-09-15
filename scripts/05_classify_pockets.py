#!/usr/bin/env python
"""
Phase 2D — parse detector output and apply the frozen protocol v1.1 pocket classes.

For every predicted pocket: detector, condition, rank, score, centroid, associated residues,
f_target, f_fusion, f_linker, class (at 0.50 / 0.70-primary / 0.90 dominance thresholds),
distance to nearest fusion-deletion boundary, and overlap with the outcome-independent
target reference site.

Outputs:
  results/pockets_classified.tsv    one row per predicted pocket
  results/pockets_classified.json
"""
import csv, json, math, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")

PRIMARY_DOMINANCE = 0.70          # protocol 6
INTERFACE_FLOOR = 0.20
LINKER_DOMINANCE = 0.50
MIN_CLASSIFIABLE = 5
SENSITIVITY_THRESHOLDS = (0.50, 0.70, 0.90)

# reference-site recovery criterion (pre-declared; DEVIATIONS B4)
REF_OVERLAP_FRACTION = 0.25
REF_OVERLAP_MIN_RESIDUES = 3


def classify(f_t, f_f, f_l, n_classifiable, dom):
    if n_classifiable < MIN_CLASSIFIABLE:
        return "UNASSIGNED"
    if f_l >= LINKER_DOMINANCE:
        return "LINKER"
    if f_f >= dom:
        return "FUSION_DOMINATED"
    if f_t >= dom:
        return "TARGET_DOMINATED"
    if f_t >= INTERFACE_FLOOR and f_f >= INTERFACE_FLOOR:
        return "INTERFACE"
    return "MIXED"


# --------------------------------------------------------------------------- parsers
def parse_p2rank(outdir, input_basename):
    path = os.path.join(outdir, input_basename + "_predictions.csv")
    if not os.path.exists(path):
        return []
    pockets = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            row = {k.strip(): (v.strip() if isinstance(v, str) else v)
                   for k, v in row.items() if k}
            res = []
            for tok in (row.get("residue_ids") or "").split():
                m = re.match(r"^([^_]*)_(-?\d+)([A-Za-z]?)$", tok)
                if m:
                    res.append(f"{m.group(2)}{m.group(3)}")
            pockets.append({
                "rank": int(float(row["rank"])),
                "score": float(row["score"]),
                "extra_score": float(row.get("probability") or "nan"),
                "centroid": (float(row["center_x"]), float(row["center_y"]),
                             float(row["center_z"])),
                "residues": res,
            })
    pockets.sort(key=lambda p: p["rank"])
    return pockets


def parse_fpocket(outdir, base):
    info = os.path.join(outdir, base + "_info.txt")
    pockdir = os.path.join(outdir, "pockets")
    if not os.path.exists(info):
        return []
    scores = {}
    cur = None
    for line in open(info, encoding="utf-8", errors="replace"):
        m = re.match(r"^Pocket\s+(\d+)\s*:", line.strip())
        if m:
            cur = int(m.group(1))
            scores[cur] = {}
            continue
        if cur is not None and ":" in line:
            k, v = line.split(":", 1)
            try:
                scores[cur][k.strip()] = float(v.strip())
            except ValueError:
                pass
    pockets = []
    for pid in sorted(scores):
        # fpocket 4.2.3 names pocket files 1-based, matching "Pocket N :" in the info file.
        # (Validated in tests/test_detector_parsers.py by alpha-sphere-count cross-check.)
        atm = os.path.join(pockdir, f"pocket{pid}_atm.pdb")
        vert = os.path.join(pockdir, f"pocket{pid}_vert.pqr")
        res = []
        if os.path.exists(atm):
            seen = set()
            for line in open(atm, encoding="utf-8", errors="replace"):
                if line.startswith(("ATOM", "HETATM")):
                    key = f"{line[22:26].strip()}{line[26].strip()}"
                    if key not in seen:
                        seen.add(key)
                        res.append(key)
        cx = cy = cz = 0.0
        n = 0
        if os.path.exists(vert):
            for line in open(vert, encoding="utf-8", errors="replace"):
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        cx += float(line[30:38]); cy += float(line[38:46])
                        cz += float(line[46:54]); n += 1
                    except ValueError:
                        pass
        centroid = (cx / n, cy / n, cz / n) if n else (float("nan"),) * 3
        s = scores[pid]
        pockets.append({
            "rank": pid,
            "score": s.get("Score", float("nan")),
            "extra_score": s.get("Druggability Score", float("nan")),
            "centroid": centroid,
            "residues": res,
            "n_alpha_spheres": s.get("Number of Alpha Spheres"),
            "volume": s.get("Volume"),
        })
    pockets.sort(key=lambda p: p["rank"])
    return pockets


# --------------------------------------------------------------------------- main
def main():
    man = {p["pdb_id"] + "_" + p.get("auth_chain", ""): p
           for p in json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    runs = json.load(open(os.path.join(RESULTS, "detector_runs.json")))["runs"]

    rows = []
    for run in runs:
        if not run["success"]:
            continue
        skey = f"{run['pdb_id']}_{run['chain']}"
        meta = man[skey]
        outdir = os.path.join(ROOT, run["archived_output"])
        base = run["key"]
        if run["detector"] == "p2rank":
            pockets = parse_p2rank(outdir, base + ".pdb")
        else:
            pockets = parse_fpocket(outdir, base)

        cmap = meta["auth_class_map"]
        refset = set(meta["reference_site_residues"])
        bnds = meta["deletion_boundaries"]

        for p in pockets:
            classes = [cmap.get(r, "UNMAPPED") for r in p["residues"]]
            known = [c for c in classes if c in ("TARGET", "FUSION", "LINKER")]
            n = len(known)
            f_t = known.count("TARGET") / n if n else 0.0
            f_f = known.count("FUSION") / n if n else 0.0
            f_l = known.count("LINKER") / n if n else 0.0

            dmin = None
            if run["condition"] == "REMOVED" and bnds and not math.isnan(p["centroid"][0]):
                dmin = min(math.dist(p["centroid"], (b["x"], b["y"], b["z"])) for b in bnds)

            inter = refset & set(p["residues"])
            ref_frac = len(inter) / len(refset) if refset else None
            recovered = bool(refset) and len(inter) >= REF_OVERLAP_MIN_RESIDUES and \
                ref_frac >= REF_OVERLAP_FRACTION

            row = {
                "pdb_id": run["pdb_id"], "chain": run["chain"],
                "fusion_partner": run["partner"], "target_accession": meta["target_accession"],
                "topology": meta["topology"], "resolution": meta["resolution"],
                "detector": run["detector"], "condition": run["condition"],
                "rank": p["rank"], "score": p["score"], "extra_score": p.get("extra_score"),
                "centroid_x": p["centroid"][0], "centroid_y": p["centroid"][1],
                "centroid_z": p["centroid"][2],
                "n_residues": len(p["residues"]), "n_classifiable": n,
                "f_target": round(f_t, 4), "f_fusion": round(f_f, 4), "f_linker": round(f_l, 4),
                "dist_to_deletion_boundary": None if dmin is None else round(dmin, 2),
                "ref_site_overlap_n": len(inter),
                "ref_site_overlap_frac": None if ref_frac is None else round(ref_frac, 4),
                "recovers_reference_site": recovered,
                "residues": ";".join(p["residues"]),
            }
            for t in SENSITIVITY_THRESHOLDS:
                row[f"class_at_{t:.2f}"] = classify(f_t, f_f, f_l, n, t)
            row["pocket_class"] = row[f"class_at_{PRIMARY_DOMINANCE:.2f}"]
            row["fusion_associated"] = row["pocket_class"] in (
                "FUSION_DOMINATED", "INTERFACE", "LINKER")
            rows.append(row)

    rows.sort(key=lambda r: (r["detector"], r["pdb_id"], r["condition"], r["rank"]))
    json.dump({"protocol_version": "1.1",
               "primary_dominance": PRIMARY_DOMINANCE,
               "interface_floor": INTERFACE_FLOOR,
               "ref_overlap_criterion": {"fraction": REF_OVERLAP_FRACTION,
                                         "min_residues": REF_OVERLAP_MIN_RESIDUES},
               "pockets": rows},
              open(os.path.join(RESULTS, "pockets_classified.json"), "w"), indent=1)

    cols = [c for c in rows[0].keys() if c != "residues"] + ["residues"] if rows else []
    with open(os.path.join(RESULTS, "pockets_classified.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    print(f"classified {len(rows)} pockets")
    for det in ("p2rank", "fpocket"):
        for cond in ("ORIGINAL", "REMOVED"):
            sub = [r for r in rows if r["detector"] == det and r["condition"] == cond]
            structs = len({r["pdb_id"] for r in sub})
            print(f"  {det:8s} {cond:9s}: {len(sub):5d} pockets over {structs} structures "
                  f"(median {len(sub)/max(structs,1):.1f}/structure)")


if __name__ == "__main__":
    main()
