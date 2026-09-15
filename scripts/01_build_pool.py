#!/usr/bin/env python
"""
Phase 2B step 1 — enumerate the full candidate pool of chimeric crystallization constructs
and record, for every candidate, whether it passes each protocol v1.1 filter and why.

Outcome-blind: uses only RCSB/SIFTS annotation and entry metadata. No pocket detector is
involved and none has been run at the time this script executes.

Outputs:
  data_manifest/pool_raw.json        one record per chimeric polymer entity, all raw fields
  data_manifest/pool_eligibility.tsv one row per entity with pass/fail + reasons
"""
import json, os, sys, time, urllib.request, urllib.error, hashlib, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
ANNO = os.path.join(ROOT, "annotations")
os.makedirs(MANI, exist_ok=True)

PARTNERS = {"P0ABE7": "BRIL", "P00720": "T4L", "P0AEX9": "MBP"}

# protocol v1.1 section 3.2 / 3.3 thresholds
MAX_RESOLUTION = 3.5
MIN_FUSION_MODELLED = 50
MIN_TARGET_MODELLED = 80
MIN_ALIGN_COVERAGE = 0.80
MAX_SEGMENT_DISORDER = 0.20
ALLOWED_METHODS = {"X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "ELECTRON CRYSTALLOGRAPHY",
                   "NEUTRON DIFFRACTION"}

GQL = "https://data.rcsb.org/graphql"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

QUERY = """
query($ids:[String!]!){
  polymer_entities(entity_ids:$ids){
    rcsb_id
    rcsb_polymer_entity { pdbx_description }
    entity_poly { rcsb_sample_sequence_length }
    rcsb_polymer_entity_align {
      reference_database_accession reference_database_name provenance_source
      aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
    }
    rcsb_polymer_entity_container_identifiers { entry_id auth_asym_ids }
    entry {
      rcsb_id
      exptl { method }
      rcsb_entry_info { resolution_combined }
      struct { title }
      nonpolymer_entities {
        nonpolymer_comp { chem_comp { id name formula_weight } }
        rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
      }
    }
    polymer_entity_instances {
      rcsb_id
      rcsb_polymer_entity_instance_container_identifiers { auth_asym_id asym_id }
      rcsb_polymer_instance_feature { type feature_positions { beg_seq_id end_seq_id } }
    }
  }
}
"""


def post(url, payload, timeout=180, tries=4):
    body = json.dumps(payload).encode()
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=body,
                                         headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=timeout))
        except Exception as exc:
            if attempt == tries - 1:
                raise
            time.sleep(3 * (attempt + 1))
            sys.stderr.write(f"  retry {attempt+1} after {type(exc).__name__}\n")


def search_entities(acc):
    d = post(SEARCH, {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 10000},
                            "results_verbosity": "compact"}})
    return d.get("result_set", [])


def ranges_from_align(align_blocks, accession):
    """Entity-numbering [start, end] inclusive ranges for one UniProt accession."""
    out = []
    for blk in align_blocks or []:
        if blk.get("reference_database_name") != "UniProt":
            continue
        if blk.get("reference_database_accession") != accession:
            continue
        for reg in blk.get("aligned_regions") or []:
            beg, ln = reg.get("entity_beg_seq_id"), reg.get("length")
            if beg is None or not ln:
                continue
            out.append([int(beg), int(beg) + int(ln) - 1])
    return sorted(out)


def expand(ranges):
    s = set()
    for a, b in ranges:
        s.update(range(a, b + 1))
    return s


def fetch_pool():
    cache = os.path.join(MANI, "pool_raw.json")
    if os.path.exists(cache):
        sys.stderr.write("using cached pool_raw.json\n")
        return json.load(open(cache))

    records = []
    for acc, tag in PARTNERS.items():
        ids = search_entities(acc)
        sys.stderr.write(f"{tag} ({acc}): {len(ids)} entities with this accession\n")
        got = []
        for i in range(0, len(ids), 25):
            chunk = ids[i:i + 25]
            d = post(GQL, {"query": QUERY, "variables": {"ids": chunk}})
            got += [r for r in d["data"]["polymer_entities"] if r]
            if (i // 25) % 5 == 0:
                sys.stderr.write(f"  {tag} {i + len(chunk)}/{len(ids)}\n")
        for r in got:
            r["_partner_accession"] = acc
            r["_partner_tag"] = tag
        records += got
    json.dump(records, open(cache, "w"))
    sys.stderr.write(f"cached {len(records)} entity records\n")
    return records


def evaluate(rec):
    """Apply protocol v1.1 3.2/3.3 filters. Returns a dict with verdict and every reason."""
    partner_acc = rec["_partner_accession"]
    out = {
        "entity_id": rec["rcsb_id"],
        "pdb_id": rec["rcsb_id"].split("_")[0],
        "fusion_partner": rec["_partner_tag"],
        "fusion_accession": partner_acc,
        "description": (rec.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
        "title": ((rec.get("entry") or {}).get("struct") or {}).get("title"),
        "exclusions": [],
    }
    align = rec.get("rcsb_polymer_entity_align") or []
    accs = sorted({b["reference_database_accession"] for b in align
                   if b.get("reference_database_name") == "UniProt"})
    out["uniprot_accessions"] = accs
    out["provenance"] = sorted({b.get("provenance_source") for b in align
                                if b.get("reference_database_name") == "UniProt"})

    # --- E1: exactly two UniProt accessions, one of which is the partner
    if len(accs) < 2:
        out["exclusions"].append("E1_not_chimeric_single_uniprot")
    elif len(accs) > 2:
        out["exclusions"].append(f"E1_multi_uniprot_n={len(accs)}")
    if partner_acc not in accs:
        out["exclusions"].append("E1_partner_accession_absent")

    target_acc = None
    if len(accs) == 2 and partner_acc in accs:
        target_acc = [a for a in accs if a != partner_acc][0]
    out["target_accession"] = target_acc

    # --- entry-level metadata
    entry = rec.get("entry") or {}
    methods = [m["method"] for m in (entry.get("exptl") or []) if m.get("method")]
    out["method"] = methods[0] if methods else None
    res = ((entry.get("rcsb_entry_info") or {}).get("resolution_combined") or [])
    out["resolution"] = float(res[0]) if res else None

    if out["method"] not in ALLOWED_METHODS:
        out["exclusions"].append(f"E2_method={out['method']}")
    if out["resolution"] is None:
        out["exclusions"].append("E3_no_resolution")
    elif out["resolution"] > MAX_RESOLUTION:
        out["exclusions"].append(f"E3_resolution={out['resolution']:.2f}>{MAX_RESOLUTION}")

    # --- representative chain: first auth_asym_id in ASCII order
    ids_blk = rec.get("rcsb_polymer_entity_container_identifiers") or {}
    chains = sorted(ids_blk.get("auth_asym_ids") or [])
    out["representative_chain"] = chains[0] if chains else None
    out["all_chains"] = chains
    if not chains:
        out["exclusions"].append("E4_no_chain")

    seq_len = (rec.get("entity_poly") or {}).get("rcsb_sample_sequence_length")
    out["entity_length"] = seq_len

    # --- segments in entity numbering
    fus_r = ranges_from_align(align, partner_acc)
    tgt_r = ranges_from_align(align, target_acc) if target_acc else []
    out["fusion_ranges"] = fus_r
    out["target_ranges"] = tgt_r
    fus_set, tgt_set = expand(fus_r), expand(tgt_r)

    if fus_set & tgt_set:
        out["exclusions"].append(f"E5_overlapping_segments_n={len(fus_set & tgt_set)}")
    if not fus_r:
        out["exclusions"].append("E5_no_fusion_segment")
    if target_acc and not tgt_r:
        out["exclusions"].append("E5_no_target_segment")

    # --- unobserved residues for the representative instance
    unobs = set()
    inst_found = False
    for inst in rec.get("polymer_entity_instances") or []:
        cid = (inst.get("rcsb_polymer_entity_instance_container_identifiers") or {})
        if cid.get("auth_asym_id") != out["representative_chain"]:
            continue
        inst_found = True
        out["instance_id"] = inst.get("rcsb_id")
        for feat in inst.get("rcsb_polymer_instance_feature") or []:
            if feat.get("type") != "UNOBSERVED_RESIDUE_XYZ":
                continue
            for pos in feat.get("feature_positions") or []:
                b, e = pos.get("beg_seq_id"), pos.get("end_seq_id")
                if b is None:
                    continue
                unobs.update(range(int(b), int(e if e is not None else b) + 1))
    if not inst_found:
        out["exclusions"].append("E4_representative_instance_not_found")
    out["n_unobserved"] = len(unobs)

    # --- modelled coverage / disorder
    if seq_len:
        all_pos = set(range(1, int(seq_len) + 1))
        modelled = all_pos - unobs
        aligned = fus_set | tgt_set
        out["n_modelled"] = len(modelled)
        out["n_aligned"] = len(aligned)
        cov = len(aligned & modelled) / len(modelled) if modelled else 0.0
        out["align_coverage_of_modelled"] = round(cov, 4)
        if cov < MIN_ALIGN_COVERAGE:
            out["exclusions"].append(f"E6_align_coverage={cov:.2f}<{MIN_ALIGN_COVERAGE}")

        fus_mod = len(fus_set & modelled)
        tgt_mod = len(tgt_set & modelled)
        out["fusion_modelled"] = fus_mod
        out["target_modelled"] = tgt_mod
        if fus_mod < MIN_FUSION_MODELLED:
            out["exclusions"].append(f"E7_fusion_modelled={fus_mod}<{MIN_FUSION_MODELLED}")
        if tgt_mod < MIN_TARGET_MODELLED:
            out["exclusions"].append(f"E8_target_modelled={tgt_mod}<{MIN_TARGET_MODELLED}")

        fus_dis = (len(fus_set) - fus_mod) / len(fus_set) if fus_set else 1.0
        tgt_dis = (len(tgt_set) - tgt_mod) / len(tgt_set) if tgt_set else 1.0
        out["fusion_disorder"] = round(fus_dis, 4)
        out["target_disorder"] = round(tgt_dis, 4)
        if fus_dis > MAX_SEGMENT_DISORDER:
            out["exclusions"].append(f"E9_fusion_disorder={fus_dis:.2f}>{MAX_SEGMENT_DISORDER}")
        if tgt_dis > MAX_SEGMENT_DISORDER:
            out["exclusions"].append(f"E9_target_disorder={tgt_dis:.2f}>{MAX_SEGMENT_DISORDER}")

        # linker = unaligned positions strictly between the two accessions' territories
        # tag    = unaligned positions outside the outermost aligned position
        if aligned:
            lo, hi = min(aligned), max(aligned)
            out["linker_positions"] = sorted(p for p in all_pos - aligned if lo < p < hi)
            out["tag_positions"] = sorted(p for p in all_pos - aligned if p < lo or p > hi)
        else:
            out["linker_positions"], out["tag_positions"] = [], sorted(all_pos)
    else:
        out["exclusions"].append("E6_no_sequence_length")
        out["linker_positions"], out["tag_positions"] = [], []

    # --- topology
    topo = None
    if fus_r and tgt_r:
        fmin, fmax = min(p for p, _ in fus_r), max(q for _, q in fus_r)
        if len(tgt_r) >= 2 and any(q < fmin for _, q in tgt_r) and any(p > fmax for p, _ in tgt_r):
            topo = "INTERNAL_INSERTION"
        elif fmax < min(p for p, _ in tgt_r):
            topo = "TERMINAL_N"
        elif fmin > max(q for _, q in tgt_r):
            topo = "TERMINAL_C"
        else:
            topo = "COMPLEX"
    out["topology"] = topo

    # --- candidate biological ligands (outcome-independent; artifact filtering happens later)
    ligs = []
    for ne in (entry.get("nonpolymer_entities") or []):
        cc = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if cc.get("id"):
            ligs.append({"id": cc["id"], "name": cc.get("name"),
                         "mw": cc.get("formula_weight"),
                         "chains": (ne.get("rcsb_nonpolymer_entity_container_identifiers")
                                    or {}).get("auth_asym_ids")})
    out["nonpolymer_entities"] = ligs

    out["eligible"] = len(out["exclusions"]) == 0
    return out


def main():
    records = fetch_pool()
    src_hash = hashlib.sha256(
        json.dumps(records, sort_keys=True).encode()).hexdigest()
    sys.stderr.write(f"pool_raw.json sha256 = {src_hash}\n")

    evaluated = [evaluate(r) for r in records]
    json.dump({"source_sha256": src_hash, "n_records": len(evaluated),
               "protocol_version": "1.1", "records": evaluated},
              open(os.path.join(MANI, "pool_eligibility.json"), "w"), indent=1)

    cols = ["entity_id", "pdb_id", "fusion_partner", "target_accession", "eligible",
            "method", "resolution", "topology", "representative_chain", "entity_length",
            "fusion_modelled", "target_modelled", "align_coverage_of_modelled",
            "fusion_disorder", "target_disorder", "description", "exclusions"]
    with open(os.path.join(MANI, "pool_eligibility.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for e in evaluated:
            row = []
            for c in cols:
                v = e.get(c)
                if isinstance(v, list):
                    v = ";".join(map(str, v))
                row.append("" if v is None else str(v).replace("\t", " "))
            fh.write("\t".join(row) + "\n")

    # summary
    print(f"total chimeric-candidate entities: {len(evaluated)}")
    for tag in ("BRIL", "T4L", "MBP"):
        sub = [e for e in evaluated if e["fusion_partner"] == tag]
        ok = [e for e in sub if e["eligible"]]
        tgts = {e["target_accession"] for e in ok}
        print(f"  {tag:5s}: {len(sub):4d} candidates -> {len(ok):4d} eligible "
              f"across {len(tgts)} distinct target accessions")
    reasons = collections.Counter(
        r.split("=")[0] for e in evaluated for r in e["exclusions"])
    print("\nexclusion reasons (entity counts, entities may fail several):")
    for k, v in reasons.most_common():
        print(f"  {k:45s} {v}")


if __name__ == "__main__":
    main()
