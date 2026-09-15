#!/usr/bin/env python
"""
Phase 5 STEP 2 — identify chimeric crystallization constructs among dataset members.

Uses the same SIFTS-based construct detection as the main study: a polymer entity is chimeric
when SIFTS maps it to two or more distinct UniProt accessions. A STANDALONE structure of
MBP / T4L / BRIL is NOT a fusion construct and is counted separately.

Recognition of the fusion partner is by UniProt accession. Chimeras whose partner is not in the
recognised crystallization-partner list are reported as OTHER_CHIMERA, not silently dropped.

Outputs: dataset_audit/entry_constructs.json  (one record per PDB entry)
"""
import collections, json, os, sys, time, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "dataset_audit")
GQL = "https://data.rcsb.org/graphql"
CHUNK = 40

# Recognised crystallization / solubility fusion partners, by UniProt accession.
PARTNERS = {
    "P0ABE7": "BRIL/cyt-b562", "P00720": "T4L", "P0AEX9": "MBP",
    "P00268": "rubredoxin", "P0AA25": "thioredoxin", "P08515": "GST-Sj26",
    "P0ABY4": "flavodoxin", "P42212": "GFP", "P0CG48": "ubiquitin",
    "P0AEX9_dup": "MBP", "Q12306": "SUMO/Smt3", "P63165": "SUMO1",
    "P06654": "protein-G B1", "P00698": "hen-lysozyme", "P69441": "adenylate-kinase",
    "P0A6F5": "GroEL", "P00644": "staph-nuclease", "P0A9X4": "MreB",
}

QUERY = """
query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    exptl { method }
    rcsb_entry_info { resolution_combined }
    struct { title }
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity { pdbx_description }
      entity_poly { rcsb_sample_sequence_length }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
      rcsb_polymer_entity_align {
        reference_database_accession reference_database_name provenance_source
        aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
      }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name formula_weight } }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}
"""


def post(payload, timeout=240, tries=4):
    body = json.dumps(payload).encode()
    for a in range(tries):
        try:
            req = urllib.request.Request(GQL, data=body,
                                         headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=timeout))
        except Exception as exc:
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))
            sys.stderr.write(f"   retry {a+1} ({type(exc).__name__})\n")


def ranges_for(align, acc):
    out = []
    for b in align or []:
        if b.get("reference_database_name") != "UniProt":
            continue
        if b.get("reference_database_accession") != acc:
            continue
        for r in b.get("aligned_regions") or []:
            beg, ln = r.get("entity_beg_seq_id"), r.get("length")
            if beg and ln:
                out.append([int(beg), int(beg) + int(ln) - 1])
    return sorted(out)


def topology(tgt, fus):
    if not tgt or not fus:
        return None
    fmin, fmax = min(p for p, _ in fus), max(q for _, q in fus)
    if len(tgt) >= 2 and any(q < fmin for _, q in tgt) and any(p > fmax for p, _ in tgt):
        return "INTERNAL"
    if fmax < min(p for p, _ in tgt):
        return "TERMINAL_N"
    if fmin > max(q for _, q in tgt):
        return "TERMINAL_C"
    return "OTHER"


def analyse_entry(e):
    rec = {"pdb_id": e["rcsb_id"],
           "method": (e.get("exptl") or [{}])[0].get("method"),
           "resolution": ((e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None])[0],
           "title": (e.get("struct") or {}).get("title"),
           "entities": [], "nonpolymer": []}
    for ne in e.get("nonpolymer_entities") or []:
        cc = ((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
        if cc.get("id"):
            rec["nonpolymer"].append({
                "id": cc["id"], "name": cc.get("name"), "mw": cc.get("formula_weight"),
                "chains": (ne.get("rcsb_nonpolymer_entity_container_identifiers")
                           or {}).get("auth_asym_ids")})
    any_chimera = False
    for pe in e.get("polymer_entities") or []:
        align = pe.get("rcsb_polymer_entity_align") or []
        accs = sorted({b["reference_database_accession"] for b in align
                       if b.get("reference_database_name") == "UniProt"})
        ent = {"entity_id": pe["rcsb_id"],
               "description": (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
               "length": (pe.get("entity_poly") or {}).get("rcsb_sample_sequence_length"),
               "chains": (pe.get("rcsb_polymer_entity_container_identifiers")
                          or {}).get("auth_asym_ids") or [],
               "uniprot": accs,
               "provenance": sorted({b.get("provenance_source") for b in align
                                     if b.get("reference_database_name") == "UniProt"})}
        known = [a for a in accs if a in PARTNERS]
        if len(accs) >= 2:
            any_chimera = True
            ent["is_chimera"] = True
            if known:
                fa = known[0]
                ta = [a for a in accs if a != fa]
                ent["fusion_accession"] = fa
                ent["fusion_partner"] = PARTNERS[fa]
                ent["target_accession"] = ta[0] if len(ta) == 1 else None
                ent["target_accessions"] = ta
                ent["fusion_ranges"] = ranges_for(align, fa)
                ent["target_ranges"] = (ranges_for(align, ta[0]) if len(ta) == 1 else None)
                ent["topology"] = (topology(ent["target_ranges"], ent["fusion_ranges"])
                                   if ent["target_ranges"] else None)
                ent["construct_class"] = ("RECOGNISED_PARTNER_CHIMERA" if len(accs) == 2
                                          else "RECOGNISED_PARTNER_MULTI_CHIMERA")
            else:
                ent["construct_class"] = "OTHER_CHIMERA"
                ent["fusion_partner"] = None
        else:
            ent["is_chimera"] = False
            if len(accs) == 1 and accs[0] in PARTNERS:
                ent["construct_class"] = "STANDALONE_PARTNER_PROTEIN"
                ent["standalone_partner"] = PARTNERS[accs[0]]
            elif not accs:
                ent["construct_class"] = "NO_UNIPROT_MAPPING"
            else:
                ent["construct_class"] = "SINGLE_UNIPROT"
        rec["entities"].append(ent)
    rec["entry_has_chimera"] = any_chimera
    rec["entry_has_recognised_partner_chimera"] = any(
        x.get("construct_class", "").startswith("RECOGNISED_PARTNER") for x in rec["entities"])
    rec["entry_has_standalone_partner"] = any(
        x.get("construct_class") == "STANDALONE_PARTNER_PROTEIN" for x in rec["entities"])
    return rec


def main():
    reg = json.load(open(os.path.join(AUD, "dataset_registry.json")))
    ids = reg["all_unique_pdb_entries"]
    cache = os.path.join(AUD, "entry_constructs.json")
    done = {}
    if os.path.exists(cache):
        done = {r["pdb_id"]: r for r in json.load(open(cache))["entries"]}
        sys.stderr.write(f"resuming: {len(done)} entries cached\n")
    todo = [i for i in ids if i not in done]
    sys.stderr.write(f"entries to query: {len(todo)} of {len(ids)}\n")
    for i in range(0, len(todo), CHUNK):
        chunk = todo[i:i + CHUNK]
        try:
            d = post({"query": QUERY, "variables": {"ids": chunk}})
            for e in d["data"]["entries"] or []:
                if e:
                    done[e["rcsb_id"]] = analyse_entry(e)
        except Exception as exc:
            sys.stderr.write(f"  CHUNK FAILED {chunk[0]}..: {type(exc).__name__}: {exc}\n")
        if (i // CHUNK) % 10 == 0:
            sys.stderr.write(f"  {i + len(chunk)}/{len(todo)}  (have {len(done)})\n")
            json.dump({"entries": list(done.values())}, open(cache, "w"))
    json.dump({"n": len(done), "entries": list(done.values())}, open(cache, "w"))
    missing = [i for i in ids if i not in done]
    print(f"resolved {len(done)}/{len(ids)} entries; unresolved {len(missing)}")
    if missing:
        print(f"  unresolved (likely obsolete/withdrawn): {missing[:20]}")
    ch = sum(1 for r in done.values() if r["entry_has_recognised_partner_chimera"])
    oc = sum(1 for r in done.values()
             if r["entry_has_chimera"] and not r["entry_has_recognised_partner_chimera"])
    sa = sum(1 for r in done.values() if r["entry_has_standalone_partner"])
    print(f"  entries with a recognised-partner chimera: {ch}")
    print(f"  entries with some other chimera:           {oc}")
    print(f"  entries containing a STANDALONE partner protein (not a construct): {sa}")


if __name__ == "__main__":
    main()
