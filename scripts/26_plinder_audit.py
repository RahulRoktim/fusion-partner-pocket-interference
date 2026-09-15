#!/usr/bin/env python
"""
Phase 5 STEP 8 — PLINDER construct audit, release 2024-06 / v2.

PLINDER is far too large to SIFTS-scan exhaustively (107,963 PDB entries). Instead the audit runs
in reverse: enumerate every PDB entry in which a RECOGNISED crystallization fusion partner appears
as part of a chimeric polymer entity (RCSB search by UniProt accession + SIFTS chimera test), then
intersect that enumeration with PLINDER's frozen split membership. That is exact for the partners
enumerated and makes no claim about partners outside the list.

Two distinctions are preserved, as required:
  * PDB ENTRY contains a fusion construct   vs
  * the PLINDER SYSTEM's RECEPTOR chain is itself the chimeric chain
and for each system whose receptor is chimeric, whether the system's ligand sits on
TARGET / FUSION / INTERFACE.

Outputs: dataset_audit/plinder_audit.json / .tsv
"""
import collections, json, os, sys, time, urllib.request
import pandas as pd
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "dataset_audit")
MEM = os.path.join(AUD, "membership")
CACHE = os.path.join(ROOT, ".structure_cache")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GQL = "https://data.rcsb.org/graphql"
CONTACT = 4.0
DOMINANCE = 0.70

PARTNERS = {"P0ABE7": "BRIL/cyt-b562", "P00720": "T4L", "P0AEX9": "MBP",
            "P00268": "rubredoxin", "P0AA25": "thioredoxin", "P08515": "GST-Sj26",
            "P0ABY4": "flavodoxin", "P42212": "GFP", "P0CG48": "ubiquitin",
            "Q12306": "SUMO/Smt3", "P63165": "SUMO1", "P06654": "protein-G B1"}

QUERY = """
query($ids:[String!]!){
  polymer_entities(entity_ids:$ids){
    rcsb_id
    rcsb_polymer_entity { pdbx_description }
    rcsb_polymer_entity_container_identifiers { auth_asym_ids }
    rcsb_polymer_entity_align {
      reference_database_accession reference_database_name
      aligned_regions { entity_beg_seq_id length }
    }
  }
}"""


def post(url, payload, timeout=240, tries=4):
    body = json.dumps(payload).encode()
    for a in range(tries):
        try:
            r = urllib.request.Request(url, data=body,
                                       headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(r, timeout=timeout))
        except Exception:
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))


def search_entities(acc):
    d = post(SEARCH, {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}},
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": 10000},
                            "results_verbosity": "compact"}})
    return d.get("result_set", [])


def ranges_for(align, acc):
    out = []
    for b in align or []:
        if b.get("reference_database_name") != "UniProt" or \
                b.get("reference_database_accession") != acc:
            continue
        for r in b.get("aligned_regions") or []:
            if r.get("entity_beg_seq_id") and r.get("length"):
                out.append([int(r["entity_beg_seq_id"]),
                            int(r["entity_beg_seq_id"]) + int(r["length"]) - 1])
    return sorted(out)


def build_chimera_index():
    cache = os.path.join(AUD, "known_chimera_index.json")
    if os.path.exists(cache):
        sys.stderr.write("using cached chimera index\n")
        return json.load(open(cache))
    idx = {}
    for acc, tag in PARTNERS.items():
        ids = search_entities(acc)
        sys.stderr.write(f"  {tag} ({acc}): {len(ids)} entities with this accession\n")
        recs = []
        for i in range(0, len(ids), 50):
            d = post(GQL, {"query": QUERY, "variables": {"ids": ids[i:i + 50]}})
            recs += [r for r in d["data"]["polymer_entities"] if r]
        n = 0
        for r in recs:
            align = r.get("rcsb_polymer_entity_align") or []
            accs = sorted({b["reference_database_accession"] for b in align
                           if b.get("reference_database_name") == "UniProt"})
            if len(accs) < 2 or acc not in accs:
                continue
            tgts = [a for a in accs if a != acc]
            entry = r["rcsb_id"].split("_")[0].upper()
            idx.setdefault(entry, []).append({
                "entity_id": r["rcsb_id"], "fusion_accession": acc, "fusion_partner": tag,
                "target_accessions": tgts,
                "target_accession": tgts[0] if len(tgts) == 1 else None,
                "chains": (r.get("rcsb_polymer_entity_container_identifiers")
                           or {}).get("auth_asym_ids") or [],
                "description": (r.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
                "fusion_ranges": ranges_for(align, acc),
                "target_ranges": ranges_for(align, tgts[0]) if len(tgts) == 1 else None})
            n += 1
        sys.stderr.write(f"      -> {n} chimeric entities\n")
    json.dump(idx, open(cache, "w"))
    return idx


def expand(rs):
    s = set()
    for a, b in rs or []:
        s.update(range(a, b + 1))
    return s


def fetch_cif(p):
    path = os.path.join(CACHE, f"{p.lower()}.cif")
    if not os.path.exists(path):
        for a in range(3):
            try:
                with urllib.request.urlopen(
                        f"https://files.rcsb.org/download/{p.lower()}.cif",
                        timeout=180) as r, open(path, "wb") as fh:
                    fh.write(r.read())
                break
            except Exception:
                if a == 2:
                    raise
                time.sleep(2)
    return path


def ligand_site(pdb, ent, lig_chain_labels):
    """Classify where the system's ligand sits, using label_asym ids from the system id."""
    st = gemmi.read_structure(fetch_cif(pdb))
    st.setup_entities()
    model = st[0]
    tgt, fus = expand(ent.get("target_ranges")), expand(ent.get("fusion_ranges"))
    chains = {c.upper() for c in ent.get("chains") or []}
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    want = {x.split(".")[-1].upper() for x in lig_chain_labels}
    res_out = []
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            if res.name in ("HOH", "DOD"):
                continue
            if want and ch.name.upper() not in want and res.subchain.upper() not in want:
                continue
            nt = nf = 0
            for atom in res:
                for m in ns.find_atoms(atom.pos, '\0', radius=CONTACT):
                    cra = m.to_cra(model)
                    ri = gemmi.find_tabulated_residue(cra.residue.name)
                    if ri is None or not ri.is_amino_acid():
                        continue
                    if atom.pos.dist(cra.atom.pos) > CONTACT:
                        continue
                    if cra.chain.name.upper() not in chains:
                        continue
                    ls = cra.residue.label_seq
                    if ls is None:
                        continue
                    if int(ls) in tgt:
                        nt += 1
                    elif int(ls) in fus:
                        nf += 1
            tot = nt + nf
            if tot < 3:
                site = "NOT_ON_THIS_ENTITY"
            elif nf / tot >= DOMINANCE:
                site = "FUSION"
            elif nt / tot >= DOMINANCE:
                site = "TARGET"
            else:
                site = "INTERFACE"
            res_out.append({"ligand": res.name, "chain": ch.name, "site": site,
                            "n_target": nt, "n_fusion": nf})
    return res_out


def main():
    print("building the recognised-partner chimera index from RCSB/SIFTS ...")
    idx = build_chimera_index()
    chim_entries = set(idx)
    print(f"known chimeric PDB entries (recognised partners): {len(chim_entries)}")

    df = pd.read_parquet(os.path.join(MEM, "plinder_2024-06_v2_split.parquet"))
    df["pdb"] = df["system_id"].str.split("__").str[0].str.upper()
    print(f"PLINDER 2024-06/v2: {len(df)} systems, {df['pdb'].nunique()} entries")
    print(f"  split sizes: {df['split'].value_counts().to_dict()}")

    hit = df[df["pdb"].isin(chim_entries)].copy()
    print(f"\nPLINDER systems whose PDB ENTRY contains a recognised fusion construct: {len(hit)}")
    print(f"  distinct entries: {hit['pdb'].nunique()}")
    print(f"  by split (systems): {hit['split'].value_counts().to_dict()}")
    print(f"  by split (entries): "
          f"{hit.groupby('split')['pdb'].nunique().to_dict()}")

    rows = []
    for _, r in hit.iterrows():
        parts = r["system_id"].split("__")
        rec_chains = parts[2].split("_") if len(parts) > 2 else []
        lig_chains = parts[3].split("_") if len(parts) > 3 else []
        ents = idx[r["pdb"]]
        rec_labels = {c.split(".")[-1].upper() for c in rec_chains}
        receptor_is_chimeric = any(
            rec_labels & {c.upper() for c in e["chains"]} for e in ents)
        rows.append({"system_id": r["system_id"], "pdb_id": r["pdb"], "split": r["split"],
                     "receptor_chains": ",".join(rec_chains),
                     "ligand_chains": ",".join(lig_chains),
                     "fusion_partners": ",".join(sorted({e["fusion_partner"] for e in ents})),
                     "receptor_is_chimeric_chain": receptor_is_chimeric,
                     "passes_validation": bool(r["system_pass_validation_criteria"])})
    print(f"\n  of these, systems whose RECEPTOR chain is itself the chimeric chain: "
          f"{sum(1 for x in rows if x['receptor_is_chimeric_chain'])}")
    rec = [x for x in rows if x["receptor_is_chimeric_chain"]]
    print(f"    by split: {collections.Counter(x['split'] for x in rec)}")
    print(f"    distinct entries: {len({x['pdb_id'] for x in rec})}")
    print(f"    partners: {collections.Counter(x['fusion_partners'] for x in rec)}")

    # ligand location, limited to the non-'removed' splits (train/val/test) to bound cost
    keep = [x for x in rec if x["split"] in ("train", "val", "test")]
    print(f"\n  resolving ligand location for {len(keep)} systems in train/val/test ...")
    by_entry = collections.defaultdict(list)
    for x in keep:
        by_entry[x["pdb_id"]].append(x)
    for i, pdb in enumerate(sorted(by_entry), 1):
        for e in idx[pdb]:
            if e.get("target_ranges") is None:
                continue
            try:
                for x in by_entry[pdb]:
                    ligs = ligand_site(pdb, e, x["ligand_chains"].split(","))
                    on = [l for l in ligs if l["site"] != "NOT_ON_THIS_ENTITY"]
                    x["ligand_sites"] = ";".join(f"{l['ligand']}:{l['site']}" for l in on)
                    x["site_class"] = (";".join(sorted({l["site"] for l in on}))
                                       if on else "NOT_ON_CHIMERIC_ENTITY")
            except Exception as exc:
                for x in by_entry[pdb]:
                    x["site_class"] = f"ERROR:{type(exc).__name__}"
            break
        if i % 15 == 0:
            sys.stderr.write(f"    {i}/{len(by_entry)}\n")

    print("\n  ligand-site location for chimeric-receptor systems in train/val/test:")
    for sp in ("train", "val", "test"):
        sub = [x for x in keep if x["split"] == sp]
        c = collections.Counter(x.get("site_class", "unresolved") for x in sub)
        print(f"    {sp:6s} n={len(sub):4d}  {dict(c)}")

    json.dump({"release": "2024-06 / v2",
               "source": "https://storage.googleapis.com/plinder/2024-06/v2/splits/split.parquet",
               "n_systems": int(len(df)), "n_entries": int(df["pdb"].nunique()),
               "split_sizes": {k: int(v) for k, v in df["split"].value_counts().items()},
               "n_systems_entry_has_construct": len(rows),
               "n_systems_receptor_is_chimeric": len(rec),
               "rows": rows},
              open(os.path.join(AUD, "plinder_audit.json"), "w"), indent=1)
    cols = sorted({k for x in rows for k in x})
    with open(os.path.join(AUD, "plinder_audit.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
