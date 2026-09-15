#!/usr/bin/env python
"""
Phase 5 STEP 3 — where does the dataset's ground-truth label sit on a chimeric member?

For every chimeric dataset member with a recognised crystallization partner, locate the
benchmark's labelled ligand(s) and classify the site they define:

  TARGET_SITE                 ligand contacts predominantly TARGET residues
  FUSION_PARTNER_SITE         ligand contacts predominantly FUSION residues   -> CASE B
  TARGET_FUSION_INTERFACE     substantial contact with both
  OTHER_CHAIN / AMBIGUOUS     ligand sits on a different chain, or cannot be resolved

Ground-truth ligand codes come from the frozen `(mlig)` membership files where the dataset has
one (MOAD-2013-derived relevant ligands). Datasets without an (mlig) file have no explicit
ligand label and are reported as LABEL_NOT_SPECIFIED.

  CASE A  a real fusion-partner cavity exists but is NOT labelled -> unlabelled construct decoy
  CASE B  the labelled ligand itself sits on the fusion partner   -> fusion pocket is ground truth

Outputs: dataset_audit/groundtruth_chimeric.json / .tsv
"""
import collections, json, os, sys, time, urllib.request
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "dataset_audit")
CACHE = os.path.join(ROOT, ".structure_cache")
CONTACT = 4.0
MIN_CONTACTS = 3
DOMINANCE = 0.70          # same convention as the main study


def fetch_cif(pdb_id):
    p = os.path.join(CACHE, f"{pdb_id.lower()}.cif")
    if not os.path.exists(p):
        for a in range(4):
            try:
                with urllib.request.urlopen(
                        f"https://files.rcsb.org/download/{pdb_id.lower()}.cif",
                        timeout=180) as r, open(p, "wb") as fh:
                    fh.write(r.read())
                break
            except Exception:
                if a == 3:
                    raise
                time.sleep(2 * (a + 1))
    return p


def expand(rs):
    s = set()
    for a, b in rs or []:
        s.update(range(a, b + 1))
    return s


def classify_ligand_site(pdb_id, ent, lig_codes):
    """Returns per-ligand-instance site classification for the chimeric entity."""
    cif = fetch_cif(pdb_id)
    st = gemmi.read_structure(cif)
    st.setup_entities()
    model = st[0]
    tgt = expand(ent.get("target_ranges"))
    fus = expand(ent.get("fusion_ranges"))
    chains = set(ent.get("chains") or [])
    out = []
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            if res.name in ("HOH", "DOD"):
                continue
            if lig_codes and res.name not in lig_codes:
                continue
            nt = nf = nother = 0
            for atom in res:
                for m in ns.find_atoms(atom.pos, '\0', radius=CONTACT):
                    cra = m.to_cra(model)
                    ri = gemmi.find_tabulated_residue(cra.residue.name)
                    if ri is None or not ri.is_amino_acid():
                        continue
                    if atom.pos.dist(cra.atom.pos) > CONTACT:
                        continue
                    if cra.chain.name not in chains:
                        nother += 1
                        continue
                    ls = cra.residue.label_seq
                    if ls is None:
                        continue
                    ls = int(ls)
                    if ls in tgt:
                        nt += 1
                    elif ls in fus:
                        nf += 1
                    else:
                        nother += 1
            tot = nt + nf
            if tot < MIN_CONTACTS:
                site = "OTHER_CHAIN_OR_AMBIGUOUS"
            elif nf / tot >= DOMINANCE:
                site = "FUSION_PARTNER_SITE"
            elif nt / tot >= DOMINANCE:
                site = "TARGET_SITE"
            elif nt and nf:
                site = "TARGET_FUSION_INTERFACE"
            else:
                site = "OTHER_CHAIN_OR_AMBIGUOUS"
            out.append({"ligand": res.name, "chain": ch.name, "seq": res.seqid.num,
                        "n_target_contacts": nt, "n_fusion_contacts": nf,
                        "n_other_contacts": nother, "site_class": site})
    return out


def main():
    reg = json.load(open(os.path.join(AUD, "dataset_registry.json")))
    ents = {r["pdb_id"]: r for r in
            json.load(open(os.path.join(AUD, "entry_constructs.json")))["entries"]}

    # ground-truth ligand codes, from the (mlig) membership files
    mlig = {}
    for name, d in reg["datasets"].items():
        if not name.endswith("_mlig"):
            continue
        base = name[:-5]
        for m in d["members"]:
            mlig.setdefault(base, {}).setdefault(m["pdb_id"], set()).update(m["ligands"])

    rows = []
    for name, d in reg["datasets"].items():
        if name.endswith("_mlig"):
            continue
        members = d["members"]
        for m in members:
            rec = ents.get(m["pdb_id"])
            if not rec:
                continue
            chim = [e for e in rec["entities"]
                    if str(e.get("construct_class", "")).startswith("RECOGNISED_PARTNER")]
            if not chim:
                continue
            for ent in chim:
                # is the benchmarked chain the chimeric one?
                benched = None
                if m["chain"]:
                    benched = m["chain"] in (ent.get("chains") or [])
                codes = mlig.get(name, {}).get(m["pdb_id"], set())
                label_src = "MLIG" if codes else "LABEL_NOT_SPECIFIED"
                rows.append({
                    "dataset": name, "role": d["role"], "pdb_id": m["pdb_id"],
                    "benchmark_chain": m["chain"], "subset": m["subset"],
                    "entity_id": ent["entity_id"],
                    "entity_chains": ",".join(ent.get("chains") or []),
                    "benchmarked_chain_is_chimeric": benched,
                    "fusion_partner": ent.get("fusion_partner"),
                    "fusion_accession": ent.get("fusion_accession"),
                    "target_accession": ent.get("target_accession"),
                    "topology": ent.get("topology"),
                    "label_source": label_src,
                    "labelled_ligands": ",".join(sorted(codes)) if codes else "",
                    "_ent": ent, "_codes": codes})

    print(f"chimeric dataset-member records to resolve: {len(rows)}")
    by_entry = collections.defaultdict(list)
    for r in rows:
        by_entry[(r["pdb_id"], r["entity_id"], frozenset(r["_codes"]))].append(r)

    cache = {}
    for i, key in enumerate(sorted(by_entry, key=lambda k: k[0]), 1):
        pdb, eid, codes = key
        ent = by_entry[key][0]["_ent"]
        try:
            cache[key] = classify_ligand_site(pdb, ent, set(codes))
        except Exception as exc:
            cache[key] = [{"error": f"{type(exc).__name__}: {exc}"}]
            sys.stderr.write(f"  {pdb} failed: {exc}\n")
        if i % 20 == 0:
            sys.stderr.write(f"  resolved {i}/{len(by_entry)}\n")

    for r in rows:
        key = (r["pdb_id"], r["entity_id"], frozenset(r["_codes"]))
        ligs = cache.get(key, [])
        r.pop("_ent"); r.pop("_codes")
        sites = [l.get("site_class") for l in ligs if "site_class" in l]
        r["n_labelled_ligand_instances"] = len(sites)
        r["site_classes"] = ";".join(sorted(set(sites))) if sites else ""
        r["n_on_target"] = sites.count("TARGET_SITE")
        r["n_on_fusion"] = sites.count("FUSION_PARTNER_SITE")
        r["n_on_interface"] = sites.count("TARGET_FUSION_INTERFACE")
        r["n_ambiguous"] = sites.count("OTHER_CHAIN_OR_AMBIGUOUS")
        r["CASE_B_fusion_is_ground_truth"] = r["n_on_fusion"] > 0
        r["CASE_A_unlabelled_fusion_decoy"] = (r["n_on_fusion"] == 0
                                               and r["label_source"] == "MLIG")
        r["ligand_detail"] = ";".join(
            f"{l['ligand']}/{l['chain']}{l['seq']}:{l['site_class']}"
            f"(T{l['n_target_contacts']}/F{l['n_fusion_contacts']})"
            for l in ligs if "site_class" in l)[:400]

    json.dump({"contact_A": CONTACT, "dominance": DOMINANCE, "rows": rows},
              open(os.path.join(AUD, "groundtruth_chimeric.json"), "w"), indent=1)
    cols = [c for c in rows[0].keys()] if rows else []
    with open(os.path.join(AUD, "groundtruth_chimeric.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    print("\n" + "=" * 96)
    print("GROUND-TRUTH LOCATION ON CHIMERIC DATASET MEMBERS")
    print("=" * 96)
    for name in reg["datasets"]:
        if name.endswith("_mlig"):
            continue
        sub = [r for r in rows if r["dataset"] == name]
        if not sub:
            continue
        b = sum(1 for r in sub if r["CASE_B_fusion_is_ground_truth"])
        a = sum(1 for r in sub if r["CASE_A_unlabelled_fusion_decoy"])
        ns = sum(1 for r in sub if r["label_source"] == "LABEL_NOT_SPECIFIED")
        print(f"\n{name} ({sub[0]['role']})")
        print(f"  chimeric member-entities: {len(sub)}")
        print(f"  partners: {dict(collections.Counter(r['fusion_partner'] for r in sub))}")
        print(f"  CASE B (labelled ligand sits on the FUSION partner): {b}")
        print(f"  CASE A (fusion cavity present, not labelled):        {a}")
        print(f"  no explicit ligand label in this dataset:            {ns}")
        if b:
            for r in sub:
                if r["CASE_B_fusion_is_ground_truth"]:
                    print(f"    CASE B: {r['pdb_id']} {r['fusion_partner']} "
                          f"{r['labelled_ligands']} -> {r['ligand_detail'][:110]}")


if __name__ == "__main__":
    main()
