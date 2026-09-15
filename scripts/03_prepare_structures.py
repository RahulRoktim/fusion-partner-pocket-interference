#!/usr/bin/env python
"""
Phase 2B/2C — fetch structures and build the paired ORIGINAL / FUSION-REMOVED inputs.

Protocol v1.1 section 4:
  * single designated chimeric chain, first model, altloc blank or A
  * all non-polymer stripped (ligands, waters, ions, lipids, detergents, cryoprotectants)
  * [A1] TAG residues deleted from BOTH conditions
  * ORIGINAL  = target + fusion + linker
  * REMOVED   = target only
  * no protonation, no minimisation, no renumbering; TARGET coordinates untouched
  * [A3] deletion boundaries recorded for the REMOVED condition

Reference sites (protocol 7.4, outcome-independent) are derived here from deposited ligand
coordinates plus a pre-declared artifact exclusion list. No detector output is consulted.

Structures are cached under .structure_cache/ and are NOT redistributed with this repository.
"""
import json, os, sys, hashlib, urllib.request, time
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
CACHE = os.path.join(ROOT, ".structure_cache")
PREP = os.path.join(ROOT, "prepared")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(PREP, exist_ok=True)

REF_SITE_DIST = 4.0          # protocol 7.4
MIN_LIGAND_MW = 100.0        # below this: ion/fragment, not a biological reference ligand
MIN_LIGAND_ATOMS = 6

# ---------------------------------------------------------------------------
# PRE-DECLARED artifact-ligand exclusion list (protocol 7.4).
# Fixed before any detector was run. Purpose: exclude crystallization and cryo additives,
# buffers, detergents and membrane mimetics from *target reference site* definition.
# ---------------------------------------------------------------------------
ARTIFACT_LIGANDS = {
    # solvent / cryo / precipitant
    "HOH", "DOD", "GOL", "EDO", "PEG", "PGE", "PG4", "1PE", "2PE", "P6G", "7PE", "12P", "15P",
    "MPD", "MRD", "DMS", "DMF", "IPA", "EOH", "MOH", "ACN", "TFA", "SCN",
    # buffers
    "TRS", "MES", "EPE", "BTB", "CIT", "FLC", "TLA", "MLA", "MLI", "ACT", "ACY", "FMT", "OXL",
    "IMD", "BME", "DTT", "DTU", "TCE", "CAC", "PIN", "HEZ", "BEZ",
    # common ions / salts
    "SO4", "PO4", "NO3", "CL", "BR", "IOD", "F", "NA", "K", "LI", "CS", "RB", "MG", "CA", "MN",
    "CD", "HG", "PB", "AU", "AG", "PT", "BA", "SR", "NH4", "AZI", "CO3", "PER", "UNX", "UNL",
    # detergents / lipids / membrane mimetics (GPCR constructs)
    "OLC", "OLA", "OLB", "MYS", "PLM", "STE", "PEE", "PCW", "PC1", "PGW", "LMT", "LMN", "LDA",
    "DDQ", "C8E", "C10", "CE9", "JEF", "BNG", "HEX", "D10", "D12", "UND", "TRD", "HP6",
    "CLR", "Y01", "CHD", "CPS", "MAU", "9PE", "PX4", "SOG", "NAG_ARTIFACT_PLACEHOLDER",
}
# Glycans are a separate category: chemically real but generally not the drug-target pocket.
GLYCANS = {"NAG", "NDG", "BMA", "MAN", "GAL", "FUC", "GLC", "XYS", "SIA", "BGC", "A2G"}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def fetch_cif(pdb_id):
    pdb_id = pdb_id.lower()
    path = os.path.join(CACHE, f"{pdb_id}.cif")
    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(url, timeout=180) as r, open(path, "wb") as out:
                    out.write(r.read())
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt + 1))
        os.chmod(path, 0o444)          # immutable raw data
    return path, sha256_file(path)


def classify_positions(rec):
    """label_seq_id -> TARGET / FUSION / LINKER / TAG, from the SIFTS segments in the manifest."""
    cls = {}
    for a, b in rec["target_ranges"]:
        for p in range(a, b + 1):
            cls[p] = "TARGET"
    for a, b in rec["fusion_ranges"]:
        for p in range(a, b + 1):
            cls[p] = "FUSION"
    for p in rec.get("linker_positions", []):
        cls[p] = "LINKER"
    for p in rec.get("tag_positions", []):
        cls[p] = "TAG"
    return cls


def prepare_one(rec):
    pdb_id = rec["pdb_id"]
    chain_id = rec["representative_chain"]
    cif_path, cif_hash = fetch_cif(pdb_id)

    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    model = st[0]                                     # first model only

    chain = None
    for ch in model:
        if ch.name == chain_id:
            chain = ch
            break
    if chain is None:
        raise RuntimeError(f"{pdb_id}: auth chain {chain_id} not found")

    pos_class = classify_positions(rec)

    # ---- walk the polymer of the chosen chain, mapping label_seq -> auth numbering
    residues = []          # (auth_seq, icode, label_seq, class, gemmi.Residue)
    unmapped = 0
    for res in chain:
        info = gemmi.find_tabulated_residue(res.name)
        is_aa = info is not None and info.is_amino_acid()
        if not is_aa:
            continue                                   # non-polymer stripped here
        lseq = res.label_seq
        if lseq is None:
            unmapped += 1
            continue
        residues.append({"auth_seq": res.seqid.num,
                         "icode": res.seqid.icode.strip(),
                         "label_seq": int(lseq),
                         "cls": pos_class.get(int(lseq), "UNMAPPED"),
                         "res": res,
                         "name": res.name})

    counts = {}
    for r in residues:
        counts[r["cls"]] = counts.get(r["cls"], 0) + 1

    # ---- reference site (outcome-independent): biological ligand contacting TARGET residues
    target_res_objs = [r["res"] for r in residues if r["cls"] == "TARGET"]
    ns = gemmi.NeighborSearch(st, 5.0).populate()
    ref_site, ref_ligands, rejected_ligands = set(), [], []
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            if res.name in ARTIFACT_LIGANDS or res.name in GLYCANS:
                rejected_ligands.append({"id": res.name, "reason": "artifact_or_glycan_list"})
                continue
            heavy = [a for a in res if a.element != gemmi.Element("H")]
            if len(heavy) < MIN_LIGAND_ATOMS:
                rejected_ligands.append({"id": res.name,
                                         "reason": f"too_few_heavy_atoms={len(heavy)}"})
                continue
            contacts = set()
            for atom in heavy:
                for m in ns.find_atoms(atom.pos, '\0', radius=REF_SITE_DIST):
                    cra = m.to_cra(model)
                    if cra.chain.name != chain_id:
                        continue
                    ls = cra.residue.label_seq
                    if ls is None:
                        continue
                    if pos_class.get(int(ls)) == "TARGET":
                        if atom.pos.dist(cra.atom.pos) <= REF_SITE_DIST:
                            contacts.add((cra.residue.seqid.num,
                                          cra.residue.seqid.icode.strip()))
            if len(contacts) >= 3:
                ref_ligands.append({"id": res.name, "chain": ch.name,
                                    "seq": res.seqid.num, "n_target_contacts": len(contacts)})
                ref_site |= contacts
            else:
                rejected_ligands.append({"id": res.name,
                                         "reason": f"only_{len(contacts)}_target_contacts"})

    # ---- write the two conditions
    def write_condition(keep_classes, out_path):
        out = gemmi.Structure()
        out.spacegroup_hm = "P 1"
        out.cell = gemmi.UnitCell()
        m = gemmi.Model("1")
        c = gemmi.Chain("A")                      # single-char chain for PDB-format safety
        n = 0
        for r in residues:
            if r["cls"] not in keep_classes:
                continue
            nr = gemmi.Residue()
            nr.name = r["name"]
            nr.seqid = gemmi.SeqId(r["auth_seq"], r["icode"] or " ")
            nr.het_flag = "A"
            for atom in r["res"]:
                if atom.altloc not in ("\0", "", "A"):
                    continue                       # altloc A (or none) only
                a = gemmi.Atom()
                a.name, a.element, a.pos = atom.name, atom.element, atom.pos
                a.occ, a.b_iso, a.altloc = atom.occ, atom.b_iso, "\0"
                nr.add_atom(a)
            if len(nr) == 0:
                continue
            c.add_residue(nr)
            n += 1
        m.add_chain(c)
        out.add_model(m)
        out.setup_entities()
        out.write_pdb(out_path)
        return n

    base = f"{pdb_id}_{chain_id}"
    orig_path = os.path.join(PREP, f"{base}_ORIGINAL.pdb")
    rem_path = os.path.join(PREP, f"{base}_REMOVED.pdb")
    n_orig = write_condition({"TARGET", "FUSION", "LINKER"}, orig_path)
    n_rem = write_condition({"TARGET"}, rem_path)

    # ---- [A3] deletion boundaries: target residues flanking each excised span
    ordered = sorted(residues, key=lambda r: (r["label_seq"],))
    kept_removed = [r for r in ordered if r["cls"] == "TARGET"]
    boundaries = []
    for i in range(len(kept_removed) - 1):
        a, b = kept_removed[i], kept_removed[i + 1]
        gap = [r for r in ordered
               if a["label_seq"] < r["label_seq"] < b["label_seq"]
               and r["cls"] in ("FUSION", "LINKER")]
        if not gap:
            continue
        for side in (a, b):
            ca = side["res"].find_atom("CA", "*")
            if ca is not None:
                boundaries.append({"auth_seq": side["auth_seq"], "label_seq": side["label_seq"],
                                   "x": ca.pos.x, "y": ca.pos.y, "z": ca.pos.z,
                                   "excised_span": [gap[0]["label_seq"], gap[-1]["label_seq"]],
                                   "n_excised": len(gap)})
    # terminal fusions: boundary is the single target residue adjacent to the excised terminus
    if kept_removed and not boundaries:
        fus_lseq = [r["label_seq"] for r in ordered if r["cls"] in ("FUSION", "LINKER")]
        if fus_lseq:
            if max(fus_lseq) < kept_removed[0]["label_seq"]:
                side = kept_removed[0]
            else:
                side = kept_removed[-1]
            ca = side["res"].find_atom("CA", "*")
            if ca is not None:
                boundaries.append({"auth_seq": side["auth_seq"], "label_seq": side["label_seq"],
                                   "x": ca.pos.x, "y": ca.pos.y, "z": ca.pos.z,
                                   "excised_span": [min(fus_lseq), max(fus_lseq)],
                                   "n_excised": len(fus_lseq),
                                   "terminal": True})

    # auth_seq -> class map, for mapping detector output back
    auth_map = {}
    for r in residues:
        key = f"{r['auth_seq']}{r['icode']}"
        auth_map[key] = r["cls"]

    return {
        "pdb_id": pdb_id, "entity_id": rec["entity_id"], "auth_chain": chain_id,
        "fusion_partner": rec["fusion_partner"], "target_accession": rec["target_accession"],
        "fusion_accession": rec["fusion_accession"], "topology": rec["topology"],
        "resolution": rec["resolution"], "method": rec["method"],
        "description": rec["description"], "title": rec["title"],
        "target_ranges_entity": rec["target_ranges"],
        "fusion_ranges_entity": rec["fusion_ranges"],
        "linker_positions_entity": rec.get("linker_positions", []),
        "tag_positions_entity": rec.get("tag_positions", []),
        "sifts_provenance": rec.get("provenance"),
        "cif_sha256": cif_hash,
        "residue_class_counts": counts,
        "n_unmapped_residues": unmapped,
        "tag_residues_deleted": [f"{r['auth_seq']}{r['icode']}" for r in residues
                                 if r["cls"] == "TAG"],
        "n_residues_original": n_orig, "n_residues_removed": n_rem,
        "original_pdb": os.path.relpath(orig_path, ROOT),
        "removed_pdb": os.path.relpath(rem_path, ROOT),
        "original_sha256": sha256_file(orig_path), "removed_sha256": sha256_file(rem_path),
        "deletion_boundaries": boundaries,
        "reference_site_residues": sorted(f"{n}{i}" for n, i in ref_site),
        "reference_site_ligands": ref_ligands,
        "reference_site_available": len(ref_site) > 0,
        "rejected_ligands": rejected_ligands,
        "auth_class_map": auth_map,
    }


def main():
    sel = json.load(open(os.path.join(MANI, "pilot_selection.json")))
    out = []
    for i, rec in enumerate(sel["selection"], 1):
        sys.stderr.write(f"[{i}/{len(sel['selection'])}] {rec['pdb_id']} "
                         f"{rec['fusion_partner']}\n")
        try:
            out.append(prepare_one(rec))
        except Exception as exc:
            out.append({"pdb_id": rec["pdb_id"], "entity_id": rec["entity_id"],
                        "fusion_partner": rec["fusion_partner"],
                        "PREPARATION_FAILED": f"{type(exc).__name__}: {exc}"})
            sys.stderr.write(f"    FAILED: {type(exc).__name__}: {exc}\n")

    json.dump({"protocol_version": "1.1", "prepared": out},
              open(os.path.join(MANI, "pilot_manifest.json"), "w"), indent=1)

    ok = [o for o in out if "PREPARATION_FAILED" not in o]
    print(f"\nprepared {len(ok)}/{len(out)} structures")
    print(f"{'pdb':6s} {'chain':5s} {'partner':6s} {'topo':18s} "
          f"{'tgt':>5s} {'fus':>5s} {'lnk':>4s} {'tag':>4s} {'bnd':>4s} refsite")
    for o in ok:
        c = o["residue_class_counts"]
        print(f"{o['pdb_id']:6s} {o['auth_chain']:5s} {o['fusion_partner']:6s} "
              f"{str(o['topology']):18s} {c.get('TARGET',0):5d} {c.get('FUSION',0):5d} "
              f"{c.get('LINKER',0):4d} {c.get('TAG',0):4d} {len(o['deletion_boundaries']):4d} "
              f"{'YES:' + ','.join(l['id'] for l in o['reference_site_ligands'][:3]) if o['reference_site_available'] else 'no'}")
    cov = sum(1 for o in ok if o["reference_site_available"])
    print(f"\nreference-site coverage: {cov}/{len(ok)} = {100*cov/max(len(ok),1):.0f}%")


if __name__ == "__main__":
    main()
