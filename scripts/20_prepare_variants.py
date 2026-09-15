#!/usr/bin/env python
"""
Pre-specified STRUCTURE-REPRESENTATION and HETATM sensitivity variants (protocol v1.3 III.6).

  assembly1 : the designated chimeric chain PLUS every other polymer chain of biological
              assembly 1, so the detector sees the biological context rather than one chain.
              Residue classes are defined only for the designated chain; all other chains are
              OTHER_CHAIN and are excluded from f_target / f_fusion / f_linker (they are not
              TARGET, FUSION or LINKER), exactly as UNMAPPED residues already are.
  ions      : the primary single chain, but ion / simple-salt HETATM groups within 5 A of the
              retained polymer are kept instead of stripped.

Both variants write a manifest in the frozen format so that 15_run_confirmatory.py picks them up
and runs the frozen detector + classification path unchanged.

Chain-qualified keys: for assembly1 the auth_class_map is keyed "<chain>:<resnum><icode>" and a
chain-aware classifier is used. Equivalence with the frozen single-chain classifier is asserted by
tests/test_variant_classifier_equivalence.py.
"""
import hashlib, importlib.util, json, os, sys
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
CACHE = os.path.join(ROOT, ".structure_cache")
VARIANT = os.environ.get("FUSIONTAG_VARIANT", "assembly1")
PREPDIR = os.path.join(ROOT, "prepared_confirmatory", VARIANT)
os.makedirs(PREPDIR, exist_ok=True)

IONS = {"ZN", "MG", "CA", "NA", "K", "MN", "FE", "FE2", "CU", "CU1", "NI", "CD", "CO", "HG",
        "CL", "BR", "IOD", "F", "SO4", "PO4", "NO3", "CS", "RB", "SR", "BA", "LI", "PB", "AU",
        "AG", "PT", "MO", "W", "V", "CR"}
ION_DIST = 5.0


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


prep = _load("prep", "03_prepare_structures.py")


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build(rec, base_manrec):
    """Returns a manifest record in the frozen format, for this variant."""
    pdb_id = rec["pdb_id"]
    chain_id = rec["representative_chain"]
    cif = os.path.join(CACHE, f"{pdb_id.lower()}.cif")
    st = gemmi.read_structure(cif)
    st.setup_entities()
    pos_class = prep.classify_positions(rec)

    note = None
    if VARIANT == "assembly1":
        if st.assemblies:
            asm_name = st.assemblies[0].name
            try:
                st.transform_to_assembly(asm_name, gemmi.HowToNameCopiedChain.AddNumber)
                st.setup_entities()
                note = f"biological assembly {asm_name}"
            except Exception as exc:
                note = f"assembly generation failed ({type(exc).__name__}); used asymmetric unit"
        else:
            note = "no assembly record; used asymmetric unit"

    model = st[0]
    # locate every chain whose name starts with the designated auth chain id (copies get suffixes)
    designated = [ch for ch in model if ch.name == chain_id]
    if not designated:
        import re as _re
        designated = [ch for ch in model
                      if _re.fullmatch(_re.escape(chain_id) + r"\d*", ch.name)]
    if not designated:
        raise RuntimeError(f"{pdb_id}: designated chain {chain_id} absent after assembly build")
    prime = designated[0]

    out_letters = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    chain_map = {}

    def letter(name):
        if name not in chain_map:
            try:
                chain_map[name] = next(out_letters)
            except StopIteration:
                chain_map[name] = "z"
        return chain_map[name]

    chain_map[prime.name] = "A"
    next(out_letters)  # consume 'A'

    keep_chains = [prime] if VARIANT != "assembly1" else \
        [prime] + [ch for ch in model if ch is not prime]

    auth_map, residues = {}, []
    for ch in keep_chains:
        lt = letter(ch.name)
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if ch is prime:
                ls = res.label_seq
                cl = pos_class.get(int(ls), "UNMAPPED") if ls is not None else "UNMAPPED"
            else:
                cl = "OTHER_CHAIN"
            residues.append({"chain": lt, "res": res, "cls": cl, "key": key})
            auth_map[f"{lt}:{key}"] = cl
            if ch is prime:
                auth_map.setdefault(key, cl)      # unqualified fallback for the primary chain

    def write(keep, path):
        out = gemmi.Structure()
        out.spacegroup_hm = "P 1"
        out.cell = gemmi.UnitCell()
        m = gemmi.Model("1")
        chains = {}
        n = 0
        for r in residues:
            if r["cls"] not in keep:
                continue
            c = chains.setdefault(r["chain"], gemmi.Chain(r["chain"]))
            nr = gemmi.Residue()
            nr.name = r["res"].name
            nr.seqid = r["res"].seqid
            nr.het_flag = "A"
            for a in r["res"]:
                if a.altloc not in ("\0", "", "A"):
                    continue
                na = gemmi.Atom()
                na.name, na.element, na.pos = a.name, a.element, a.pos
                na.occ, na.b_iso, na.altloc = a.occ, a.b_iso, "\0"
                nr.add_atom(na)
            if len(nr):
                c.add_residue(nr)
                n += 1
        # ions variant: retain ion groups close to retained polymer
        if VARIANT == "ions":
            pts = [a.pos for r in residues if r["cls"] in keep for a in r["res"]]
            ic = gemmi.Chain("I")
            for ch in model:
                for res in ch:
                    if res.name not in IONS:
                        continue
                    if any(a.pos.dist(p) <= ION_DIST for a in res for p in pts):
                        nr = gemmi.Residue()
                        nr.name = res.name
                        nr.seqid = res.seqid
                        nr.het_flag = "H"
                        for a in res:
                            na = gemmi.Atom()
                            na.name, na.element, na.pos = a.name, a.element, a.pos
                            na.occ, na.b_iso = a.occ, a.b_iso
                            nr.add_atom(na)
                        ic.add_residue(nr)
            if len(ic):
                chains["I"] = ic
        for c in chains.values():
            m.add_chain(c)
        out.add_model(m)
        out.setup_entities()
        out.write_pdb(path)
        return n

    base = f"{pdb_id}_{chain_id}"
    op = os.path.join(PREPDIR, f"{base}_ORIGINAL.pdb")
    rp = os.path.join(PREPDIR, f"{base}_REMOVED.pdb")
    keep_o = {"TARGET", "FUSION", "LINKER"} | ({"OTHER_CHAIN"} if VARIANT == "assembly1" else set())
    keep_r = {"TARGET"} | ({"OTHER_CHAIN"} if VARIANT == "assembly1" else set())
    no = write(keep_o, op)
    nr_ = write(keep_r, rp)

    out = dict(base_manrec)
    out.update({
        "variant": VARIANT,
        "variant_note": note,
        "n_chains_kept": len(chain_map),
        "chain_rename_map": chain_map,
        "original_pdb": os.path.relpath(op, ROOT),
        "removed_pdb": os.path.relpath(rp, ROOT),
        "original_sha256": sha256_file(op), "removed_sha256": sha256_file(rp),
        "n_residues_original": no, "n_residues_removed": nr_,
        "auth_class_map": auth_map,
    })
    return out


def main():
    structures = json.load(open(os.path.join(MANI,
                                             "confirmatory_set_final.json")))["structures"]
    basem = {p["pdb_id"]: p for p in
             json.load(open(os.path.join(MANI,
                                         "confirmatory_manifest_primary.json")))["prepared"]
             if "PREPARATION_FAILED" not in p}
    out = []
    for i, rec in enumerate(structures, 1):
        if rec["pdb_id"] not in basem:
            continue
        try:
            out.append(build(rec, basem[rec["pdb_id"]]))
        except Exception as exc:
            out.append({"pdb_id": rec["pdb_id"], "entity_id": rec["entity_id"],
                        "fusion_partner": rec["fusion_partner"],
                        "target_accession": rec["target_accession"],
                        "PREPARATION_FAILED": f"{type(exc).__name__}: {exc}"})
            sys.stderr.write(f"  FAIL {rec['pdb_id']}: {type(exc).__name__}: {exc}\n")
        if i % 25 == 0:
            sys.stderr.write(f"  {VARIANT} prepared {i}/{len(structures)}\n")
    json.dump({"protocol_version": "1.3", "variant": VARIANT, "prepared": out},
              open(os.path.join(MANI, f"confirmatory_manifest_{VARIANT}.json"), "w"), indent=1)
    ok = [o for o in out if "PREPARATION_FAILED" not in o]
    print(f"{VARIANT}: prepared {len(ok)}/{len(out)}")
    if VARIANT == "assembly1":
        import collections
        c = collections.Counter(o["n_chains_kept"] for o in ok)
        print(f"  chains per structure: {dict(sorted(c.items()))}")
        notes = collections.Counter(o["variant_note"] for o in ok if o["variant_note"])
        print(f"  notes: {dict(notes)}")
    if VARIANT == "ions":
        import re
        n_with = 0
        for o in ok:
            txt = open(os.path.join(ROOT, o["original_pdb"]), encoding="utf-8").read()
            if "HETATM" in txt:
                n_with += 1
        print(f"  structures with retained ion HETATM: {n_with}/{len(ok)}")


if __name__ == "__main__":
    main()
