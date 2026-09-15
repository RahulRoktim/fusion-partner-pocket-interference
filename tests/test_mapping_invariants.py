#!/usr/bin/env python
"""
Phase 2C — parser / mapping validation. MUST pass before the 96-run pilot executes.

Validates the residue-class mapping and the paired-condition construction on three
structures chosen to cover the three construct geometries:

  5IU7  BRIL internal insertion into A2A adenosine receptor ICL3
  3ODU  T4 lysozyme internal insertion into CXCR4
  3L2J  MBP N-terminal fusion to the PTH1R extracellular domain

Invariants checked:
  I1  target residues map to the SIFTS target segments
  I2  fusion residues map to the SIFTS fusion segment
  I3  linker residues are the unaligned positions BETWEEN the two accessions
  I4  TAG residues are absent from BOTH written conditions
  I5  the REMOVED condition alters no TARGET atom coordinate
  I6  residue numbering stays traceable to the deposited auth numbering
  I7  written conditions contain exactly the intended residue classes
  I8  detector-style residue IDs map back onto the class map

Run:  python tests/test_mapping_invariants.py
"""
import json, os, sys, importlib.util

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(ROOT, "scripts", filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


build = _load("build_pool", "01_build_pool.py")
prep = _load("prepare", "03_prepare_structures.py")

import gemmi

CASES = [
    {"entity": "5IU7_1", "partner": "P0ABE7", "tag": "BRIL",
     "expect_topology": "INTERNAL_INSERTION"},
    {"entity": "3ODU_1", "partner": "P00720", "tag": "T4L",
     "expect_topology": "INTERNAL_INSERTION"},
    {"entity": "3L2J_1", "partner": "P0AEX9", "tag": "MBP",
     "expect_topology": "TERMINAL_N"},
]

FAILURES = []
CHECKS = [0]


def check(cond, label, detail=""):
    CHECKS[0] += 1
    if cond:
        print(f"    PASS  {label}")
    else:
        print(f"    FAIL  {label}  {detail}")
        FAILURES.append(f"{label} {detail}")


def read_pdb_residues(path):
    """auth-numbering residue keys and atom coordinates from a written condition."""
    st = gemmi.read_structure(path)
    out = {}
    for ch in st[0]:
        for res in ch:
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            out[key] = {a.name: (round(a.pos.x, 3), round(a.pos.y, 3), round(a.pos.z, 3))
                        for a in res}
    return out


def main():
    print("=" * 78)
    print("PHASE 2C — PARSER / MAPPING VALIDATION")
    print("=" * 78)

    for case in CASES:
        print(f"\n--- {case['entity']}  ({case['tag']}) ---")
        data = build.post(build.GQL, {"query": build.QUERY,
                                      "variables": {"ids": [case["entity"]]}})
        raw = data["data"]["polymer_entities"][0]
        raw["_partner_accession"] = case["partner"]
        raw["_partner_tag"] = case["tag"]
        rec = build.evaluate(raw)

        print(f"    accessions={rec['uniprot_accessions']} topology={rec['topology']} "
              f"chain={rec['representative_chain']} res={rec['resolution']}")
        print(f"    target_ranges={rec['target_ranges']}")
        print(f"    fusion_ranges={rec['fusion_ranges']}")
        print(f"    linker={rec['linker_positions']}  tag(n)={len(rec['tag_positions'])}")

        check(rec["topology"] == case["expect_topology"],
              f"topology == {case['expect_topology']}", f"got {rec['topology']}")

        out = prep.prepare_one(rec)
        cmap = out["auth_class_map"]
        counts = out["residue_class_counts"]
        print(f"    classes: {counts}")

        # ---- I1/I2: entity-position classification round-trips through label_seq
        pos_class = prep.classify_positions(rec)
        tgt_positions = {p for p, c in pos_class.items() if c == "TARGET"}
        fus_positions = {p for p, c in pos_class.items() if c == "FUSION"}
        expect_t = set()
        for a, b in rec["target_ranges"]:
            expect_t |= set(range(a, b + 1))
        expect_f = set()
        for a, b in rec["fusion_ranges"]:
            expect_f |= set(range(a, b + 1))
        check(tgt_positions == expect_t, "I1 target positions == SIFTS target segments")
        check(fus_positions == expect_f, "I2 fusion positions == SIFTS fusion segment")
        check(not (expect_t & expect_f), "I1/I2 target and fusion segments disjoint")

        # ---- I3: linker lies strictly between the outermost aligned positions
        aligned = expect_t | expect_f
        lo, hi = min(aligned), max(aligned)
        lnk = set(rec["linker_positions"])
        check(all(lo < p < hi for p in lnk), "I3 linker strictly interior to aligned span")
        check(not (lnk & aligned), "I3 linker disjoint from aligned segments")
        tagp = set(rec["tag_positions"])
        check(all(p < lo or p > hi for p in tagp), "I3 tag positions are terminal only")

        # ---- I4/I7: written conditions hold exactly the intended classes
        orig = read_pdb_residues(os.path.join(ROOT, out["original_pdb"]))
        rem = read_pdb_residues(os.path.join(ROOT, out["removed_pdb"]))
        orig_classes = {cmap.get(k, "UNMAPPED") for k in orig}
        rem_classes = {cmap.get(k, "UNMAPPED") for k in rem}
        check(orig_classes <= {"TARGET", "FUSION", "LINKER"},
              "I7 ORIGINAL contains only TARGET/FUSION/LINKER", f"got {orig_classes}")
        check(rem_classes <= {"TARGET"},
              "I7 REMOVED contains only TARGET", f"got {rem_classes}")

        tag_keys = set(out["tag_residues_deleted"])
        check(not (tag_keys & set(orig)), "I4 TAG absent from ORIGINAL",
              f"leaked {sorted(tag_keys & set(orig))[:5]}")
        check(not (tag_keys & set(rem)), "I4 TAG absent from REMOVED",
              f"leaked {sorted(tag_keys & set(rem))[:5]}")

        # ---- I5: no TARGET atom coordinate changes between conditions
        moved = []
        for key, atoms in rem.items():
            if key not in orig:
                moved.append((key, "missing_in_original"))
                continue
            for aname, xyz in atoms.items():
                if orig[key].get(aname) != xyz:
                    moved.append((key, aname))
        check(not moved, "I5 REMOVED alters no TARGET atom coordinate",
              f"{len(moved)} discrepancies e.g. {moved[:3]}")
        check(len(rem) == counts.get("TARGET", -1),
              "I5 REMOVED residue count == TARGET count",
              f"{len(rem)} vs {counts.get('TARGET')}")

        # ---- I6: numbering traceable back to the deposited file
        cif = os.path.join(prep.CACHE, f"{rec['pdb_id'].lower()}.cif")
        st = gemmi.read_structure(cif)
        st.setup_entities()
        dep = set()
        for ch in st[0]:
            if ch.name != rec["representative_chain"]:
                continue
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info is not None and info.is_amino_acid():
                    dep.add(f"{res.seqid.num}{res.seqid.icode.strip()}")
        check(set(orig) <= dep, "I6 ORIGINAL numbering is a subset of deposited auth numbering",
              f"novel: {sorted(set(orig) - dep)[:5]}")
        check(set(rem) <= dep, "I6 REMOVED numbering is a subset of deposited auth numbering")

        # ---- I8: detector-style residue IDs map back
        sample = sorted(orig)[:25]
        unmapped = [k for k in sample if k not in cmap]
        check(not unmapped, "I8 written residue IDs resolve in auth_class_map",
              f"unmapped {unmapped[:5]}")

        # ---- deletion boundaries exist where an excision happened
        if rec["topology"] == "INTERNAL_INSERTION":
            check(len(out["deletion_boundaries"]) >= 2,
                  "A3 internal insertion yields >=2 deletion boundary residues",
                  f"got {len(out['deletion_boundaries'])}")
        else:
            check(len(out["deletion_boundaries"]) >= 1,
                  "A3 terminal fusion yields >=1 deletion boundary residue",
                  f"got {len(out['deletion_boundaries'])}")

    print("\n" + "=" * 78)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} FAILED of {CHECKS[0]} checks")
        for f in FAILURES:
            print("  -", f)
        sys.exit(1)
    print(f"RESULT: all {CHECKS[0]} checks PASSED — mapping validated, pilot may proceed")


if __name__ == "__main__":
    main()
