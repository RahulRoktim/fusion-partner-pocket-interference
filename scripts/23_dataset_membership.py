#!/usr/bin/env python
"""
Phase 5 STEP 1 — freeze dataset versions and canonical membership.

Canonical membership comes from the published `.ds` files, NOT from directory listings. The
p2rank-datasets README is explicit that a directory may hold more PDB files than the dataset
defines (holo4k/ has 4543 files; holo4k.ds defines 4009 proteins, and 4009 is the published n).

Roles are taken from the P2Rank paper and the dataset repository README, recorded verbatim in
the output so the provenance claim is auditable.

Outputs: dataset_audit/dataset_registry.json
"""
import hashlib, json, os, re, sys, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD = os.path.join(ROOT, "dataset_audit")
MEM = os.path.join(AUD, "membership")
ACCESS_DATE = "2026-09-15"
BASE = "https://github.com/rdk/p2rank-datasets"

# role, evidence for the role, identifier granularity
SPEC = {
    "CHEN11": {
        "file": "chen11.ds", "role": "TRAIN",
        "role_evidence": "P2Rank (Krivak & Hoksza 2018): the distributed default model is "
                         "trained on CHEN11. Repository README lists CHEN11 as a main protein "
                         "set used 'for training and evaluation'.",
        "granularity": "PDB entry + chain"},
    "JOINED": {
        "file": "joined.ds", "role": "DEVELOPMENT/VALIDATION",
        "role_evidence": "P2Rank (Krivak & Hoksza 2018): JOINED (ASTEX + B48 + U48 + DT198 + "
                         "B210) is used for model selection / hyperparameter choice.",
        "granularity": "PDB entry"},
    "COACH420": {
        "file": "coach420.ds", "role": "TEST",
        "role_evidence": "P2Rank (Krivak & Hoksza 2018): COACH420 is an evaluation/test set.",
        "granularity": "PDB entry + chain"},
    "HOLO4K": {
        "file": "holo4k.ds", "role": "TEST",
        "role_evidence": "P2Rank (Krivak & Hoksza 2018): HOLO4K is an evaluation/test set. "
                         "README: 'Disjunct with CHEN11 and JOINED.'",
        "granularity": "PDB entry"},
    "FPTRAIN": {
        "file": "fptrain.ds", "role": "TRAIN (fpocket scoring function)",
        "role_evidence": "p2rank-datasets README, verbatim: 'FPTRAIN: dataset used by Fpocket "
                         "for training its pocket scoring function'.",
        "granularity": "PDB entry"},
    "COACH420_mlig": {
        "file": "coach420(mlig).ds", "role": "TEST (ground-truth ligand codes)",
        "role_evidence": "README: '(mlig) datasets contain explicitly specified relevant "
                         "ligands. Valid ligand codes come from MOAD 2013 database.'",
        "granularity": "PDB entry + chain + ligand codes"},
    "HOLO4K_mlig": {
        "file": "holo4k(mlig).ds", "role": "TEST (ground-truth ligand codes)",
        "role_evidence": "README, as above.", "granularity": "PDB entry + ligand codes"},
    "JOINED_mlig": {
        "file": "joined(mlig).ds", "role": "DEVELOPMENT/VALIDATION (ground-truth ligand codes)",
        "role_evidence": "README, as above.", "granularity": "PDB entry + ligand codes"},
}


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def parse(path):
    """Returns list of {pdb_id, chain, subset, ligands, raw}."""
    out = []
    for line in open(path, encoding="utf-8", errors="replace"):
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("PARAM.") or s.startswith("HEADER:"):
            continue
        parts = s.split()
        rel = parts[0]
        ligs = parts[1].split(",") if len(parts) > 1 else []
        stem = os.path.basename(rel)
        if not stem.lower().endswith(".pdb"):
            continue
        stem = stem[:-4]
        subset = os.path.dirname(rel).split("/")[-1] if "/" in os.path.dirname(rel) else None
        # chen11 style: a.001.001.001_1s69a  -> pdb 1s69 chain a
        m = re.match(r"^[a-z]\.[\d.]+_([0-9a-zA-Z]{4})([A-Za-z0-9]?)$", stem)
        if m:
            pdb, ch = m.group(1), (m.group(2) or None)
        elif re.fullmatch(r"[0-9a-zA-Z]{4}_[A-Za-z0-9]", stem):
            pdb, ch = stem[:4], stem[5]          # DT198 style: 1dug_A
        elif re.fullmatch(r"[0-9a-zA-Z]{4}[A-Za-z0-9]", stem):
            pdb, ch = stem[:4], stem[4]
        elif re.fullmatch(r"[0-9a-zA-Z]{4}", stem):
            pdb, ch = stem, None
        else:
            pdb, ch = stem[:4], None
        out.append({"pdb_id": pdb.upper(), "chain": ch, "subset": subset,
                    "ligands": ligs, "raw": rel})
    return out


def main():
    reg = {"frozen": ACCESS_DATE, "source_repository": BASE,
           "canonical_membership_note":
               "Membership taken from the published .ds files, not directory listings. "
               "The repository README warns that directories contain more PDB files than the "
               "datasets define (holo4k/ has 4543 files vs holo4k.ds 4009 entries).",
           "datasets": {}}
    for name, sp in SPEC.items():
        path = os.path.join(MEM, sp["file"])
        if not os.path.exists(path):
            print(f"  MISSING {sp['file']}")
            continue
        members = parse(path)
        entries = {m["pdb_id"] for m in members}
        chains = {(m["pdb_id"], m["chain"]) for m in members if m["chain"]}
        subs = {}
        for m in members:
            if m["subset"]:
                subs[m["subset"]] = subs.get(m["subset"], 0) + 1
        reg["datasets"][name] = {
            "dataset_name": name, "version": "p2rank-datasets master @ " + ACCESS_DATE,
            "source_url": f"{BASE}/blob/master/{sp['file']}",
            "membership_file": sp["file"],
            "sha256": sha256(path), "access_date": ACCESS_DATE,
            "role": sp["role"], "role_evidence": sp["role_evidence"],
            "identifier_granularity": sp["granularity"],
            "n_members": len(members), "n_unique_pdb_entries": len(entries),
            "n_unique_chains": len(chains) or None,
            "constituent_subsets": subs or None,
            "members": members}
        print(f"{name:16s} role={sp['role']:34s} members={len(members):5d} "
              f"entries={len(entries):5d} chains={len(chains) or '-':>5} "
              f"sha256={reg['datasets'][name]['sha256'][:12]}")
        if subs:
            print(f"                 constituents: {subs}")

    allpdb = sorted({m["pdb_id"] for d in reg["datasets"].values() for m in d["members"]})
    reg["all_unique_pdb_entries"] = allpdb
    reg["n_all_unique_pdb_entries"] = len(allpdb)
    print(f"\nunique PDB entries across all frozen membership lists: {len(allpdb)}")
    json.dump(reg, open(os.path.join(AUD, "dataset_registry.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
