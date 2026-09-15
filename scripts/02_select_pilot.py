#!/usr/bin/env python
"""
Phase 2B step 2 — draw the frozen n=24 pilot.

Selection chain, exactly as protocol v1.1 section 9:
  enumerate -> inclusion/exclusion filters -> one structure per (partner, target accession)
  [highest resolution, PDB-ID tie break] -> lexicographic target ordering
  -> fixed seed 20260915 -> sample 8 per partner

Outcome-blind: no detector has been run; no structure has been visually inspected.
Seeding is per-partner (`Random(f"20260915-{PARTNER}")`) so that each stratum's draw does not
depend on the order in which partners happen to be processed.

Outputs:
  data_manifest/pilot_selection.json   selection record incl. the per-partner candidate lists
  data_manifest/collapsed_pool.tsv     the one-per-(partner,target) pool the draw sampled from
"""
import json, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
SEED_BASE = "20260915"
N_PER_PARTNER = 8
PARTNER_ORDER = ["BRIL", "T4L", "MBP"]


def main():
    data = json.load(open(os.path.join(MANI, "pool_eligibility.json")))
    recs = data["records"]
    eligible = [r for r in recs if r["eligible"]]
    print(f"eligible entities: {len(eligible)} / {len(recs)} candidates")

    # --- collapse: one structure per (partner, target accession)
    by_key = {}
    for r in eligible:
        key = (r["fusion_partner"], r["target_accession"])
        by_key.setdefault(key, []).append(r)

    collapsed = {}
    for key, group in by_key.items():
        # highest resolution == numerically smallest; tie broken by lowest PDB ID
        group_sorted = sorted(group, key=lambda r: (r["resolution"], r["pdb_id"]))
        winner = group_sorted[0]
        winner = dict(winner)
        winner["n_structures_for_this_target"] = len(group)
        winner["collapsed_alternatives"] = [g["pdb_id"] for g in group_sorted[1:]]
        collapsed[key] = winner

    print("\ncollapsed pool (one structure per partner x target):")
    for tag in PARTNER_ORDER:
        n = sum(1 for k in collapsed if k[0] == tag)
        print(f"  {tag:5s}: {n} distinct target accessions")

    # --- draw
    selection, draw_log = [], {}
    for tag in PARTNER_ORDER:
        accs = sorted(k[1] for k in collapsed if k[0] == tag)   # lexicographic ordering
        rng = random.Random(f"{SEED_BASE}-{tag}")
        if len(accs) < N_PER_PARTNER:
            raise SystemExit(f"FATAL: only {len(accs)} eligible targets for {tag}")
        drawn = sorted(rng.sample(accs, N_PER_PARTNER))
        draw_log[tag] = {"seed": f"{SEED_BASE}-{tag}",
                         "n_candidate_targets": len(accs),
                         "candidate_targets": accs,
                         "drawn_targets": drawn}
        for acc in drawn:
            selection.append(collapsed[(tag, acc)])
        print(f"  {tag}: drew 8 of {len(accs)} -> {drawn}")

    # cross-partner duplicate targets (allowed; strata are separate) — flagged for transparency
    tgt_counts = {}
    for s in selection:
        tgt_counts[s["target_accession"]] = tgt_counts.get(s["target_accession"], 0) + 1
    dupes = {k: v for k, v in tgt_counts.items() if v > 1}
    if dupes:
        print(f"\nNOTE: target accession(s) drawn in more than one partner stratum: {dupes}")

    out = {"protocol_version": "1.1",
           "seed_base": SEED_BASE,
           "n_per_partner": N_PER_PARTNER,
           "selection_rule": "eligible -> one per (partner,target) by (resolution asc, pdb_id asc)"
                             " -> lexicographic target sort -> Random(f'{seed}-{partner}')"
                             ".sample(targets, 8)",
           "source_pool_sha256": data["source_sha256"],
           "draw_log": draw_log,
           "selection": selection}
    json.dump(out, open(os.path.join(MANI, "pilot_selection.json"), "w"), indent=1)

    cols = ["fusion_partner", "target_accession", "pdb_id", "entity_id", "representative_chain",
            "resolution", "method", "topology", "n_structures_for_this_target"]
    with open(os.path.join(MANI, "collapsed_pool.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for key in sorted(collapsed):
            r = collapsed[key]
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    print(f"\nPILOT n = {len(selection)}")
    print(f"{'partner':6s} {'pdb':6s} {'chain':5s} {'res':>5s}  {'topology':18s} target")
    for s in selection:
        print(f"{s['fusion_partner']:6s} {s['pdb_id']:6s} {s['representative_chain']:5s} "
              f"{s['resolution']:5.2f}  {str(s['topology']):18s} {s['target_accession']} "
              f"{(s['description'] or '')[:48]}")


if __name__ == "__main__":
    main()
