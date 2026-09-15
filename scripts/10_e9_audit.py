#!/usr/bin/env python
"""
Phase 3A.2 — audit of exclusion criterion E9 (segment disorder).

E9 excluded 1,223 / 1,971 candidate entities and therefore materially determines the study
population. This script characterises it and independently verifies its implementation against
deposited coordinates on an outcome-blind random sample.

Implementation under audit (protocol v1.1 section 3.2, as coded in 01_build_pool.py):

    fusion_disorder = (|FUS_aligned| - |FUS_aligned AND modelled|) / |FUS_aligned|
    target_disorder = (|TGT_aligned| - |TGT_aligned AND modelled|) / |TGT_aligned|
    EXCLUDE if fusion_disorder > 0.20  OR  target_disorder > 0.20

  where *_aligned are SIFTS-aligned entity (label_seq) positions for that accession, and
  `modelled` = all entity positions MINUS positions reported by the RCSB instance feature
  UNOBSERVED_RESIDUE_XYZ for the representative chain.

Outputs:
  results/e9_audit/e9_distributions.tsv
  results/e9_audit/e9_manual_validation.tsv
  results/e9_audit/E9_AUDIT.md data appendix (printed)
"""
import collections, json, math, os, random, sys, urllib.request, time
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
OUT = os.path.join(ROOT, "results", "e9_audit")
CACHE = os.path.join(ROOT, ".structure_cache")
os.makedirs(OUT, exist_ok=True)
VALIDATION_N = 30
SEED = "E9AUDIT-20260915"


def pct(xs, q):
    if not xs:
        return None
    s = sorted(xs)
    return round(s[min(len(s) - 1, int(q * len(s)))], 3)


def summarise(name, xs):
    if not xs:
        return f"  {name:28s} n=0"
    s = sorted(xs)
    return (f"  {name:28s} n={len(s):4d}  min={s[0]:.2f}  p25={pct(s,.25):.2f}  "
            f"median={pct(s,.5):.2f}  p75={pct(s,.75):.2f}  p90={pct(s,.90):.2f}  max={s[-1]:.2f}")


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
                time.sleep(3 * (a + 1))
    return p


def expand(ranges):
    s = set()
    for a, b in ranges:
        s.update(range(a, b + 1))
    return s


def main():
    data = json.load(open(os.path.join(MANI, "pool_eligibility.json")))
    recs = data["records"]

    def has(r, code):
        return any(e.startswith(code) for e in r["exclusions"])

    e9 = [r for r in recs if has(r, "E9")]
    non_e9 = [r for r in recs if not has(r, "E9")]
    print("=" * 100)
    print("E9 AUDIT — segment-disorder exclusion")
    print("=" * 100)
    print(f"\ncandidates: {len(recs)}   E9-excluded: {len(e9)}   not E9: {len(non_e9)}")

    # ---- counts by partner
    print("\nE9 exclusions by partner:")
    for pt in ("BRIL", "T4L", "MBP"):
        sub = [r for r in recs if r["fusion_partner"] == pt]
        s9 = [r for r in sub if has(r, "E9")]
        ft = sum(1 for r in s9 if any(e.startswith("E9_fusion") for e in r["exclusions"]))
        tt = sum(1 for r in s9 if any(e.startswith("E9_target") for e in r["exclusions"]))
        print(f"  {pt:5s}: {len(s9):4d}/{len(sub):4d} = {100*len(s9)/len(sub):4.1f}%  "
              f"(fusion-side {ft}, target-side {tt}, both {ft+tt-len(s9)})")

    # ---- E9 as the SOLE reason (the decision-relevant number)
    sole = [r for r in e9 if all(e.startswith("E9") for e in r["exclusions"])]
    print(f"\nE9 as the ONLY exclusion reason: {len(sole)} entities "
          f"({len(sole)/len(e9)*100:.1f}% of E9 exclusions)")
    print("  -> these are the entities that E9 alone removes from the study population")
    for pt in ("BRIL", "T4L", "MBP"):
        ss = [r for r in sole if r["fusion_partner"] == pt]
        tg = {r["target_accession"] for r in ss}
        print(f"     {pt:5s}: {len(ss):4d} entities across {len(tg)} distinct target accessions")

    # ---- overlap with other reasons
    print("\noverlap of E9 with other exclusion codes (entities failing both):")
    co = collections.Counter()
    for r in e9:
        for e in r["exclusions"]:
            c = e.split("=")[0]
            if not c.startswith("E9"):
                co[c] += 1
    for k, v in co.most_common():
        print(f"  {k:38s} {v:5d}  ({100*v/len(e9):.1f}% of E9-excluded)")

    # ---- distributions
    print("\ndisorder-fraction distributions:")
    print(summarise("fusion_disorder (all)", [r["fusion_disorder"] for r in recs
                                              if r.get("fusion_disorder") is not None]))
    print(summarise("target_disorder (all)", [r["target_disorder"] for r in recs
                                              if r.get("target_disorder") is not None]))
    print(summarise("target_disorder (E9 only)", [r["target_disorder"] for r in e9
                                                  if r.get("target_disorder") is not None]))
    print(summarise("target_disorder (passed)", [r["target_disorder"] for r in non_e9
                                                 if r.get("target_disorder") is not None]))

    # how far above the line are the failures?
    over = [r["target_disorder"] for r in e9
            if r.get("target_disorder") is not None and r["target_disorder"] > 0.20]
    near = sum(1 for x in over if x <= 0.30)
    print(f"\n  of {len(over)} target-side E9 failures: {near} ({100*near/max(len(over),1):.0f}%) "
          f"fall between 0.20 and 0.30 — i.e. marginally over the line")
    for thr in (0.25, 0.30, 0.35, 0.40, 0.50):
        n = sum(1 for r in recs if (r.get("target_disorder") or 0) <= thr
                and (r.get("fusion_disorder") or 0) <= thr)
        print(f"  entities passing a disorder threshold of {thr:.2f}: {n}")

    # ---- length / method / resolution / topology comparisons
    print("\nlength distributions, pass vs E9-fail:")
    print(summarise("target_modelled  (pass)", [r["target_modelled"] for r in non_e9
                                                if r.get("target_modelled") is not None]))
    print(summarise("target_modelled  (E9)", [r["target_modelled"] for r in e9
                                              if r.get("target_modelled") is not None]))
    print(summarise("fusion_modelled  (pass)", [r["fusion_modelled"] for r in non_e9
                                                if r.get("fusion_modelled") is not None]))
    print(summarise("fusion_modelled  (E9)", [r["fusion_modelled"] for r in e9
                                              if r.get("fusion_modelled") is not None]))
    print(summarise("resolution       (pass)", [r["resolution"] for r in non_e9 if r.get("resolution")]))
    print(summarise("resolution       (E9)", [r["resolution"] for r in e9 if r.get("resolution")]))

    print("\nexperimental method, pass vs E9-fail:")
    for grp, name in ((non_e9, "pass"), (e9, "E9-fail")):
        c = collections.Counter(r["method"] for r in grp)
        tot = sum(c.values())
        print(f"  {name:8s} " + "  ".join(f"{k}={v} ({100*v/tot:.0f}%)"
                                          for k, v in c.most_common(3)))

    print("\ntopology, pass vs E9-fail:")
    for grp, name in ((non_e9, "pass"), (e9, "E9-fail")):
        c = collections.Counter(r["topology"] for r in grp)
        tot = sum(c.values())
        print(f"  {name:8s} " + "  ".join(f"{k}={v} ({100*v/tot:.0f}%)"
                                          for k, v in c.most_common(5)))
    # differential effect
    print("\n  E9 failure rate within each topology:")
    for topo in ("INTERNAL_INSERTION", "TERMINAL_N", "TERMINAL_C", "COMPLEX", "None"):
        sub = [r for r in recs if str(r["topology"]) == topo]
        if not sub:
            continue
        k = sum(1 for r in sub if has(r, "E9"))
        print(f"    {topo:20s} {k}/{len(sub)} = {100*k/len(sub):.1f}%")

    with open(os.path.join(OUT, "e9_distributions.tsv"), "w", encoding="utf-8") as fh:
        cols = ["entity_id", "fusion_partner", "target_accession", "topology", "method",
                "resolution", "entity_length", "target_modelled", "fusion_modelled",
                "target_disorder", "fusion_disorder", "align_coverage_of_modelled",
                "e9_excluded", "e9_sole_reason", "exclusions"]
        fh.write("\t".join(cols) + "\n")
        soleset = {r["entity_id"] for r in sole}
        for r in recs:
            row = dict(r)
            row["e9_excluded"] = has(r, "E9")
            row["e9_sole_reason"] = r["entity_id"] in soleset
            row["exclusions"] = ";".join(r["exclusions"])
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")

    # ------------------------------------------------------------------ manual validation
    print("\n" + "=" * 100)
    print(f"INDEPENDENT VALIDATION — {VALIDATION_N} E9-excluded entities, "
          f"outcome-blind random sample (seed {SEED})")
    print("Recomputes disorder directly from deposited coordinates, bypassing the RCSB")
    print("UNOBSERVED_RESIDUE_XYZ feature that the pipeline used.")
    print("=" * 100)
    rng = random.Random(SEED)
    # Validate only on entities E9 can meaningfully act on: genuinely chimeric, both segments
    # present. Entities with an empty target segment receive target_disorder = 1.0 from the
    # `else 1.0` fallback in 01_build_pool.py and are E1 exclusions, not E9 exclusions.
    pool = sorted([r for r in e9
                   if len(r["uniprot_accessions"]) == 2 and r["target_accession"]
                   and r["target_ranges"] and r["fusion_ranges"]],
                  key=lambda r: r["entity_id"])
    print(f"  (validation pool restricted to {len(pool)} genuinely chimeric E9-excluded "
          f"entities, of {len(e9)} E9-flagged)")
    sample = rng.sample(pool, min(VALIDATION_N, len(pool)))

    vrows, agree, disagree = [], 0, 0
    for i, r in enumerate(sample, 1):
        try:
            cif = fetch_cif(r["pdb_id"])
            st = gemmi.read_structure(cif)
            st.setup_entities()
            chain = r["representative_chain"]
            modelled = set()
            auth_seen = []
            for ch in st[0]:
                if ch.name != chain:
                    continue
                for res in ch:
                    info = gemmi.find_tabulated_residue(res.name)
                    if info is not None and info.is_amino_acid() and res.label_seq is not None:
                        modelled.add(int(res.label_seq))
                        auth_seen.append(res.seqid.num)
            tgt = expand(r["target_ranges"]); fus = expand(r["fusion_ranges"])
            td = (len(tgt) - len(tgt & modelled)) / len(tgt) if tgt else None
            fd = (len(fus) - len(fus & modelled)) / len(fus) if fus else None
            ok = (td is not None and abs(td - r["target_disorder"]) < 0.02) and \
                 (fd is not None and abs(fd - r["fusion_disorder"]) < 0.02)
            agree += ok; disagree += (not ok)

            # where do the missing target residues sit?
            miss = sorted(tgt - modelled)
            lo, hi = (min(tgt), max(tgt)) if tgt else (0, 0)
            n_term = sum(1 for p in miss if p < lo + 0.10 * (hi - lo + 1))
            c_term = sum(1 for p in miss if p > hi - 0.10 * (hi - lo + 1))
            internal = len(miss) - n_term - c_term
            # author-numbering gaps that are NOT label_seq gaps (numbering-gap confusion check)
            auth_gaps = sum(1 for a, b in zip(sorted(auth_seen), sorted(auth_seen)[1:])
                            if b - a > 1)
            vrows.append({
                "entity_id": r["entity_id"], "partner": r["fusion_partner"],
                "topology": r["topology"], "chain": chain,
                "api_target_disorder": r["target_disorder"],
                "coord_target_disorder": None if td is None else round(td, 4),
                "api_fusion_disorder": r["fusion_disorder"],
                "coord_fusion_disorder": None if fd is None else round(fd, 4),
                "agreement": ok,
                "n_missing_target": len(miss),
                "missing_at_segment_N_term": n_term,
                "missing_internal": internal,
                "missing_at_segment_C_term": c_term,
                "author_numbering_gaps": auth_gaps,
                "unaligned_construct_residues": len(r.get("tag_positions", []))
                                                + len(r.get("linker_positions", [])),
            })
            print(f"  [{i:2d}/{len(sample)}] {r['entity_id']:8s} {r['fusion_partner']:5s} "
                  f"api tgt={r['target_disorder']:.3f}/fus={r['fusion_disorder']:.3f}  "
                  f"coord tgt={td:.3f}/fus={fd:.3f}  {'OK' if ok else 'MISMATCH'}  "
                  f"missing N/int/C = {n_term}/{internal}/{c_term}")
        except Exception as exc:
            print(f"  [{i:2d}] {r['entity_id']} FAILED: {type(exc).__name__}: {exc}")
            vrows.append({"entity_id": r["entity_id"], "agreement": None,
                          "error": f"{type(exc).__name__}: {exc}"})

    cols = sorted({k for v in vrows for k in v})
    with open(os.path.join(OUT, "e9_manual_validation.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for v in vrows:
            fh.write("\t".join(str(v.get(c, "")) for c in cols) + "\n")

    print(f"\n  implementation agreement: {agree}/{agree+disagree} "
          f"({100*agree/max(agree+disagree,1):.0f}%)")
    ok_rows = [v for v in vrows if v.get("agreement") is not None]
    if ok_rows:
        tot_miss = sum(v["n_missing_target"] for v in ok_rows)
        tn = sum(v["missing_at_segment_N_term"] for v in ok_rows)
        ti = sum(v["missing_internal"] for v in ok_rows)
        tc = sum(v["missing_at_segment_C_term"] for v in ok_rows)
        print(f"  location of missing target residues across the sample: "
              f"N-term {tn} ({100*tn/max(tot_miss,1):.0f}%), "
              f"internal {ti} ({100*ti/max(tot_miss,1):.0f}%), "
              f"C-term {tc} ({100*tc/max(tot_miss,1):.0f}%)")
        ng = sum(1 for v in ok_rows if v["author_numbering_gaps"] > 0)
        print(f"  entities with author-numbering gaps (potential numbering-gap confusion): "
              f"{ng}/{len(ok_rows)} — note the pipeline uses label_seq, not author numbering, "
              f"so these cannot be mistaken for disorder")


if __name__ == "__main__":
    main()
