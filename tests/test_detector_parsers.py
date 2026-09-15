#!/usr/bin/env python
"""
Phase 2C — detector output parser validation. MUST pass before the 96-run pilot.

Checks that the P2Rank and fpocket parsers recover the right pockets, in the right rank
order, with the right residues, and that every reported residue maps back onto the
entity-position class map.

Critical invariant D3: fpocket's `pockets/pocketN_*` files are matched to the correct
"Pocket N :" block in the info file. fpocket 4.2.3 numbers these files 1-based; an off-by-one
here silently shifts every pocket's residue set by one rank, which would corrupt every
downstream classification. Verified independently via the alpha-sphere count.

Run:  python tests/test_detector_parsers.py
"""
import importlib.util, json, os, re, subprocess, sys, shutil, math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(os.environ.get("TEMP", "/tmp"), "fusiontag_parser_test")
JAVA_HOME = r"C:\Program Files\Eclipse Adoptium\jdk-17.0.19.10-hotspot"
CASE_PDB = "5IU7"
CASE_CHAIN = "A"


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(ROOT, "scripts", filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


build = _load("build_pool", "01_build_pool.py")
prep = _load("prepare", "03_prepare_structures.py")
cls = _load("classify", "05_classify_pockets.py")
runner = _load("runner", "common_runner.py")

FAILURES, CHECKS = [], [0]


def check(cond, label, detail=""):
    CHECKS[0] += 1
    print(("    PASS  " if cond else "    FAIL  ") + label + ("  " + detail if detail and not cond else ""))
    if not cond:
        FAILURES.append(f"{label} {detail}")


def win_to_wsl(p):
    p = os.path.abspath(p).replace("\\", "/")
    return "/mnt/" + p[0].lower() + p[2:] if p[1:3] == ":/" else p


def main():
    print("=" * 78)
    print("PHASE 2C — DETECTOR PARSER VALIDATION")
    print("=" * 78)
    os.makedirs(WORK, exist_ok=True)

    # ---- prepare the case structure through the real pipeline
    raw = build.post(build.GQL, {"query": build.QUERY,
                                 "variables": {"ids": [f"{CASE_PDB}_1"]}})["data"]["polymer_entities"][0]
    raw["_partner_accession"], raw["_partner_tag"] = "P0ABE7", "BRIL"
    rec = build.evaluate(raw)
    out = prep.prepare_one(rec)
    cmap = out["auth_class_map"]

    base = f"{CASE_PDB}_{CASE_CHAIN}_ORIGINAL"
    inp = os.path.join(WORK, base + ".pdb")
    shutil.copyfile(os.path.join(ROOT, out["original_pdb"]), inp)

    # ---- run both detectors exactly as the pilot will
    print(f"\n--- running detectors on {base} ---")
    p2dir = os.path.join(WORK, "p2rank_out")
    subprocess.run(runner.p2rank_cmd(inp, p2dir), capture_output=True, text=True,
                   env=runner.p2rank_env(), timeout=1800, check=True)
    fargv, fenv = runner.fpocket_cmd(inp)
    subprocess.run(fargv, capture_output=True, text=True, env=fenv, timeout=1800, check=True)
    fpdir = os.path.join(WORK, base + "_out")

    p2 = cls.parse_p2rank(p2dir, base + ".pdb")
    fp = cls.parse_fpocket(fpdir, base)
    print(f"    P2Rank pockets: {len(p2)}   fpocket pockets: {len(fp)}")

    # ---- D1: pockets parsed at all
    check(len(p2) > 0, "D1 P2Rank pockets parsed")
    check(len(fp) > 0, "D1 fpocket pockets parsed")

    # ---- D2: ranks are 1..N, contiguous, ascending; scores non-increasing
    for name, pk in (("P2Rank", p2), ("fpocket", fp)):
        ranks = [p["rank"] for p in pk]
        check(ranks == list(range(1, len(pk) + 1)),
              f"D2 {name} ranks contiguous 1..N", f"got {ranks[:6]}")
        scores = [p["score"] for p in pk]
        check(all(scores[i] >= scores[i + 1] - 1e-9 for i in range(len(scores) - 1)),
              f"D2 {name} scores non-increasing with rank", f"got {scores[:5]}")

    # ---- D3: fpocket file<->info correspondence via independent alpha-sphere count
    info_spheres = {}
    cur = None
    for line in open(os.path.join(fpdir, base + "_info.txt"), encoding="utf-8"):
        m = re.match(r"^Pocket\s+(\d+)\s*:", line.strip())
        if m:
            cur = int(m.group(1))
        elif cur and "Number of Alpha Spheres" in line:
            info_spheres[cur] = int(float(line.split(":")[1].strip()))
    mismatch = []
    for pid, n_expected in sorted(info_spheres.items())[:10]:
        vert = os.path.join(fpdir, "pockets", f"pocket{pid}_vert.pqr")
        n_actual = sum(1 for l in open(vert, encoding="utf-8", errors="replace")
                       if l.startswith(("ATOM", "HETATM"))) if os.path.exists(vert) else -1
        if n_actual != n_expected:
            mismatch.append((pid, n_expected, n_actual))
    check(not mismatch,
          "D3 fpocket pocketN files correspond to 'Pocket N' info blocks (alpha-sphere count)",
          f"mismatches (pocket, info, file): {mismatch[:4]}")

    # ---- D4: every pocket has residues and a finite centroid
    for name, pk in (("P2Rank", p2), ("fpocket", fp)):
        empty = [p["rank"] for p in pk if not p["residues"]]
        check(not empty, f"D4 {name} every pocket has >=1 residue", f"empty ranks {empty[:6]}")
        nan = [p["rank"] for p in pk if any(math.isnan(c) for c in p["centroid"])]
        check(not nan, f"D4 {name} every pocket has a finite centroid", f"NaN ranks {nan[:6]}")

    # ---- D5: all residue IDs map back onto the entity class map
    for name, pk in (("P2Rank", p2), ("fpocket", fp)):
        allres = [r for p in pk for r in p["residues"]]
        miss = sorted({r for r in allres if r not in cmap})
        check(not miss, f"D5 {name} residue IDs resolve in auth_class_map",
              f"{len(miss)} unmapped e.g. {miss[:5]}")
        check(len(allres) > 0, f"D5 {name} produced residues at all")

    # ---- D6: both detectors can see fusion residues at all (mapping is two-sided)
    for name, pk in (("P2Rank", p2), ("fpocket", fp)):
        seen = {cmap.get(r) for p in pk for r in p["residues"]}
        check("TARGET" in seen, f"D6 {name} reports TARGET residues somewhere")
        check("FUSION" in seen, f"D6 {name} reports FUSION residues somewhere",
              f"classes seen: {seen}")

    # ---- D7: classification is total (no pocket falls through the rules)
    for name, pk in (("P2Rank", p2), ("fpocket", fp)):
        bad = []
        for p in pk:
            known = [cmap.get(r) for r in p["residues"]]
            known = [c for c in known if c in ("TARGET", "FUSION", "LINKER")]
            n = len(known)
            c = cls.classify(known.count("TARGET") / n if n else 0,
                             known.count("FUSION") / n if n else 0,
                             known.count("LINKER") / n if n else 0, n, 0.70)
            if c not in ("TARGET_DOMINATED", "FUSION_DOMINATED", "LINKER",
                         "INTERFACE", "MIXED", "UNASSIGNED"):
                bad.append((p["rank"], c))
        check(not bad, f"D7 {name} every pocket receives a valid class", str(bad[:4]))

    print("\n" + "=" * 78)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} FAILED of {CHECKS[0]} checks")
        for f in FAILURES:
            print("  -", f)
        sys.exit(1)
    print(f"RESULT: all {CHECKS[0]} checks PASSED — detector parsing validated")


if __name__ == "__main__":
    main()
