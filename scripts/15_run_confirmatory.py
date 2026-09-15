#!/usr/bin/env python
"""
Confirmatory execution — protocol v1.3, frozen state
environment/PROTOCOL_FREEZE.json.

Drives the frozen pipeline over the 128-structure target-disjoint confirmatory set:
  prepare (ORIGINAL + FUSION-REMOVED) -> P2Rank 2.5.1 + fpocket 4.2.3 -> frozen classification.

No parameter tuning. No per-structure rescue. Technical failures stay in the manifest.

Variants (pre-specified structure-representation / HETATM sensitivities) are produced by
FUSIONTAG_VARIANT: "primary" (default), "assembly1", "ions".

Outputs under results/confirmatory/<variant>/.
"""
import collections, hashlib, importlib.util, json, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = os.environ.get("FUSIONTAG_VARIANT", "primary")
OUTDIR = os.path.join(ROOT, "results", "confirmatory", VARIANT)
PREPDIR = os.path.join(ROOT, "prepared_confirmatory", VARIANT)
WORK = os.path.join(r"C:\Users\user\AppData\Local\Temp\claude"
                    r"\C--AI-PROJECTS-Softwares-Fusion-Tag-Hazard"
                    r"\5c4a9bb4-b9dd-40f9-b5cc-b8b58139430b\scratchpad", "conf_" + VARIANT)
for d in (OUTDIR, PREPDIR, WORK):
    os.makedirs(d, exist_ok=True)


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


prep = _load("prep", "03_prepare_structures.py")
runner = _load("runner", "common_runner.py")
cls = _load("cls", "05_classify_pockets.py")
prep.PREP = PREPDIR


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_tree(p):
    h = hashlib.sha256()
    for root, dirs, files in os.walk(p):
        dirs.sort()
        for f in sorted(files):
            fp = os.path.join(root, f)
            h.update(os.path.relpath(fp, p).replace("\\", "/").encode())
            h.update(sha256_file(fp).encode())
    return h.hexdigest()


# ---------------------------------------------------------------- 1. prepare
def do_prepare(structures):
    out_path = os.path.join(MANI, f"confirmatory_manifest_{VARIANT}.json")
    if os.path.exists(out_path):
        sys.stderr.write("using cached prepared manifest\n")
        return json.load(open(out_path))["prepared"]
    out = []
    for i, rec in enumerate(structures, 1):
        try:
            out.append(prep.prepare_one(rec))
        except Exception as exc:
            out.append({"pdb_id": rec["pdb_id"], "entity_id": rec["entity_id"],
                        "fusion_partner": rec["fusion_partner"],
                        "target_accession": rec["target_accession"],
                        "PREPARATION_FAILED": f"{type(exc).__name__}: {exc}"})
            sys.stderr.write(f"  PREP FAIL {rec['pdb_id']}: {type(exc).__name__}: {exc}\n")
        if i % 20 == 0:
            sys.stderr.write(f"  prepared {i}/{len(structures)}\n")
    json.dump({"protocol_version": "1.3", "variant": VARIANT, "prepared": out},
              open(out_path, "w"), indent=1)
    return out


# ---------------------------------------------------------------- 2. run detectors
def do_run(prepared):
    out_path = os.path.join(OUTDIR, "detector_runs.json")
    if os.path.exists(out_path):
        sys.stderr.write("using cached detector runs\n")
        return json.load(open(out_path))["runs"]

    ok = [p for p in prepared if "PREPARATION_FAILED" not in p]
    staged = []
    for p in ok:
        for cond, key in (("ORIGINAL", "original_pdb"), ("REMOVED", "removed_pdb")):
            src = os.path.join(ROOT, p[key])
            base = f"{p['pdb_id']}_{p['auth_chain']}_{cond}"
            d = os.path.join(WORK, base)
            os.makedirs(d, exist_ok=True)
            dst = os.path.join(d, base + ".pdb")
            shutil.copyfile(src, dst)
            staged.append({"pdb_id": p["pdb_id"], "chain": p["auth_chain"],
                           "partner": p["fusion_partner"],
                           "target": p["target_accession"], "condition": cond,
                           "key": base, "input": dst, "input_sha256": sha256_file(dst)})
    sys.stderr.write(f"staged {len(staged)} condition-instances\n")
    runs = []

    # ---- fpocket, one WSL batch
    script = os.path.join(WORK, "_fp_batch.sh")
    lines = ["#!/bin/bash", "set -u", f'FP="{runner.FPOCKET_WSL}"', ""]
    for s in staged:
        wp = runner.win_to_wsl(s["input"])
        lines += [f'T0=$(date +%s.%N)',
                  f'"$FP" -f "{wp}" >/dev/null 2>"{wp}.fperr"', 'RC=$?',
                  f'T1=$(date +%s.%N)',
                  f'echo "###RC {s["key"]} $RC $(echo "$T1 - $T0" | bc)"']
    open(script, "w", newline="\n").write("\n".join(lines) + "\n")
    env = dict(os.environ, MSYS2_ARG_CONV_EXCL="*", MSYS_NO_PATHCONV="1")
    t0 = time.time()
    proc = subprocess.run(["wsl.exe", "-d", runner.WSL_DISTRO, "--", "bash",
                           runner.win_to_wsl(script)],
                          capture_output=True, text=True, env=env, timeout=14400)
    sys.stderr.write(f"fpocket batch done in {time.time()-t0:.0f}s\n")
    fp = {}
    for line in proc.stdout.replace("\r", "").splitlines():
        if line.startswith("###RC "):
            _, k, rc, rt = line.split()
            fp[k] = (int(rc), round(float(rt), 2))
    for s in staged:
        rc, rt = fp.get(s["key"], (-999, None))
        d = os.path.join(os.path.dirname(s["input"]), s["key"] + "_out")
        good = rc == 0 and os.path.isdir(d)
        rec = dict(s, detector="fpocket", detector_version=runner.FPOCKET_VERSION,
                   command=f'{runner.FPOCKET_WSL} -f {runner.win_to_wsl(s["input"])}',
                   returncode=rc, runtime_s=rt, output_dir=d if good else None,
                   output_sha256=sha256_tree(d) if good else None, success=bool(good))
        ef = s["input"] + ".fperr"
        if os.path.exists(ef):
            rec["stderr_tail"] = open(ef, encoding="utf-8", errors="replace").read()[-2000:]
        runs.append(rec)
    sys.stderr.write(f"  fpocket ok {sum(1 for r in runs if r['success'])}/{len(staged)}\n")

    # ---- P2Rank, one JVM per run
    for i, s in enumerate(staged, 1):
        od = os.path.join(os.path.dirname(s["input"]), "p2rank_out")
        cmd = runner.p2rank_cmd(s["input"], od)
        t0 = time.time()
        try:
            pr = subprocess.run(cmd, capture_output=True, text=True,
                                env=runner.p2rank_env(), timeout=3600)
            rc, so, se = pr.returncode, pr.stdout, pr.stderr
        except subprocess.TimeoutExpired:
            rc, so, se = -998, "", "TIMEOUT"
        rt = round(time.time() - t0, 2)
        pred = os.path.join(od, os.path.basename(s["input"]) + "_predictions.csv")
        good = rc == 0 and os.path.exists(pred)
        runs.append(dict(s, detector="p2rank", detector_version=runner.P2RANK_VERSION,
                         command=" ".join(cmd), returncode=rc, runtime_s=rt,
                         output_dir=od if good else None,
                         output_sha256=sha256_tree(od) if good else None,
                         success=bool(good),
                         stdout_tail=so[-1500:], stderr_tail=se[-1500:]))
        if i % 25 == 0 or not good:
            sys.stderr.write(f"  p2rank {i}/{len(staged)} rc={rc} {rt}s ok={good}\n")
    sys.stderr.write(f"  p2rank ok "
                     f"{sum(1 for r in runs if r['detector']=='p2rank' and r['success'])}"
                     f"/{len(staged)}\n")

    # ---- archive raw outputs
    raw = os.path.join(OUTDIR, "raw_detector")
    for r in runs:
        if not r["success"]:
            continue
        dest = os.path.join(raw, r["detector"], r["key"])
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(r["output_dir"], dest)
        r["archived_output"] = os.path.relpath(dest, ROOT)

    json.dump({"protocol_version": "1.3", "variant": VARIANT,
               "expected_runs": len(ok) * 4, "actual_runs": len(runs), "runs": runs},
              open(out_path, "w"), indent=1)
    return runs


# ---------------------------------------------------------------- 3. classify
def do_classify(prepared, runs):
    man = {p["pdb_id"] + "_" + p.get("auth_chain", ""): p
           for p in prepared if "PREPARATION_FAILED" not in p}
    rows = []
    for run in runs:
        if not run["success"]:
            continue
        meta = man[f"{run['pdb_id']}_{run['chain']}"]
        od = os.path.join(ROOT, run["archived_output"])
        base = run["key"]
        pockets = (cls.parse_p2rank(od, base + ".pdb") if run["detector"] == "p2rank"
                   else cls.parse_fpocket(od, base))
        cmap = meta["auth_class_map"]
        refset = set(meta["reference_site_residues"])
        bnds = meta["deletion_boundaries"]
        import math
        for p in pockets:
            classes = [cmap.get(r, "UNMAPPED") for r in p["residues"]]
            known = [c for c in classes if c in ("TARGET", "FUSION", "LINKER")]
            n = len(known)
            ft = known.count("TARGET") / n if n else 0.0
            ff = known.count("FUSION") / n if n else 0.0
            fl = known.count("LINKER") / n if n else 0.0
            dmin = None
            if run["condition"] == "REMOVED" and bnds and not math.isnan(p["centroid"][0]):
                dmin = min(math.dist(p["centroid"], (b["x"], b["y"], b["z"])) for b in bnds)
            inter = refset & set(p["residues"])
            rf = len(inter) / len(refset) if refset else None
            row = {"pdb_id": run["pdb_id"], "chain": run["chain"],
                   "fusion_partner": run["partner"], "target_accession": meta["target_accession"],
                   "topology": meta["topology"], "resolution": meta["resolution"],
                   "detector": run["detector"], "condition": run["condition"],
                   "rank": p["rank"], "score": p["score"], "extra_score": p.get("extra_score"),
                   "centroid_x": p["centroid"][0], "centroid_y": p["centroid"][1],
                   "centroid_z": p["centroid"][2],
                   "n_residues": len(p["residues"]), "n_classifiable": n,
                   "f_target": round(ft, 4), "f_fusion": round(ff, 4), "f_linker": round(fl, 4),
                   "dist_to_deletion_boundary": None if dmin is None else round(dmin, 2),
                   "ref_site_overlap_n": len(inter),
                   "ref_site_overlap_frac": None if rf is None else round(rf, 4),
                   "recovers_reference_site": bool(refset) and len(inter) >= 3
                                              and rf >= 0.25,
                   "residues": ";".join(p["residues"])}
            for t in cls.SENSITIVITY_THRESHOLDS:
                row[f"class_at_{t:.2f}"] = cls.classify(ft, ff, fl, n, t)
            row["pocket_class"] = row["class_at_0.70"]
            row["fusion_associated"] = row["pocket_class"] in ("FUSION_DOMINATED", "INTERFACE",
                                                               "LINKER")
            rows.append(row)
    rows.sort(key=lambda r: (r["detector"], r["pdb_id"], r["condition"], r["rank"]))
    json.dump({"protocol_version": "1.3", "variant": VARIANT,
               "primary_dominance": cls.PRIMARY_DOMINANCE, "pockets": rows},
              open(os.path.join(OUTDIR, "pockets_classified.json"), "w"), indent=1)
    cols = [c for c in rows[0].keys() if c != "residues"] + ["residues"]
    with open(os.path.join(OUTDIR, "pockets_classified.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    return rows


def main():
    frz = json.load(open(os.path.join(ROOT, "environment", "PROTOCOL_FREEZE.json")))
    print(f"protocol freeze: {frz['freeze']}  sha256(PROTOCOL.md)="
          f"{frz['sha256']['PROTOCOL.md'][:16]}...  variant={VARIANT}")
    structures = json.load(open(os.path.join(MANI,
                                             "confirmatory_set_final.json")))["structures"]
    print(f"confirmatory structures: {len(structures)}")
    prepared = do_prepare(structures)
    ok = [p for p in prepared if "PREPARATION_FAILED" not in p]
    print(f"prepared {len(ok)}/{len(prepared)}")
    for f in [p for p in prepared if "PREPARATION_FAILED" in p]:
        print(f"  PREPARATION FAILURE (retained in manifest): {f['pdb_id']} "
              f"{f['PREPARATION_FAILED']}")
    runs = do_run(prepared)
    print(f"runs {len(runs)} (expected {len(ok)*4})")
    for det in ("p2rank", "fpocket"):
        for cond in ("ORIGINAL", "REMOVED"):
            sub = [r for r in runs if r["detector"] == det and r["condition"] == cond]
            print(f"  {det:8s} {cond:9s}: {sum(1 for r in sub if r['success'])}/{len(sub)}")
    fails = [r for r in runs if not r["success"]]
    if fails:
        print("TECHNICAL FAILURES (retained, not substituted):")
        for f in fails:
            print(f"  {f['detector']} {f['key']} rc={f['returncode']}")
    rows = do_classify(prepared, runs)
    print(f"classified {len(rows)} pockets")
    cov = sum(1 for p in ok if p["reference_site_available"])
    print(f"reference-site coverage: {cov}/{len(ok)} = {100*cov/len(ok):.1f}%")


if __name__ == "__main__":
    main()
