#!/usr/bin/env python
"""
Amendment J — prepare, run and classify ONLY the newly admitted E9 = 0.35 structures, then merge
with the frozen primary outputs to produce a complete sensitivity result set.

Structures already in the E9 = 0.20 primary cohort are NEVER re-run: their prepared inputs,
detector runs and classified pockets are copied verbatim from the frozen primary artifacts.

Emits artifacts in exactly the layout the frozen analysis script expects, so
16_confirmatory_analysis.py runs unchanged with FUSIONTAG_VARIANT=e9_035:
  data_manifest/confirmatory_manifest_e9_035.json
  results/confirmatory/e9_035/detector_runs.json
  results/confirmatory/e9_035/pockets_classified.json / .tsv
"""
import collections, hashlib, importlib.util, json, math, os, shutil, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
PRIM = os.path.join(ROOT, "results", "confirmatory", "primary")
OUT = os.path.join(ROOT, "results", "confirmatory", "e9_035")
PREPDIR = os.path.join(ROOT, "prepared_confirmatory", "e9_035")
WORK_ROOT = os.path.abspath(os.path.expanduser(
    os.environ.get("FUSIONTAG_WORK") or os.path.join(tempfile.gettempdir(), "fusiontag_hazard")
))
WORK = os.path.join(WORK_ROOT, "e9_035")
for d in (OUT, PREPDIR, WORK):
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
    for r, dirs, files in os.walk(p):
        dirs.sort()
        for f in sorted(files):
            fp = os.path.join(r, f)
            h.update(os.path.relpath(fp, p).replace("\\", "/").encode())
            h.update(sha256_file(fp).encode())
    return h.hexdigest()


def main():
    sens = json.load(open(os.path.join(MANI, "e9_035_set_final.json")))["structures"]
    plan = json.load(open(os.path.join(MANI, "e9_035_run_plan.json")))
    reuse_ids, new_ids = set(plan["reuse"]), list(plan["new_runs"])
    prim_man = {p["pdb_id"]: p for p in
                json.load(open(os.path.join(MANI,
                                            "confirmatory_manifest_primary.json")))["prepared"]
                if "PREPARATION_FAILED" not in p}
    prim_runs = json.load(open(os.path.join(PRIM, "detector_runs.json")))["runs"]
    prim_pockets = json.load(open(os.path.join(PRIM, "pockets_classified.json")))["pockets"]

    print(f"sensitivity cohort: {len(sens)} structures")
    print(f"  reuse frozen primary output for : {len(reuse_ids)}")
    print(f"  NEW structures to run           : {len(new_ids)} -> {new_ids}")
    print(f"  expected NEW detector runs      : {len(new_ids)*4}")

    # ---------------- prepare new structures with the frozen preparation rules
    cache = os.path.join(MANI, "confirmatory_manifest_e9_035.json")
    new_prepared = []
    if os.path.exists(cache):
        print("cached merged manifest present; reusing")
        merged = json.load(open(cache))["prepared"]
        new_prepared = [p for p in merged if p["pdb_id"] in set(new_ids)]
    else:
        by_pdb = {s["pdb_id"]: s for s in sens}
        for i, pid in enumerate(new_ids, 1):
            try:
                new_prepared.append(prep.prepare_one(by_pdb[pid]))
            except Exception as exc:
                new_prepared.append({"pdb_id": pid,
                                     "entity_id": by_pdb[pid]["entity_id"],
                                     "fusion_partner": by_pdb[pid]["fusion_partner"],
                                     "target_accession": by_pdb[pid]["target_accession"],
                                     "PREPARATION_FAILED": f"{type(exc).__name__}: {exc}"})
                sys.stderr.write(f"  PREP FAIL {pid}: {exc}\n")
            sys.stderr.write(f"  prepared {i}/{len(new_ids)} {pid}\n")
        merged = [prim_man[p] for p in sorted(reuse_ids)] + new_prepared
        json.dump({"protocol_version": "1.3", "variant": "e9_035",
                   "note": "125 records copied verbatim from the frozen primary manifest; "
                           "16 newly prepared.",
                   "prepared": merged}, open(cache, "w"), indent=1)
    okprep = [p for p in new_prepared if "PREPARATION_FAILED" not in p]
    print(f"prepared new: {len(okprep)}/{len(new_prepared)}")

    # ---------------- run detectors on the new structures only
    runs_path = os.path.join(OUT, "detector_runs.json")
    if os.path.exists(runs_path):
        print("cached detector runs present; reusing")
        allruns = json.load(open(runs_path))["runs"]
        newruns = [r for r in allruns if r["pdb_id"] in set(new_ids)]
    else:
        staged = []
        for p in okprep:
            for cond, key in (("ORIGINAL", "original_pdb"), ("REMOVED", "removed_pdb")):
                base = f"{p['pdb_id']}_{p['auth_chain']}_{cond}"
                d = os.path.join(WORK, base)
                os.makedirs(d, exist_ok=True)
                dst = os.path.join(d, base + ".pdb")
                shutil.copyfile(os.path.join(ROOT, p[key]), dst)
                staged.append({"pdb_id": p["pdb_id"], "chain": p["auth_chain"],
                               "partner": p["fusion_partner"],
                               "target": p["target_accession"], "condition": cond,
                               "key": base, "input": dst,
                               "input_sha256": sha256_file(dst)})
        print(f"staged {len(staged)} condition-instances")
        newruns = []

        # fpocket, one WSL batch
        script = os.path.join(WORK, "_fp_batch.sh")
        lines = ["#!/bin/bash", "set -u", f'FP="{runner.FPOCKET_WSL}"', ""]
        for s in staged:
            wp = runner.win_to_wsl(s["input"])
            lines += ['T0=$(date +%s.%N)',
                      f'"$FP" -f "{wp}" >/dev/null 2>"{wp}.fperr"', 'RC=$?',
                      'T1=$(date +%s.%N)',
                      f'echo "###RC {s["key"]} $RC $(echo "$T1 - $T0" | bc)"']
        open(script, "w", newline="\n").write("\n".join(lines) + "\n")
        env = dict(os.environ, MSYS2_ARG_CONV_EXCL="*", MSYS_NO_PATHCONV="1")
        proc = subprocess.run(["wsl.exe", "-d", runner.WSL_DISTRO, "--", "bash",
                               runner.win_to_wsl(script)],
                              capture_output=True, text=True, env=env, timeout=7200)
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
                rec["stderr_tail"] = open(ef, encoding="utf-8",
                                          errors="replace").read()[-2000:]
            newruns.append(rec)
        print(f"  fpocket ok {sum(1 for r in newruns if r['success'])}/{len(staged)}")

        # P2Rank, one JVM per run
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
            newruns.append(dict(s, detector="p2rank", detector_version=runner.P2RANK_VERSION,
                                command=" ".join(cmd), returncode=rc, runtime_s=rt,
                                output_dir=od if good else None,
                                output_sha256=sha256_tree(od) if good else None,
                                success=bool(good),
                                stdout_tail=so[-1500:], stderr_tail=se[-1500:]))
            if i % 8 == 0 or not good:
                print(f"  p2rank {i}/{len(staged)} rc={rc} {rt}s ok={good}")
        print(f"  p2rank ok "
              f"{sum(1 for r in newruns if r['detector']=='p2rank' and r['success'])}/{len(staged)}")

        raw = os.path.join(OUT, "raw_detector")
        for r in newruns:
            if not r["success"]:
                continue
            dest = os.path.join(raw, r["detector"], r["key"])
            if os.path.exists(dest):
                shutil.rmtree(dest)
            shutil.copytree(r["output_dir"], dest)
            r["archived_output"] = os.path.relpath(dest, ROOT)

        reused_runs = [r for r in prim_runs if r["pdb_id"] in reuse_ids]
        allruns = reused_runs + newruns
        json.dump({"protocol_version": "1.3", "variant": "e9_035",
                   "n_reused_from_primary": len(reused_runs), "n_new": len(newruns),
                   "expected_new_runs": len(new_ids) * 4,
                   "runs": allruns}, open(runs_path, "w"), indent=1)

    print(f"\nNEW detector runs: expected {len(new_ids)*4}, "
          f"executed {len(newruns)}, successful {sum(1 for r in newruns if r['success'])}, "
          f"failed {sum(1 for r in newruns if not r['success'])}")

    # ---------------- classify new, merge with frozen primary rows
    pk_path = os.path.join(OUT, "pockets_classified.json")
    if not os.path.exists(pk_path):
        man = {p["pdb_id"]: p for p in new_prepared if "PREPARATION_FAILED" not in p}
        rows = []
        for run in newruns:
            if not run["success"]:
                continue
            meta = man[run["pdb_id"]]
            od = os.path.join(ROOT, run["archived_output"])
            base = run["key"]
            pockets = (cls.parse_p2rank(od, base + ".pdb") if run["detector"] == "p2rank"
                       else cls.parse_fpocket(od, base))
            cmap = meta["auth_class_map"]
            refset = set(meta["reference_site_residues"])
            bnds = meta["deletion_boundaries"]
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
                       "fusion_partner": run["partner"],
                       "target_accession": meta["target_accession"],
                       "topology": meta["topology"], "resolution": meta["resolution"],
                       "detector": run["detector"], "condition": run["condition"],
                       "rank": p["rank"], "score": p["score"],
                       "extra_score": p.get("extra_score"),
                       "centroid_x": p["centroid"][0], "centroid_y": p["centroid"][1],
                       "centroid_z": p["centroid"][2],
                       "n_residues": len(p["residues"]), "n_classifiable": n,
                       "f_target": round(ft, 4), "f_fusion": round(ff, 4),
                       "f_linker": round(fl, 4),
                       "dist_to_deletion_boundary": None if dmin is None else round(dmin, 2),
                       "ref_site_overlap_n": len(inter),
                       "ref_site_overlap_frac": None if rf is None else round(rf, 4),
                       "recovers_reference_site": bool(refset) and len(inter) >= 3
                                                  and rf >= 0.25,
                       "residues": ";".join(p["residues"])}
                for t in cls.SENSITIVITY_THRESHOLDS:
                    row[f"class_at_{t:.2f}"] = cls.classify(ft, ff, fl, n, t)
                row["pocket_class"] = row["class_at_0.70"]
                row["fusion_associated"] = row["pocket_class"] in (
                    "FUSION_DOMINATED", "INTERFACE", "LINKER")
                rows.append(row)
        reused_rows = [r for r in prim_pockets if r["pdb_id"] in reuse_ids]
        allrows = reused_rows + rows
        allrows.sort(key=lambda r: (r["detector"], r["pdb_id"], r["condition"], r["rank"]))
        json.dump({"protocol_version": "1.3", "variant": "e9_035",
                   "n_reused_pockets": len(reused_rows), "n_new_pockets": len(rows),
                   "primary_dominance": cls.PRIMARY_DOMINANCE, "pockets": allrows},
                  open(pk_path, "w"), indent=1)
        cols = [c for c in allrows[0].keys() if c != "residues"] + ["residues"]
        with open(os.path.join(OUT, "pockets_classified.tsv"), "w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in allrows:
                fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                                   for c in cols) + "\n")
        print(f"classified: {len(reused_rows)} reused + {len(rows)} new = {len(allrows)} pockets")

    d = json.load(open(pk_path))
    man = json.load(open(cache))["prepared"]
    ok = [p for p in man if "PREPARATION_FAILED" not in p]
    cov = sum(1 for p in ok if p["reference_site_available"])
    print(f"merged sensitivity set: {len(ok)} structures, {d['n_reused_pockets']} reused + "
          f"{d['n_new_pockets']} new pockets")
    print(f"reference-site coverage: {cov}/{len(ok)} = {100*cov/len(ok):.1f}%")


if __name__ == "__main__":
    main()
