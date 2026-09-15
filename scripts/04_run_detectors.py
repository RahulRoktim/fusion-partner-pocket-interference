#!/usr/bin/env python
"""
Phase 2D — run P2Rank and fpocket on every (structure x condition).

24 structures x 2 conditions x 2 detectors = 96 runs expected.

No per-structure tuning. Every run records: full command, return code, runtime, detector
version, input SHA-256 and output SHA-256. Technical failures are recorded, never silently
replaced or retried with different settings.

Execution note: detectors run in a space-free scratch workspace because the repository path
contains spaces and fpocket is invoked through WSL. Raw outputs are copied back into
results/raw_detector/ afterwards.
"""
import json, os, sys, shutil, subprocess, hashlib, time, importlib.util

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    "runner", os.path.join(ROOT, "scripts", "common_runner.py"))
runner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(runner)
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")
RAWOUT = os.path.join(RESULTS, "raw_detector")
LOGS = os.path.join(ROOT, "logs")
WORK = os.environ.get(
    "FUSIONTAG_WORK",
    r"C:\Users\user\AppData\Local\Temp\claude"
    r"\C--AI-PROJECTS-Softwares-Fusion-Tag-Hazard"
    r"\5c4a9bb4-b9dd-40f9-b5cc-b8b58139430b\scratchpad\runs")

FPOCKET_WSL = runner.FPOCKET_WSL
P2RANK_VERSION = runner.P2RANK_VERSION
FPOCKET_VERSION = runner.FPOCKET_VERSION

for d in (RAWOUT, LOGS, WORK):
    os.makedirs(d, exist_ok=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def sha256_tree(path):
    """Stable hash over a directory's file contents."""
    h = hashlib.sha256()
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for f in sorted(files):
            fp = os.path.join(root, f)
            h.update(os.path.relpath(fp, path).replace("\\", "/").encode())
            h.update(sha256_file(fp).encode())
    return h.hexdigest()


win_to_wsl = runner.win_to_wsl


def run_p2rank(inp, outdir):
    cmd = runner.p2rank_cmd(inp, outdir)
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env=runner.p2rank_env(), timeout=1800)
    return {"command": " ".join(cmd), "returncode": proc.returncode,
            "runtime_s": round(time.time() - t0, 2),
            "stdout": proc.stdout[-8000:], "stderr": proc.stderr[-8000:]}


def run_fpocket_batch(jobs):
    """Run all fpocket jobs in one WSL invocation; returns {job_key: record}."""
    script = os.path.join(WORK, "_fpocket_batch.sh")
    lines = ["#!/bin/bash", "set -u", 'FP="%s"' % FPOCKET_WSL, ""]
    for key, inp in jobs:
        wp = win_to_wsl(inp)
        lines += [
            f'echo "###JOB {key}"',
            'T0=$(date +%s.%N)',
            f'"$FP" -f "{wp}" > /dev/null 2>"{wp}.fperr"',
            'RC=$?',
            'T1=$(date +%s.%N)',
            f'echo "###RC {key} $RC $(echo "$T1 - $T0" | bc)"',
        ]
    with open(script, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    env = dict(os.environ, MSYS2_ARG_CONV_EXCL="*", MSYS_NO_PATHCONV="1")
    proc = subprocess.run(["wsl.exe", "-d", "Ubuntu", "--", "bash", win_to_wsl(script)],
                          capture_output=True, text=True, env=env, timeout=7200)
    out = {}
    for line in proc.stdout.replace("\r", "").splitlines():
        if line.startswith("###RC "):
            _, key, rc, rt = line.split()
            out[key] = {"returncode": int(rc), "runtime_s": round(float(rt), 2)}
    return out, proc


def main():
    man = json.load(open(os.path.join(MANI, "pilot_manifest.json")))
    prepared = [p for p in man["prepared"] if "PREPARATION_FAILED" not in p]
    print(f"structures to run: {len(prepared)}")

    # --- stage inputs into the space-free workspace
    staged = []
    for p in prepared:
        for cond, key in (("ORIGINAL", "original_pdb"), ("REMOVED", "removed_pdb")):
            src = os.path.join(ROOT, p[key])
            base = f"{p['pdb_id']}_{p['auth_chain']}_{cond}"
            dstdir = os.path.join(WORK, base)
            os.makedirs(dstdir, exist_ok=True)
            dst = os.path.join(dstdir, base + ".pdb")
            shutil.copyfile(src, dst)
            staged.append({"pdb_id": p["pdb_id"], "chain": p["auth_chain"],
                           "partner": p["fusion_partner"], "condition": cond,
                           "key": base, "input": dst,
                           "input_sha256": sha256_file(dst)})

    runs = []

    # ---------------- fpocket (single WSL batch) ----------------
    print(f"\nfpocket: {len(staged)} runs (one WSL batch)")
    jobs = [(s["key"], s["input"]) for s in staged]
    fp_res, fp_proc = run_fpocket_batch(jobs)
    for s in staged:
        r = fp_res.get(s["key"], {"returncode": -999, "runtime_s": None})
        outdir = os.path.join(os.path.dirname(s["input"]), s["key"] + "_out")
        ok = os.path.isdir(outdir) and r["returncode"] == 0
        rec = dict(s, detector="fpocket", detector_version=FPOCKET_VERSION,
                   command=f'{FPOCKET_WSL} -f {win_to_wsl(s["input"])}',
                   returncode=r["returncode"], runtime_s=r["runtime_s"],
                   output_dir=outdir if ok else None,
                   output_sha256=sha256_tree(outdir) if ok else None,
                   success=bool(ok))
        errf = s["input"] + ".fperr"
        if os.path.exists(errf):
            rec["stderr"] = open(errf, encoding="utf-8", errors="replace").read()[-4000:]
        runs.append(rec)
    nok = sum(1 for r in runs if r["detector"] == "fpocket" and r["success"])
    print(f"  fpocket succeeded: {nok}/{len(staged)}")
    with open(os.path.join(LOGS, "fpocket_batch.log"), "w", encoding="utf-8") as fh:
        fh.write(fp_proc.stdout + "\n===STDERR===\n" + fp_proc.stderr)

    # ---------------- P2Rank (one JVM per run) ----------------
    print(f"\nP2Rank: {len(staged)} runs")
    for i, s in enumerate(staged, 1):
        outdir = os.path.join(os.path.dirname(s["input"]), "p2rank_out")
        try:
            r = run_p2rank(s["input"], outdir)
        except subprocess.TimeoutExpired:
            r = {"command": "timeout", "returncode": -998, "runtime_s": None,
                 "stdout": "", "stderr": "TIMEOUT"}
        pred = os.path.join(outdir, os.path.basename(s["input"]) + "_predictions.csv")
        ok = r["returncode"] == 0 and os.path.exists(pred)
        runs.append(dict(s, detector="p2rank", detector_version=P2RANK_VERSION,
                         command=r["command"], returncode=r["returncode"],
                         runtime_s=r["runtime_s"],
                         output_dir=outdir if ok else None,
                         output_sha256=sha256_tree(outdir) if ok else None,
                         success=bool(ok),
                         stdout_tail=r["stdout"][-2000:], stderr_tail=r["stderr"][-2000:]))
        if i % 8 == 0 or not ok:
            print(f"  [{i}/{len(staged)}] {s['key']}  rc={r['returncode']} "
                  f"{r['runtime_s']}s ok={ok}")
    nok = sum(1 for r in runs if r["detector"] == "p2rank" and r["success"])
    print(f"  P2Rank succeeded: {nok}/{len(staged)}")

    # ---------------- persist raw outputs into the repository ----------------
    print("\ncopying raw detector outputs into results/raw_detector/ ...")
    for r in runs:
        if not r["success"]:
            continue
        dest = os.path.join(RAWOUT, r["detector"], r["key"])
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(r["output_dir"], dest)
        r["archived_output"] = os.path.relpath(dest, ROOT)

    json.dump({"protocol_version": "1.1",
               "expected_runs": len(prepared) * 2 * 2,
               "actual_runs": len(runs),
               "runs": runs},
              open(os.path.join(RESULTS, "detector_runs.json"), "w"), indent=1)

    print(f"\nTOTAL RUNS: {len(runs)} (expected {len(prepared)*2*2})")
    for det in ("p2rank", "fpocket"):
        sub = [r for r in runs if r["detector"] == det]
        print(f"  {det}: {sum(1 for r in sub if r['success'])}/{len(sub)} succeeded")
    fails = [r for r in runs if not r["success"]]
    if fails:
        print("\nFAILURES (recorded, not replaced):")
        for f in fails:
            print(f"  {f['detector']:8s} {f['key']:28s} rc={f['returncode']}")


if __name__ == "__main__":
    main()
