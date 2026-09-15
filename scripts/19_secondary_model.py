#!/usr/bin/env python
"""
SECONDARY multivariable model — protocol v1.3 III.5. Kept secondary throughout.

Outcome: rank-1 pocket is fusion-associated (one row per structure x detector).
Covariates: fusion partner, topology, target size, resolution, experimental method, detector,
segment disorder. Clustering on target accession via cluster-robust covariance.

[F5] Complete / quasi-complete separation is checked BEFORE fitting. On separation the
pre-specified fallback is Firth penalised logistic regression (implemented here, since no Firth
package is installed); if that also fails the script reports that the conventional model is not
estimable rather than quietly substituting something else.

Only one interaction may be tested: detector x partner, and only if estimable.
Modelling strategy is never chosen by which yields smaller p-values.
"""
import collections, json, math, os
import numpy as np
import statsmodels.api as sm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = os.environ.get("FUSIONTAG_VARIANT", "primary")
RES = os.path.join(ROOT, "results", "confirmatory", VARIANT)
PARTNERS = ["BRIL", "T4L", "MBP"]


def firth_logit(X, y, max_iter=200, tol=1e-8):
    """Firth penalised logistic regression (Jeffreys-prior penalty), Newton-Raphson."""
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = X @ beta
        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -35, 35)))
        w = mu * (1 - mu)
        Xw = X * np.sqrt(w)[:, None]
        I = Xw.T @ Xw
        try:
            Iinv = np.linalg.pinv(I)
        except np.linalg.LinAlgError:
            return None, None
        H = (Xw @ Iinv @ Xw.T).diagonal()
        U = X.T @ (y - mu + H * (0.5 - mu))
        step = Iinv @ U
        # step halving for stability
        for _ in range(30):
            nb = beta + step
            if np.all(np.isfinite(nb)):
                break
            step /= 2
        if np.max(np.abs(nb - beta)) < tol:
            beta = nb
            break
        beta = nb
    eta = X @ beta
    mu = 1.0 / (1.0 + np.exp(-np.clip(eta, -35, 35)))
    w = mu * (1 - mu)
    Xw = X * np.sqrt(w)[:, None]
    cov = np.linalg.pinv(Xw.T @ Xw)
    return beta, np.sqrt(np.diag(cov))


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       f"confirmatory_manifest_{VARIANT}.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    pool = {r["pdb_id"]: r for r in
            json.load(open(os.path.join(MANI, "confirmatory_set_final.json")))["structures"]}
    P = json.load(open(os.path.join(RES, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    rows = []
    for det in ("p2rank", "fpocket"):
        for s in sorted(man):
            pk = idx.get((det, s, "ORIGINAL"))
            if not pk:
                continue
            m, pl = man[s], pool[s]
            rows.append({
                "pdb_id": s, "cluster": m["target_accession"], "detector": det,
                "y": 1 if pk[0]["fusion_associated"] else 0,
                "partner": m["fusion_partner"], "topology": m["topology"],
                "target_size": m["residue_class_counts"].get("TARGET", 0),
                "resolution": m["resolution"],
                "cryoem": 1 if m["method"] == "ELECTRON MICROSCOPY" else 0,
                "target_disorder": pl.get("target_disorder") or 0.0,
                "fusion_disorder": pl.get("fusion_disorder") or 0.0,
            })
    print("=" * 96)
    print(f"SECONDARY MULTIVARIABLE MODEL (variant '{VARIANT}') — kept secondary")
    print("=" * 96)
    print(f"rows = {len(rows)} (structure x detector), "
          f"clusters = {len({r['cluster'] for r in rows})}, "
          f"events = {sum(r['y'] for r in rows)}")

    # ---------------- separation check
    print("\nSEPARATION CHECK (per categorical level; a level with all-0 or all-1 outcome")
    print("produces complete separation and an infinite coefficient)")
    sep = []
    for var in ("partner", "topology", "detector"):
        for lev in sorted({r[var] for r in rows}):
            sub = [r["y"] for r in rows if r[var] == lev]
            if sub and (sum(sub) == 0 or sum(sub) == len(sub)):
                sep.append((var, lev, sum(sub), len(sub)))
            print(f"  {var:9s} {str(lev):20s} events {sum(sub):3d}/{len(sub):3d}"
                  + ("   <-- SEPARATED" if sub and (sum(sub) == 0 or sum(sub) == len(sub)) else ""))
    # detector x partner cells
    print("  detector x partner cells:")
    for det in ("p2rank", "fpocket"):
        for pt in PARTNERS:
            sub = [r["y"] for r in rows if r["detector"] == det and r["partner"] == pt]
            flag = sub and (sum(sub) == 0 or sum(sub) == len(sub))
            if flag:
                sep.append((f"{det}x{pt}", "cell", sum(sub), len(sub)))
            print(f"    {det:8s} {pt:5s} events {sum(sub):3d}/{len(sub):3d}"
                  + ("   <-- SEPARATED" if flag else ""))
    separated = bool(sep)
    print(f"\n  separation detected: {separated}" + (f"  {sep}" if separated else ""))

    # ---------------- design matrix
    def design(include_interaction):
        cols, names = [], []
        cols.append(np.ones(len(rows))); names.append("intercept")
        for pt in PARTNERS[1:]:                      # BRIL reference
            cols.append(np.array([1.0 if r["partner"] == pt else 0.0 for r in rows]))
            names.append(f"partner[{pt}]")
        topos = sorted({r["topology"] for r in rows})
        for t in topos[1:]:
            cols.append(np.array([1.0 if r["topology"] == t else 0.0 for r in rows]))
            names.append(f"topology[{t}]")
        cols.append(np.array([r["target_size"] / 100.0 for r in rows]))
        names.append("target_size/100")
        cols.append(np.array([r["resolution"] for r in rows])); names.append("resolution")
        cols.append(np.array([float(r["cryoem"]) for r in rows])); names.append("cryoEM")
        cols.append(np.array([r["target_disorder"] for r in rows]))
        names.append("target_disorder")
        cols.append(np.array([1.0 if r["detector"] == "fpocket" else 0.0 for r in rows]))
        names.append("detector[fpocket]")
        if include_interaction:
            for pt in PARTNERS[1:]:
                cols.append(np.array([1.0 if (r["detector"] == "fpocket"
                                              and r["partner"] == pt) else 0.0 for r in rows]))
                names.append(f"fpocket x {pt}")
        X = np.column_stack(cols)
        return X, names

    y = np.array([r["y"] for r in rows], dtype=float)
    clusters = np.array([r["cluster"] for r in rows])

    for interaction in (False, True):
        tag = "with detector x partner interaction" if interaction else "main effects only"
        print("\n" + "-" * 96)
        print(f"MODEL: {tag}")
        print("-" * 96)
        X, names = design(interaction)
        rank = np.linalg.matrix_rank(X)
        if rank < X.shape[1]:
            print(f"  design matrix is rank deficient ({rank} < {X.shape[1]}): NOT ESTIMABLE")
            continue
        if not separated:
            try:
                mod = sm.GLM(y, X, family=sm.families.Binomial())
                fit = mod.fit(cov_type="cluster", cov_kwds={"groups": clusters})
                print(f"  standard logistic GLM, cluster-robust SE on target accession "
                      f"({len(set(clusters))} clusters)")
                print(f"  {'term':26s} {'coef':>9s} {'SE':>8s} {'OR':>9s} "
                      f"{'95% CI':>22s} {'p':>9s}")
                for i, nm in enumerate(names):
                    b, se = fit.params[i], fit.bse[i]
                    lo, hi = b - 1.96 * se, b + 1.96 * se
                    print(f"  {nm:26s} {b:9.3f} {se:8.3f} {math.exp(b):9.3f} "
                          f"[{math.exp(lo):8.3f},{math.exp(hi):8.3f}] {fit.pvalues[i]:9.3g}")
                continue
            except Exception as exc:
                print(f"  standard GLM failed: {type(exc).__name__}: {exc}")
        print("  separation present (or GLM failed) -> pre-specified fallback: "
              "FIRTH penalised logistic regression")
        beta, se = firth_logit(X, y)
        if beta is None:
            print("  Firth did not converge: the conventional model is NOT ESTIMABLE and no "
                  "robust alternative is reported for this specification.")
            continue
        print(f"  NOTE: Firth SEs here are model-based, not cluster-robust; the partner-stratified")
        print(f"        estimates in the primary analysis carry the clustered inference.")
        print(f"  {'term':26s} {'coef':>9s} {'SE':>8s} {'OR':>9s} {'95% CI':>22s}")
        for i, nm in enumerate(names):
            b, s_ = beta[i], se[i]
            print(f"  {nm:26s} {b:9.3f} {s_:8.3f} {math.exp(np.clip(b,-30,30)):9.3f} "
                  f"[{math.exp(np.clip(b-1.96*s_,-30,30)):8.3f},"
                  f"{math.exp(np.clip(b+1.96*s_,-30,30)):8.3f}]")

    json.dump({"variant": VARIANT, "n_rows": len(rows),
               "n_clusters": len({r["cluster"] for r in rows}),
               "separation_detected": separated, "separated_cells": sep, "rows": rows},
              open(os.path.join(RES, "secondary_model.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
