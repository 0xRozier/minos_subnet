#!/usr/bin/env python3
"""How each GATK parameter moves the score, from every evaluation so far.

Rounds differ a lot in difficulty, so each result is first turned into a
delta against the mean of all configs on the same round; a parameter's effect
is then comparable across rounds. A random forest is fitted per target
(score, core, false positives per target) and two things are reported per
parameter:

  * permutation importance -- how much the forest's prediction degrades
    when that parameter is shuffled (share of explained variance);
  * partial dependence -- the predicted delta at low / mid / high values of
    the parameter, everything else as observed. This gives the direction.

Correlated parameters share importance, and a parameter that never varied in
the data gets none: read the table as "what the data we have says", not as
the true effect.

    python scripts/param_influence.py [run_dir ...]   # default: all of results/
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance

REPO = Path(__file__).resolve().parent.parent
TARGETS = ("score", "core", "fp_per_target")
MARGIN = 3.0  # points below the round's best config


def load(paths):
    recs = {}
    for p in paths:
        for line in Path(p).read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if r.get("error") or r.get("score") is None or r.get("core") is None:
                continue
            recs[(r["config_hash"], r["sample_id"])] = r
    return list(recs.values())


def main(argv):
    paths = argv or sorted(str(p) for p in (REPO / "results").glob("*/evals.jsonl"))
    paths += [str(p) for p in (REPO / "optimize_runs").glob("*/evals.jsonl") if not argv]
    recs = load(paths)
    # Keep viable configs only (within MARGIN of the best on their round):
    # catastrophic ones (-20 to -45 points) would otherwise drown the fine
    # effects that matter near the optimum.
    best = defaultdict(float)
    for r in recs:
        best[r["sample_id"]] = max(best[r["sample_id"]], r["score"])
    recs = [r for r in recs if r["score"] >= best[r["sample_id"]] - MARGIN]
    by_round = defaultdict(list)
    for r in recs:
        by_round[r["sample_id"]].append(r)

    params = sorted({k for r in recs for k, v in r["options"].items()
                     if isinstance(v, (int, float, bool, str))})
    enums = {k: sorted({str(r["options"].get(k)) for r in recs}) for k in params
             if any(isinstance(r["options"].get(k), str) for r in recs)}
    varying = [k for k in params
               if len({json.dumps(r["options"].get(k)) for r in recs}) > 1]

    def row(o):
        out = []
        for k in varying:
            v = o.get(k)
            if k in enums:
                out.extend(1.0 if str(v) == e else 0.0 for e in enums[k])
            else:
                out.append(float(v))
        return out

    cols = [f"{k}={e}" if k in enums else k for k in varying for e in (enums.get(k) or [None])]
    X = np.array([row(r["options"]) for r in recs])
    print(f"{len(recs)} evaluations, {len(by_round)} rounds, {len(varying)} varying parameters\n")

    for target in TARGETS:
        means = {sid: np.mean([r[target] for r in rs]) for sid, rs in by_round.items()}
        y = np.array([r[target] - means[r["sample_id"]] for r in recs])
        rf = RandomForestRegressor(n_estimators=400, min_samples_leaf=3, random_state=0, n_jobs=-1)
        rf.fit(X, y)
        imp = permutation_importance(rf, X, y, n_repeats=10, random_state=0, n_jobs=-1)
        total = max(imp.importances_mean.clip(min=0).sum(), 1e-12)
        order = np.argsort(-imp.importances_mean)
        print(f"=== {target} (delta vs round mean; forest R2 on fit = {rf.score(X, y):.2f}) ===")
        print(f"{'parameter':<52}{'share':>7}   low -> mid -> high (predicted delta)")
        for i in order[:12]:
            share = imp.importances_mean[i] / total
            if share < 0.01:
                break
            vals = np.unique(X[:, i])
            probe = [vals.min(), np.median(vals), vals.max()] if len(vals) > 2 else list(vals)
            pd = []
            for v in probe:
                Xp = X.copy()
                Xp[:, i] = v
                pd.append(rf.predict(Xp).mean())
            fmt = "f" if target != "score" else "f"
            curve = " -> ".join(f"{p:+.3{fmt}}@{v:g}" for p, v in zip(pd, probe))
            print(f"{cols[i]:<52}{share:>6.0%}   {curve}")
        print()


if __name__ == "__main__":
    main(sys.argv[1:])
