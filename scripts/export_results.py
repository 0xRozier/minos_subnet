#!/usr/bin/env python3
"""Export every scored (config, round) result to results/ so it lives in git.

optimize_runs/ is git-ignored and the machine it runs on may be recycled, so
this copies the raw evaluation cache and writes two readable tables:

  results/<run>/evals.jsonl   every (config, round) record, verbatim
  results/<run>/scores.csv    one row per (config, round): score, core,
                              germline, FP/target, per-class F1, gap to the
                              network top-1 / top-20
  results/<run>/configs.csv   one row per config: name, rounds scored, mean,
                              mean gap to top-1, and its full parameters

Config names come from any .conf under configs/ whose effective options hash
to the same value; Optuna trials are named trial<N> from study.db.

    python scripts/export_results.py            # run = real
"""

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO))

import optimize_config as oc  # noqa: E402

CLASSES = ["indel_1bp", "snp_het", "indel_2_3bp", "indel_8bp", "indel_4_7bp", "snp_hom"]


def names_by_hash(ranges, run_dir: Path):
    names = {oc.config_hash("gatk", oc.base_options("gatk", ranges)): "gatk_defaults"}
    import logging

    logging.getLogger("utils.config_loader").setLevel(logging.ERROR)
    # Your own configs name a hash first, then the generated candidates.
    confs = [c for c in (REPO / "configs").rglob("*.conf")
             if c.parent.name != "configs" or c.name.startswith("gatk")]
    confs.sort(key=lambda c: (0 if c.parent.name == "user" else 1 if c.parent.name == "candidates" else 2, str(c)))
    for conf in confs:
        try:
            opts = oc.base_options("gatk", ranges)
            for k, v in oc.parse_conf("gatk", conf).items():
                if k in ranges:
                    opts[k] = v if oc.in_range(v, ranges[k]) else ranges[k]["default"]
        except Exception:  # noqa: BLE001 - not a GATK config
            continue
        names.setdefault(oc.config_hash("gatk", opts), str(conf.relative_to(REPO / "configs")))
    db = run_dir / "study.db"
    if db.exists():
        try:
            import optuna

            optuna.logging.set_verbosity(optuna.logging.WARNING)
            study = optuna.load_study(study_name="minos-gatk", storage=f"sqlite:///{db}")
            for t in study.trials:
                h = t.user_attrs.get("config_hash")
                if h:
                    names.setdefault(h, f"trial{t.number}")
        except Exception:  # noqa: BLE001 - names are a convenience
            pass
    return names


def main(run: str = "real") -> None:
    run_dir = REPO / "optimize_runs" / run
    out = REPO / "results" / run
    out.mkdir(parents=True, exist_ok=True)
    live = json.loads((REPO / "optimize_runs" / "parameter_ranges.json").read_text())
    ranges = oc.effective_ranges("gatk", live)
    names = names_by_hash(ranges, run_dir)

    manifest = {}
    for m in sorted((REPO / "datasets" / "rounds").glob("manifest*.json")):
        for e in json.loads(m.read_text()):
            manifest[e["sample_id"]] = e

    recs = {}
    for line in (run_dir / "evals.jsonl").read_text().splitlines():
        try:
            r = json.loads(line)
        except json.JSONDecodeError:
            continue
        recs[(r["config_hash"], r["sample_id"])] = r
    (out / "evals.jsonl").write_text("".join(json.dumps(r, default=str) + "\n" for r in recs.values()))

    by_conf = defaultdict(list)
    with (out / "scores.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["config", "hash", "round", "region", "score", "core", "germline", "fp_per_target",
                    *[f"f1_{c}" for c in CLASSES], "network_top1", "network_top20", "gap_top1", "error"])
        for (h, sid), r in sorted(recs.items(), key=lambda kv: (names.get(kv[0][0], kv[0][0]), kv[0][1])):
            m = manifest.get(sid, {})
            top, r20 = m.get("network_top"), m.get("network_r20")
            cf = r.get("class_f1") or {}
            gap = round(r["score"] - top, 4) if top is not None and not r.get("error") else ""
            w.writerow([names.get(h, ""), h, sid, r.get("region"), r.get("score"), r.get("core"),
                        r.get("germline"), r.get("fp_per_target"), *[cf.get(c) for c in CLASSES],
                        top, r20, gap, r.get("error") or ""])
            by_conf[h].append((r, top))

    with (out / "configs.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["config", "hash", "rounds", "mean_score", "mean_gap_top1", "options"])
        rows = []
        for h, items in by_conf.items():
            scores = [r["score"] for r, _ in items]
            gaps = [r["score"] - t for r, t in items if t is not None]
            rows.append([names.get(h, ""), h, len(items), round(sum(scores) / len(scores), 4),
                         round(sum(gaps) / len(gaps), 4) if gaps else "",
                         json.dumps(items[0][0].get("options"), sort_keys=True)])
        for row in sorted(rows, key=lambda r: (-r[2], -r[3])):
            w.writerow(row)
    print(f"{len(recs)} results, {len(by_conf)} configs -> {out}")


if __name__ == "__main__":
    main(*sys.argv[1:])
