#!/usr/bin/env bash
# Keep a long real-round job (JOB=search, default, or JOB=ablation) running
# across container restarts.
#
# Idempotent: starts the Docker daemon if it is down, then starts the search
# if it is not already running. The search resumes from optimize_runs/real
# (cached results + Optuna study), so nothing finished is recomputed.
# Prints one status line, so a scheduled check-in costs almost nothing.
set -u
cd "$(dirname "$0")/.."

# Two overlapping guards must never start two searches on the same study.
exec 9>/tmp/minos_resume_search.lock
flock 9

JOB=${JOB:-search}
mkdir -p optimize_runs/real

case "$JOB" in
  search)
    LOG=optimize_runs/real/search.log
    PATTERN="optimize_config.py search"
    DONE="Best config written"
    CMD=(python3 scripts/optimize_config.py search --tool gatk
         --manifest datasets/rounds/manifest.json --holdout 2
         --run-dir optimize_runs/real --trials "${TRIALS:-40}"
         --conf configs/candidates/user/soft_thr65.conf configs/candidates/user/soft_base.conf
         --params min_pruning,standard_min_confidence_threshold_for_calling,heterozygosity,min_dangling_branch_length,min_mapping_quality_score,min_base_quality_score,base_quality_score_threshold,pcr_indel_model,phred_scaled_global_read_mismapping_rate,pair_hmm_gap_continuation_penalty,pruning_lod_threshold)
    ;;
  ablation)
    # Phase 3: soft_thr65 with one parameter (or one group) reset to its
    # default, on 4 rounds. Shares the cache with the search.
    LOG=optimize_runs/real/ablation.log
    PATTERN="optimize_config.py baseline"
    DONE="=> def_phred_scaled"
    CMD=(python3 scripts/optimize_config.py baseline --tool gatk
         --manifest datasets/rounds/manifest_ablation.json --holdout 0
         --run-dir optimize_runs/real
         --conf configs/candidates/user/soft_thr65.conf configs/candidates/ablation/*.conf)
    ;;
  *) echo "STATUS unknown JOB=$JOB"; exit 2 ;;
esac

if ! docker info >/dev/null 2>&1; then
    sudo setsid nohup dockerd >/tmp/dockerd.log 2>&1 </dev/null 9>&- &
    for _ in $(seq 60); do docker info >/dev/null 2>&1 && break; sleep 1; done
fi
docker info >/dev/null 2>&1 || { echo "STATUS docker-failed"; exit 1; }

if pgrep -f "$PATTERN" >/dev/null; then
    state=running
elif grep -q "$DONE" "$LOG" 2>/dev/null; then
    state=finished
else
    setsid nohup "${CMD[@]}" >>"$LOG" 2>&1 </dev/null 9>&- &
    state=restarted
fi

if [ "$JOB" = ablation ]; then
    echo "STATUS $state | $(grep -c '\] chr.*\[[0-9]*s\]' "$LOG" 2>/dev/null) runs done"
    exit 0
fi

python3 - "$state" <<'EOF'
import sqlite3, sys
state = sys.argv[1]
try:
    db = sqlite3.connect("optimize_runs/real/study.db")
    rows = db.execute("SELECT state, COUNT(*) FROM trials GROUP BY state").fetchall()
    best = db.execute("SELECT MAX(v.value) FROM trial_values v JOIN trials t "
                      "ON t.trial_id = v.trial_id WHERE t.state = 'COMPLETE'").fetchone()[0]
    counts = ", ".join(f"{s.lower()}={n}" for s, n in rows)
    print(f"STATUS {state} | {counts} | best={best:.2f}" if best else f"STATUS {state} | {counts}")
except Exception as e:  # noqa: BLE001 - a status line must never fail the check
    print(f"STATUS {state} | study unreadable ({e})")
EOF
