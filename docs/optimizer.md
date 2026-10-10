# Config optimizer

`scripts/optimize_config.py` searches for the variant-caller config with the
best **average** v2 score across many practice samples, scoring each run
exactly the way a validator scores a live round (template → hap.py/vcfeval →
difficulty-weighted v2).

## Why not just tune with `--practice`?

* **One sample is noise.** Rounds rotate regions; a config that wins one
  sample often loses the next. The optimizer ranks configs by their mean over
  8+ samples, then re-checks the finalists on held-out samples it never tuned
  on.
* **The live platform bounds are narrower than the code.** The platform
  resets a value outside `GET /scoring/parameter-ranges` to its default, so a
  config tuned outside them (e.g. GATK
  `standard_min_confidence_threshold_for_calling < 30`, `min_pruning = 1`,
  `max_num_haplotypes_in_population > 128`) is **not** what validators run.
  The optimizer only searches inside the live bounds; `check` tells you if
  your current config has this problem.
* **v2 is not v1.** v2 is `0.70 × core + 0.30 × germline`, where core weights
  `indel_1bp` at 0.40, `snp_het` 0.18, `indel_2_3bp` 0.16, `indel_8bp` 0.14,
  `indel_4_7bp` 0.10 and `snp_hom` 0.02. Short indels are 66% of core, so
  indel accuracy matters far more than SNPs. Germline is
  `exp(-fp_per_target / 8)`, where `fp_per_target` counts every false
  positive in the whole region against the ~100 injected targets. On a GATK
  default run on chr21 it was 4.3 (germline 0.58), i.e. about 12 points lost
  to false positives — precision is a major lever, not a detail.

## Usage

```bash
pip install optuna            # only needed for `search`

# 1. Is my current config what validators actually run?
python scripts/optimize_config.py check --tool gatk --conf configs/gatk.conf

# 2. Compare tools on the same samples before tuning one of them
python scripts/optimize_config.py baseline --tool gatk
python scripts/optimize_config.py baseline --tool deepvariant

# 3. Search (resumable: re-run the same command to continue after a stop)
python scripts/optimize_config.py search --tool gatk --trials 60

# 4. Leaderboard + write configs/<tool>.optimized.conf
python scripts/optimize_config.py report --tool gatk
cp configs/gatk.optimized.conf configs/gatk.conf
```

Useful options:

| Option | Default | Meaning |
|---|---|---|
| `--n-samples` | 8 | Tuning samples (more = slower but less overfitting) |
| `--holdout` | 3 | Samples kept aside for the final check |
| `--chroms` | all | Restrict to e.g. `chr18,chr19` |
| `--params` | `curated` | `all`, or a comma-separated list of parameters to search |
| `--risk` | 0 | Objective `mean - risk × stdev`; `0.5` favours consistent configs |
| `--threads` | all cores | Threads per variant-caller run |
| `--no-prune` | off | Score every trial on every sample (no early stopping) |
| `--per-chrom` | off | Also write `configs/<tool>_<chrom>.conf` where one clearly wins |

Everything lives in `optimize_runs/<tool>/`: the chosen samples
(`samples.json`), every (config, sample) result (`evals.jsonl`) and the Optuna
study (`study.db`). Delete the directory to start over with new samples.

## Cost

Each trial runs the caller and hap.py once per tuning sample. With GATK on a
5 Mb region that is a few minutes per sample on 4 cores, so 60 trials × 8
samples is a long job. Pruning abandons a trial once its running mean falls
below the median of previous trials at the same point, which typically saves
a third to a half of that. Start with `--trials 30`, look at `report`, then
continue with a higher `--trials`.

## Caveats

* Practice samples cover chr18–chr22, while live rounds draw their region from
  `chromosome_rotation` in `/scoring/network-config` (currently mostly chr12
  and chr13). A config that wins on average over many chr18–22 regions is the
  best available estimate for those chromosomes; per-chromosome configs for
  chr12/chr13 cannot be tuned from practice data.
* Validator hardware differs from yours. The report shows each finalist's
  slowest runtime; a config near the timeout on your box may time out on a
  validator and score nothing.
