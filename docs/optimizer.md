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
| `--manifest` | — | JSON list of local samples (e.g. real rounds) instead of practice ones |
| `--params` | `curated` | `all`, or a comma-separated list of parameters to search |
| `--risk` | 0 | Objective `mean - risk × stdev`; `0.5` favours consistent configs |
| `--threads` | all cores | Threads per variant-caller run |
| `--no-prune` | off | Score every trial on every sample (no early stopping) |
| `--per-chrom` | off | Also write `configs/<tool>_<chrom>.conf` where one clearly wins |

Everything lives in `optimize_runs/<tool>/`: the chosen samples
(`samples.json`), every (config, sample) result (`evals.jsonl`) and the Optuna
study (`study.db`). Delete the directory to start over with new samples.

## Tuning on real network rounds (recommended)

Round files are public once a round closes (`/verification/round/{id}`), and
every round's leaderboard is public once identities are revealed. One command
downloads the most recent finalized rounds, checks each file against the hash
the platform committed to before the round, and writes a manifest that
includes the top-1 and rank-20 scores to beat on each:

```bash
python scripts/optimize_config.py fetch-rounds --count 8          # chr12,chr13 by default
python scripts/optimize_config.py baseline --manifest datasets/rounds/manifest.json \
    --holdout 2 --run-dir optimize_runs/real --conf configs/gatk.conf
python scripts/optimize_config.py search   --manifest datasets/rounds/manifest.json \
    --holdout 2 --run-dir optimize_runs/real --trials 40
```

The per-sample report then shows `vs top1` and `vs top20` for every round.

To start from a known-good config, pass it first with `--conf`: every
parameter outside `--params` keeps its value. `--enqueue` queues targeted
variations of it, run before the sampler proposes anything (they also count
as its startup trials). Adding them to a search that is already running:

```bash
python scripts/optimize_config.py search ... --enqueue-only \
    --enqueue standard_min_confidence_threshold_for_calling=75.0 min_pruning=3
```

## Using your own rounds (`--manifest`)

Practice samples only cover chr18–chr22. If you have real network rounds on
disk (BAM + truth + mutations), describe them in a JSON list and pass it with
`--manifest`; the search then tunes on exactly the chromosomes live rounds use.

```json
[
  {"region": "chr13:22051001-27704000",
   "bam": "rounds/2026-10-09T10/input.bam",
   "truth": "rounds/2026-10-09T10/truth.vcf.gz",
   "mutations": "rounds/2026-10-09T10/mutations.vcf.gz",
   "network_top": 85.07},
  {"region": "chr12:82292001-87615000", "bam": "...", "truth": "...", "mutations": "...",
   "split": "holdout"}
]
```

Paths are relative to the manifest. `network_top` (optional, 0–100) is the
best network score on that round; the per-sample report then shows your gap
to it. Without `split`, the last `--holdout` entries are held out. Use a
separate `--run-dir` per manifest.

## Findings worth knowing

* `emit_ref_confidence=GVCF` writes ~900k records; hap.py fails and the round
  scores nothing. The optimizer keeps it at `NONE`.
* GATK refuses `base_quality_score_threshold` below 6 and crashes, although
  the platform accepts 0–50. Above ~35 it can drop every variant on
  chromosomes whose base qualities are capped.
* `dont_use_soft_clipped_bases` is fixed to `false` live: a `true` that helps
  locally is not what validators run.
* `contamination_fraction_to_filter=0.05` gained +0.2 on a default GATK run
  (no effect at 0.01–0.02).

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
  and chr13). Tune on real chr12/chr13 rounds with `--manifest` when you have
  them; otherwise a config that wins on average over many chr18–22 regions is
  the best available estimate.
* The search is a heuristic: it finds the best config among those it tries,
  validated on held-out samples. It cannot prove that no better config
  exists.
* Validator hardware differs from yours. The report shows each finalist's
  slowest runtime; a config near the timeout on your box may time out on a
  validator and score nothing.
