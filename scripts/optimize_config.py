#!/usr/bin/env python3
"""Search for the variant-caller config that scores best on Minos.

Every candidate config is run on several practice samples and scored exactly
the way a validator scores a live round (template -> hap.py/vcfeval ->
difficulty-weighted v2 score). The objective is the MEAN score across
samples, because a single sample rewards configs that happen to suit one
region; a config only wins here if it wins on average.

What it does that the one-shot ``--practice`` mode does not:

* Bounds come from the LIVE ``/scoring/parameter-ranges`` endpoint,
  intersected with ``templates/tool_params.py``. The platform resets a value
  outside the live range to its default without telling you, so a config
  tuned outside them is not the config validators actually run.
* Bayesian search (Optuna TPE) with pruning: a candidate that is clearly
  behind after the first samples is abandoned before it costs a full pass.
* Every (config, sample) result is cached, and the study lives in SQLite, so
  an interrupted search resumes where it stopped.
* The best candidates are re-scored on held-out samples the search never
  saw, so the config you ship is not just the luckiest one on the tuning set.

Subcommands::

    # Is my config actually what validators will run?
    python scripts/optimize_config.py check --tool gatk --conf configs/gatk.conf

    # Score one or more configs on the same samples (e.g. gatk vs deepvariant)
    python scripts/optimize_config.py baseline --tool gatk --conf configs/gatk.conf

    # The search itself (resumable; re-run the same command to continue)
    python scripts/optimize_config.py search --tool gatk --trials 60

    # Re-print the leaderboard and write the best config
    python scripts/optimize_config.py report --tool gatk

Requires Docker, the Python deps from requirements.txt, ``optuna`` (only for
``search``), and network access to the platform for practice samples.
"""

from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import math
import os
import random
import secrets
import shutil
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

DEFAULT_PLATFORM_URL = "https://api.theminos.ai"
REF_BASE = os.getenv("REF_S3_BASE", "https://api.theminos.ai/reference")
SDF_FILES = ("done", "mainIndex", "nameIndex0", "namedata0", "namepointer0",
             "progress", "seqdata0", "seqpointer0", "sequenceIndex0", "summary.txt")

# Parameters that are never searched: they are either fixed by the platform,
# change the output format rather than accuracy, or are simply wrong for
# Minos BAMs (Illumina WGS, diploid).
FIXED_PARAMS = {
    # emit_ref_confidence=GVCF was measured: ~900k records, hap.py fails and
    # the round scores nothing. BP_RESOLUTION writes even more.
    "gatk": {"emit_ref_confidence": "NONE", "sample_ploidy": 2},
    "deepvariant": {"model_type": "WGS"},
    "bcftools": {"ploidy": "GRCh38", "variants_only": True,
                 "multiallelic_caller": True, "consensus_caller": False},
}

# Parameters searched by default. They are the ones that move SNP/indel
# precision and recall; the rest mostly trade runtime. ``--params all``
# searches every tunable parameter instead.
CURATED_PARAMS = {
    "gatk": [
        "pcr_indel_model",
        "standard_min_confidence_threshold_for_calling",
        "min_base_quality_score",
        "min_mapping_quality_score",
        "base_quality_score_threshold",
        "min_pruning",
        "pruning_lod_threshold",
        "recover_all_dangling_branches",
        "min_dangling_branch_length",
        "active_probability_threshold",
        "assembly_region_padding",
        "max_assembly_region_size",
        "heterozygosity",
        "indel_heterozygosity",
        "phred_scaled_global_read_mismapping_rate",
        "pair_hmm_gap_continuation_penalty",
        "max_reads_per_alignment_start",
        "contamination_fraction_to_filter",
    ],
    "deepvariant": [
        "vsc_min_fraction_indels",
        "vsc_min_fraction_snps",
        "vsc_min_count_snps",
        "vsc_min_count_indels",
        "min_mapping_quality",
        "min_base_quality",
        "realign_reads",
        "normalize_reads",
        "max_reads_per_partition",
        "qual_filter",
        "multi_allelic_qual_filter",
        "use_multiallelic_model",
    ],
    "bcftools": [
        "min_MQ", "min_BQ", "max_BQ", "delta_BQ", "adjust_MQ", "full_BAQ",
        "open_prob", "ext_prob", "gap_frac", "tandem_qual", "indel_bias",
        "min_ireads", "indels_cns", "prior", "max_depth", "max_idepth",
    ],
}

# Lower bounds the TOOL enforces but the platform does not: a value below
# them passes validation and then crashes the caller, so the round scores
# nothing. Found by running the caller, not by reading the ranges.
TOOL_MINIMUMS = {
    "gatk": {
        # GATK: "BASE_QUALITY_SCORE_THRESHOLD must be greater than or equal to 6"
        "base_quality_score_threshold": 6,
    },
}

# Searched on a log scale: their useful values span orders of magnitude.
LOG_SCALE_PARAMS = {
    "heterozygosity", "indel_heterozygosity", "active_probability_threshold",
    "adaptive_pruning_initial_error_rate", "gap_frac", "prior",
}

# Hand-picked starting points evaluated before the sampler takes over, so
# the search never ends below a known-good config.
SEED_CONFIGS = {
    "gatk": [
        {"pcr_indel_model": "NONE"},
        {"pcr_indel_model": "NONE", "recover_all_dangling_branches": True,
         "min_dangling_branch_length": 2, "assembly_region_padding": 150},
    ],
    "deepvariant": [
        {"vsc_min_fraction_indels": 0.06, "vsc_min_fraction_snps": 0.12},
    ],
    "bcftools": [],
}


# ---------------------------------------------------------------------------
# Search space
# ---------------------------------------------------------------------------

def tool_params(tool: str) -> Dict[str, Dict[str, Any]]:
    from templates import tool_params as tp
    return {
        "gatk": tp.GATK_QUALITY_PARAMS,
        "deepvariant": tp.DEEPVARIANT_QUALITY_PARAMS,
        "bcftools": tp.BCFTOOLS_QUALITY_PARAMS,
    }[tool]


def fetch_live_ranges(platform_url: str, cache_path: Path) -> Optional[Dict[str, Any]]:
    """GET /scoring/parameter-ranges; fall back to the last cached copy."""
    import httpx

    try:
        r = httpx.get(f"{platform_url.rstrip('/')}/scoring/parameter-ranges", timeout=30)
        r.raise_for_status()
        data = r.json()
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(data, indent=1))
        return data
    except Exception as e:  # noqa: BLE001 - fall back to the cache below
        if cache_path.exists():
            print(f"WARNING: live parameter ranges unavailable ({e}); using cached copy "
                  f"{cache_path}", flush=True)
            return json.loads(cache_path.read_text())
        print(f"WARNING: live parameter ranges unavailable ({e}) and no cache; using the "
              f"code ranges, which may be WIDER than what the platform accepts", flush=True)
        return None


def effective_ranges(tool: str, live: Optional[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Per-parameter definition the platform will actually honour.

    Intersection of the code definition (what the template validates) and the
    live definition (what the platform keeps rather than resetting to default).
    A parameter the live endpoint does not list is dropped: the platform would
    not forward it.
    """
    code = tool_params(tool)
    live_params = None
    if live:
        tool_live = (live.get("tools") or {}).get(tool)
        if tool_live:
            live_params = {p["name"]: p for p in tool_live.get("parameters", [])}

    out: Dict[str, Dict[str, Any]] = {}
    for name, cdef in code.items():
        d = {k: cdef.get(k) for k in ("type", "min", "max", "default", "allowed_values")}
        if live_params is not None:
            ldef = live_params.get(name)
            if ldef is None:
                continue
            if d["type"] in ("int", "float"):
                if ldef.get("min") is not None:
                    d["min"] = max(d["min"], ldef["min"])
                if ldef.get("max") is not None:
                    d["max"] = min(d["max"], ldef["max"])
                if d["type"] == "int":
                    d["min"], d["max"] = int(math.ceil(d["min"])), int(math.floor(d["max"]))
            if ldef.get("allowed_values") is not None:
                allowed = ldef["allowed_values"]
                if d["type"] == "enum":
                    d["allowed_values"] = [v for v in d["allowed_values"] if v in allowed]
                else:
                    d["allowed_values"] = list(allowed)
            if ldef.get("default") is not None:
                d["default"] = ldef["default"]
        tool_min = TOOL_MINIMUMS.get(tool, {}).get(name)
        if tool_min is not None:
            d["min"] = max(d["min"], tool_min)
        out[name] = d
    return out


def in_range(value: Any, d: Dict[str, Any]) -> bool:
    allowed = d.get("allowed_values")
    if allowed is not None and d["type"] != "enum" and d["type"] != "bool":
        return value in allowed
    if d["type"] == "bool":
        if not isinstance(value, bool):
            return False
        return allowed is None or value in allowed
    if d["type"] == "enum":
        return str(value) in [str(v) for v in allowed]
    if d["type"] == "int":
        return type(value) is int and d["min"] <= value <= d["max"]
    if d["type"] == "float":
        return (type(value) is not bool and isinstance(value, (int, float))
                and math.isfinite(value) and d["min"] <= value <= d["max"])
    return False


def is_tunable(d: Dict[str, Any]) -> bool:
    if d["type"] == "bool":
        return d.get("allowed_values") is None or len(d["allowed_values"]) > 1
    if d["type"] == "enum" or d.get("allowed_values") is not None:
        return len(d.get("allowed_values") or []) > 1
    return d["max"] > d["min"]


def clamp(value: Any, d: Dict[str, Any]) -> Any:
    """Nearest in-range value (used for defaults and seeds, never silently on user input)."""
    if in_range(value, d):
        return value
    if d["type"] == "bool" or d.get("allowed_values") is not None:
        allowed = d.get("allowed_values") or [False, True]
        return d["default"] if d["default"] in allowed else allowed[0]
    if d["type"] == "int":
        return int(min(max(int(round(value)), d["min"]), d["max"]))
    return float(min(max(float(value), d["min"]), d["max"]))


def round_float(x: float) -> float:
    return float(f"{x:.4g}")


def base_options(tool: str, ranges: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Defaults (clamped into the live ranges) plus the fixed parameters."""
    opts = {name: clamp(d["default"], d) for name, d in ranges.items()}
    for name, value in FIXED_PARAMS.get(tool, {}).items():
        if name in ranges:
            opts[name] = clamp(value, ranges[name])
    return opts


def search_params(tool: str, ranges: Dict[str, Dict[str, Any]], which: str) -> List[str]:
    names = list(ranges) if which == "all" else (
        CURATED_PARAMS[tool] if which == "curated" else [p.strip() for p in which.split(",")]
    )
    unknown = [n for n in names if n not in tool_params(tool)]
    if unknown:
        raise SystemExit(f"Unknown {tool} parameter(s): {', '.join(unknown)}")
    fixed = FIXED_PARAMS.get(tool, {})
    return [n for n in names if n in ranges and n not in fixed and is_tunable(ranges[n])]


def suggest(trial, name: str, d: Dict[str, Any]) -> Any:
    if d["type"] == "bool":
        return trial.suggest_categorical(name, d.get("allowed_values") or [False, True])
    if d["type"] == "enum" or d.get("allowed_values") is not None:
        return trial.suggest_categorical(name, list(d["allowed_values"]))
    if d["type"] == "int":
        return trial.suggest_int(name, d["min"], d["max"])
    log = name in LOG_SCALE_PARAMS and d["min"] > 0
    return round_float(trial.suggest_float(name, d["min"], d["max"], log=log))


def consistent(tool: str, opts: Dict[str, Any]) -> Dict[str, Any]:
    """Fix combinations that are individually in range but jointly invalid."""
    if tool == "gatk" and {"min_assembly_region_size", "max_assembly_region_size"} <= set(opts):
        opts["min_assembly_region_size"] = min(opts["min_assembly_region_size"],
                                               opts["max_assembly_region_size"])
    return opts


def config_hash(tool: str, options: Dict[str, Any]) -> str:
    blob = json.dumps({"tool": tool, "options": options}, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# Config files
# ---------------------------------------------------------------------------

def parse_conf(tool: str, path: Path) -> Dict[str, Any]:
    """Parse a .conf exactly as the live miner does before submitting it."""
    from utils.config_loader import extract_tool_options

    return extract_tool_options(tool, config_file=path)


def format_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return repr(round_float(v))
    return str(v)


def write_conf(path: Path, tool: str, options: Dict[str, Any], header: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {h}" for h in header] + [""]
    lines += [f"{k}={format_value(v)}" for k, v in sorted(options.items())]
    path.write_text("\n".join(lines) + "\n")


def check_options(tool: str, options: Dict[str, Any],
                  ranges: Dict[str, Dict[str, Any]]) -> Tuple[List[str], List[str]]:
    """(errors, warnings). Errors fail the round; warnings get reset by the platform."""
    from templates.tool_params import validate_and_build_flags

    errors = list(validate_and_build_flags(tool, options)["errors"])
    warnings = []
    for name, value in options.items():
        d = ranges.get(name)
        if name not in tool_params(tool):
            continue  # already an error above
        tool_min = TOOL_MINIMUMS.get(tool, {}).get(name)
        if tool_min is not None and isinstance(value, (int, float)) and value < tool_min:
            errors.append(f"{name}={value!r}: the platform accepts it but {tool} refuses "
                          f"values below {tool_min} and crashes, so the round scores nothing")
            continue
        if d is None:
            warnings.append(f"{name}={value!r}: not accepted by the platform, it is dropped")
        elif not in_range(value, d):
            bounds = (d.get("allowed_values") if d.get("allowed_values") is not None
                      else f"[{d['min']}, {d['max']}]")
            warnings.append(f"{name}={value!r}: outside the live range {bounds}; the platform "
                            f"resets it to {d['default']!r}")
    return errors, warnings


# ---------------------------------------------------------------------------
# Samples and reference data
# ---------------------------------------------------------------------------

def _import_miner():
    """Import neurons.miner without letting bittensor parse our argv."""
    saved = sys.argv
    sys.argv = [saved[0]]
    try:
        import neurons.miner as miner  # noqa: F401
    finally:
        sys.argv = saved
    return miner


def _fetch(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` atomically (the CDN refuses urllib's user agent)."""
    import httpx

    tmp = dest.with_name(dest.name + ".part")
    with httpx.stream("GET", url, follow_redirects=True, timeout=120) as r:
        r.raise_for_status()
        with tmp.open("wb") as f:
            for chunk in r.iter_bytes(1 << 20):
                f.write(chunk)
    tmp.replace(dest)


def ensure_reference(chrom: str) -> Tuple[Path, Path]:
    """Fetch the FASTA (+fai, dict) and RTG SDF for ``chrom`` if missing."""
    ref_dir = REPO_ROOT / "datasets" / "reference" / chrom
    ref_dir.mkdir(parents=True, exist_ok=True)
    fa = ref_dir / f"{chrom}.fa"
    for ext in ("fa", "fa.fai", "dict"):
        dest = ref_dir / f"{chrom}.{ext}"
        if not dest.exists() or dest.stat().st_size == 0:
            print(f"   Fetching {chrom}.{ext}...", flush=True)
            _fetch(f"{REF_BASE}/{chrom}/{chrom}.{ext}", dest)
    sdf = ref_dir / f"{chrom}.sdf"
    if not (sdf / "seqdata0").exists():
        print(f"   Fetching {chrom} RTG SDF...", flush=True)
        tmp_sdf = ref_dir / f"{chrom}.sdf.part"
        shutil.rmtree(tmp_sdf, ignore_errors=True)
        tmp_sdf.mkdir()
        for f in SDF_FILES:
            _fetch(f"{REF_BASE}/{chrom}/{chrom}.sdf/{f}", tmp_sdf / f)
        shutil.rmtree(sdf, ignore_errors=True)
        tmp_sdf.replace(sdf)
    return fa, sdf


async def _practice_client(platform_url: str):
    from bittensor_wallet import Keypair
    from utils.platform_client import MinerPlatformClient, PlatformConfig

    keypair = Keypair.create_from_uri(f"//practice-{secrets.token_hex(4)}")
    return MinerPlatformClient(
        keypair=keypair,
        config=PlatformConfig(base_url=platform_url, timeout=float(os.getenv("PLATFORM_TIMEOUT", "60"))),
        demo=False,
    )


def pick_samples(listing: List[Dict[str, Any]], n: int, chroms: Optional[List[str]],
                 seed: int) -> List[Dict[str, Any]]:
    """Stratified pick: spread across chromosomes and difficulties.

    Hard and medium samples are drawn first in each pass because they are the
    ones that separate good configs from great ones; easy samples alone tend
    to tie.
    """
    rng = random.Random(seed)
    pool = [s for s in listing if not chroms or s.get("chromosome") in chroms]
    order = {"hard": 0, "medium": 1, None: 2, "easy": 3}
    buckets: Dict[Tuple[str, Any], List[Dict[str, Any]]] = defaultdict(list)
    for s in pool:
        buckets[(s.get("chromosome"), s.get("difficulty"))].append(s)
    for b in buckets.values():
        rng.shuffle(b)
    keys = sorted(buckets, key=lambda k: (order.get(k[1], 2), rng.random()))
    picked: List[Dict[str, Any]] = []
    while len(picked) < n and any(buckets.values()):
        used_chroms = set()
        for k in keys:
            if len(picked) >= n:
                break
            if buckets[k] and k[0] not in used_chroms:
                picked.append(buckets[k].pop())
                used_chroms.add(k[0])
    return picked


def load_manifest(path: Path, holdout: int) -> List[Dict[str, Any]]:
    """Samples you already have on disk, e.g. real network rounds.

    ``path`` is a JSON list; each entry needs ``region``, ``bam``, ``truth``
    and ``mutations`` (paths relative to the manifest are fine). Optional:
    ``sample_id`` (default: the region), ``split`` ("tune"/"holdout"; by
    default the last ``holdout`` entries are held out), ``difficulty``,
    ``reference`` and ``sdf`` (default: fetched into datasets/reference), and
    ``network_top`` — the score (0-100) the network's best miner got on that
    round, which the report then compares against.
    """
    from neurons import safe_chrom

    entries = json.loads(path.read_text())
    if not isinstance(entries, list) or not entries:
        raise SystemExit(f"{path}: expected a non-empty JSON list of samples")
    base = path.parent
    samples = []
    for i, e in enumerate(entries):
        missing = [k for k in ("region", "bam", "truth", "mutations") if not e.get(k)]
        if missing:
            raise SystemExit(f"{path}: entry {i} is missing {', '.join(missing)}")
        chrom = safe_chrom(e["region"])
        if chrom is None:
            raise SystemExit(f"{path}: entry {i} has an unsupported region {e['region']!r}")
        s = dict(e)
        for k in ("bam", "truth", "mutations", "reference", "sdf"):
            if s.get(k):
                p = Path(s[k]).expanduser()
                s[k] = str(p if p.is_absolute() else (base / p).resolve())
                if not Path(s[k]).exists():
                    raise SystemExit(f"{path}: entry {i}: {k} not found: {s[k]}")
        if not (s.get("reference") and s.get("sdf")):
            s["reference"], s["sdf"] = (str(p) for p in ensure_reference(chrom))
        s.setdefault("sample_id", e["region"])
        s.setdefault("chromosome", chrom)
        s.setdefault("split", "holdout" if i >= len(entries) - holdout else "tune")
        samples.append(s)
    if not any(s["split"] == "tune" for s in samples):
        raise SystemExit(f"{path}: no tuning samples left; lower --holdout")
    miner = _import_miner()
    for s in samples:
        if not miner._ensure_bam_index(Path(s["bam"])):
            raise SystemExit(f"Could not index {s['bam']}")
    return samples


def prepare_samples(args, run_dir: Path) -> List[Dict[str, Any]]:
    """Choose (once, then reuse) and download the practice samples.

    The choice is persisted so a resumed search keeps scoring the same
    samples in the same order; otherwise the study's numbers would not be
    comparable from one session to the next.
    """
    if args.manifest:
        return load_manifest(Path(args.manifest), args.holdout)

    samples_file = run_dir / "samples.json"
    if samples_file.exists():
        chosen = json.loads(samples_file.read_text())
    else:
        async def _list():
            client = await _practice_client(args.platform_url)
            return (await client.list_practice_samples()).get("samples", [])

        listing = asyncio.run(_list())
        if not listing:
            raise SystemExit("No practice samples available on the platform.")
        chroms = args.chroms.split(",") if args.chroms else None
        total = args.n_samples + args.holdout
        chosen_meta = pick_samples(listing, total, chroms, args.seed)
        if len(chosen_meta) < total:
            raise SystemExit(f"Only {len(chosen_meta)} practice samples match; need {total}.")
        # Holdout samples are taken from the end so the tuning set is still
        # spread across chromosomes.
        chosen = []
        for i, s in enumerate(chosen_meta):
            chosen.append({
                "sample_id": s["sample_id"], "chromosome": s.get("chromosome"),
                "region": s.get("region"), "difficulty": s.get("difficulty"),
                "split": "tune" if i < args.n_samples else "holdout",
            })
        run_dir.mkdir(parents=True, exist_ok=True)
        samples_file.write_text(json.dumps(chosen, indent=1))

    miner = _import_miner()
    from neurons import safe_chrom
    from utils.path_utils import safe_round_dir_name

    async def _fetch(sample_id):
        client = await _practice_client(args.platform_url)
        return await client.get_practice_sample(sample_id)

    for s in chosen:
        chrom = safe_chrom(s["region"])
        if chrom is None:
            raise SystemExit(f"Sample {s['sample_id']} has an unsupported region {s['region']!r}")
        s["reference"], s["sdf"] = (str(p) for p in ensure_reference(chrom))
        files = None
        local_dir = REPO_ROOT / "datasets" / "practice" / safe_round_dir_name(s["sample_id"])
        if (local_dir / "input.bam").exists() and (local_dir / "truth.vcf.gz.tbi").exists():
            files = {"bam": local_dir / "input.bam", "truth": local_dir / "truth.vcf.gz",
                     "mutations": local_dir / "mutations.vcf.gz"
                     if (local_dir / "mutations.vcf.gz").exists() else None}
        if files is None:
            print(f"\n   Downloading practice sample {s['sample_id']} ({s['region']})", flush=True)
            files = miner._download_practice_files(asyncio.run(_fetch(s["sample_id"])))
            if files is None:
                raise SystemExit(f"Could not download sample {s['sample_id']}")
        if files.get("mutations") is None:
            raise SystemExit(f"Sample {s['sample_id']} has no mutations VCF; v2 cannot be "
                             f"computed without it. Drop it from {samples_file}.")
        if not miner._ensure_bam_index(Path(files["bam"])):
            raise SystemExit(f"Could not index {files['bam']}")
        s["bam"], s["truth"], s["mutations"] = str(files["bam"]), str(files["truth"]), str(files["mutations"])
    return chosen


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

# Errors that say the machine failed, not the config. A run that hits one is
# not recorded, so it cannot be mistaken for a config that scores zero.
INFRA_ERRORS = ("rc=125 ", "Docker not running", "Docker not found", "docker pull",
                "Cannot connect to the Docker daemon", "No space left on device")


def _dir_name(sample_id: str) -> str:
    """A path component for any sample id (round ids contain ':' and '+',
    which break docker -v mount specs)."""
    from utils.path_utils import safe_round_dir_name

    return safe_round_dir_name(sample_id)


class EvalCache:
    """Append-only JSONL of (config, sample) results. Safe to share across runs."""

    def __init__(self, path: Path):
        self.path = path
        self.data: Dict[Tuple[str, str], Dict[str, Any]] = {}
        if path.exists():
            for line in path.read_text().splitlines():
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue  # a line cut short by a crash
                self.data[(rec["config_hash"], rec["sample_id"])] = rec

    def get(self, h: str, sample_id: str) -> Optional[Dict[str, Any]]:
        return self.data.get((h, sample_id))

    def put(self, rec: Dict[str, Any]) -> None:
        self.data[(rec["config_hash"], rec["sample_id"])] = rec
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a") as f:
            f.write(json.dumps(rec, default=str) + "\n")


def score_sample(tool: str, options: Dict[str, Any], sample: Dict[str, Any],
                 work_dir: Path, threads: int, timeout: int, keep: bool) -> Dict[str, Any]:
    """Run one config on one sample and score it as a validator would (v2)."""
    from templates import load_template
    from utils.scoring import (AdvancedScorer, HappyScorer, difficulty_class_counts,
                               parse_happy_vcf)

    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_vcf = work_dir / "output.vcf.gz"
    rec: Dict[str, Any] = {"score": 0.0, "error": None}

    t0 = time.time()
    result = load_template(tool).variant_call(
        bam_path=Path(sample["bam"]), reference_path=Path(sample["reference"]),
        output_vcf_path=out_vcf, region=sample["region"],
        config={f"{tool}_options": options, "threads": threads, "timeout": timeout,
                "ref_build": "GRCh38"},
    )
    rec["runtime_s"] = round(time.time() - t0, 1)
    rec["variant_count"] = result.get("variant_count", 0)
    if not result.get("success") or not out_vcf.exists():
        # A validator records nothing for this miner: as good as a zero.
        rec["error"] = f"variant calling failed: {str(result.get('error'))[:500]}"
        return rec

    metrics = HappyScorer().score_vcf(
        truth_vcf=sample["truth"], query_vcf=str(out_vcf),
        reference_fasta=sample["reference"], region=sample["region"],
        reference_sdf=sample["sdf"], mutations_vcf=sample["mutations"],
    )
    if metrics is None:
        rec["error"] = "hap.py returned no metrics"
        return rec

    records = parse_happy_vcf(metrics["happy_vcf_path"], truth_vcf_path=sample["truth"])
    if not records:
        rec["error"] = "hap.py VCF could not be parsed"
        return rec
    counts = difficulty_class_counts(records)
    score = AdvancedScorer.compute_score_v2(metrics, counts)
    if score is None:
        rec["error"] = "v2 score unavailable"
        return rec

    from utils.scoring import GERMLINE_FP_SCALE, difficulty_weighted_f1, plausibility_gate

    def _f1(c):
        tp, fp, fn = c["tp"], c["fp"], c["fn"]
        return round(2 * tp / (2 * tp + fp + fn), 4) if (2 * tp + fp + fn) else None

    rec.update({
        "score": round(score, 4),
        "core": round(difficulty_weighted_f1(counts) or 0.0, 4),
        "germline": round(math.exp(-float(metrics["fp_per_target"]) / GERMLINE_FP_SCALE), 4),
        "fp_per_target": round(float(metrics["fp_per_target"]), 4),
        "gate_ok": plausibility_gate(metrics),
        "class_f1": {k: _f1(c) for k, c in counts.items()},
        "class_counts": counts,
    })
    if not keep:
        shutil.rmtree(work_dir, ignore_errors=True)
    return rec


class Evaluator:
    def __init__(self, tool: str, run_dir: Path, samples: List[Dict[str, Any]], args):
        self.tool = tool
        self.run_dir = run_dir
        self.samples = samples
        self.args = args
        self.cache = EvalCache(run_dir / "evals.jsonl")

    def evaluate(self, options: Dict[str, Any], sample: Dict[str, Any], label: str = "") -> Dict[str, Any]:
        h = config_hash(self.tool, options)
        cached = self.cache.get(h, sample["sample_id"])
        if cached is not None:
            return cached
        print(f"   [{label or h}] {sample['region']} ({sample.get('difficulty') or '-'}) ...",
              end="", flush=True)
        rec = score_sample(self.tool, options, sample,
                           self.run_dir / "work" / h / _dir_name(sample["sample_id"]),
                           self.args.threads, self.args.timeout, self.args.keep_outputs)
        infra = next((m for m in INFRA_ERRORS if m in (rec.get("error") or "")), None)
        if infra:
            # The machine failed, not the config: never cache this as a score.
            print(f" ABORTED ({rec['error'][:200]})", flush=True)
            raise SystemExit(f"Infrastructure failure ({infra}); fix it and re-run the same "
                             f"command — finished results are kept.")
        rec.update({"config_hash": h, "sample_id": sample["sample_id"], "tool": self.tool,
                    "options": options, "region": sample["region"],
                    "chromosome": sample.get("chromosome"), "split": sample.get("split"),
                    "ts": time.time()})
        self.cache.put(rec)
        tail = f" {rec['score']:.2f}" if not rec["error"] else f" FAILED ({rec['error'][:120]})"
        print(f"{tail}  [{rec.get('runtime_s', 0):.0f}s]", flush=True)
        return rec

    def evaluate_all(self, options, samples, label="") -> List[Dict[str, Any]]:
        return [self.evaluate(options, s, label) for s in samples]


def objective_value(scores: List[float], risk: float) -> float:
    """Mean score, optionally minus ``risk`` x stdev (prefers consistent configs)."""
    if not scores:
        return 0.0
    mean = statistics.fmean(scores)
    if risk and len(scores) > 1:
        mean -= risk * statistics.pstdev(scores)
    return mean


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def _setup(args) -> Tuple[Path, Dict[str, Dict[str, Any]]]:
    run_dir = Path(args.run_dir or REPO_ROOT / "optimize_runs" / args.tool).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    live = fetch_live_ranges(args.platform_url, run_dir.parent / "parameter_ranges.json")
    return run_dir, effective_ranges(args.tool, live)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _network_scores(platform_url: str, round_id: str) -> Tuple[Optional[float], Optional[float]]:
    """(top-1, rank-20) combined_final of a finalized round, on the 0-100 scale."""
    import httpx

    r = httpx.get(f"{platform_url}/scoring/rounds/{round_id}/leaderboard", timeout=60)
    r.raise_for_status()
    scores = sorted((e["combined_final"] for e in r.json().get("entries", [])
                     if isinstance(e.get("combined_final"), (int, float))), reverse=True)
    if not scores:
        return None, None
    top20 = scores[19] if len(scores) >= 20 else scores[-1]
    return round(100 * scores[0], 4), round(100 * top20, 4)


def cmd_fetch_rounds(args) -> int:
    """Download recent finalized network rounds and write a manifest for them.

    Round files are public once a round closes (/verification/round/{id}),
    and are checked against the hashes the platform committed to before the
    round. The leaderboard gives the score to beat on exactly that data.
    """
    import httpx

    url = args.platform_url.rstrip("/")
    chroms = set((args.chroms or "chr12,chr13").split(","))
    rounds = httpx.get(f"{url}/scoring/rounds", params={"limit": 100}, timeout=60).json()["rounds"]
    picked = [r for r in rounds if r.get("is_finalized") and r.get("identity_revealed")
              and r.get("region", "").split(":")[0] in chroms][:args.count]
    if not picked:
        raise SystemExit(f"No finalized rounds on {', '.join(sorted(chroms))}")

    out_dir = Path(args.out).resolve()
    entries = []
    for r in picked:
        rid = r["round_id"]
        ver = httpx.get(f"{url}/verification/round/{rid}", timeout=60).json()
        files = {f["name"]: f["url"] for f in (ver.get("reveal") or {}).get("files", [])}
        hashes = ver.get("file_hashes") or {}
        if not {"input.bam", "truth.vcf.gz", "mutations.vcf.gz"} <= set(files):
            print(f"   {rid}: round files not revealed, skipped")
            continue
        rdir = out_dir / (ver["reveal"].get("prefix", rid).rsplit("/", 1)[-1])
        rdir.mkdir(parents=True, exist_ok=True)
        print(f"   {rid}  {r['region']}", flush=True)
        for name, key in (("input.bam", "bam"), ("truth.vcf.gz", "truth_vcf"),
                          ("mutations.vcf.gz", "mutations_vcf")):
            dest = rdir / name
            if not (dest.exists() and hashes.get(key) and _sha256(dest) == hashes[key]):
                _fetch(files[name], dest)
                if hashes.get(key) and _sha256(dest) != hashes[key]:
                    dest.unlink()
                    raise SystemExit(f"{rid}: {name} does not match its committed hash")
            idx = name + (".bai" if name.endswith(".bam") else ".tbi")
            if idx in files and not (rdir / idx).exists():
                _fetch(files[idx], rdir / idx)
        top, r20 = _network_scores(url, rid)
        entries.append({
            "sample_id": rid, "region": r["region"], "chromosome": r["region"].split(":")[0],
            "bam": str(rdir / "input.bam"), "truth": str(rdir / "truth.vcf.gz"),
            "mutations": str(rdir / "mutations.vcf.gz"),
            "network_top": top, "network_r20": r20,
        })
    manifest = out_dir / "manifest.json"
    manifest.write_text(json.dumps(entries, indent=1))
    print(f"\n{len(entries)} rounds -> {manifest}")
    print(f"Next: python scripts/optimize_config.py baseline --manifest {manifest}")
    return 0


def cmd_check(args) -> int:
    _, ranges = _setup(args)
    status = 0
    for conf in args.conf:
        opts = parse_conf(args.tool, Path(conf))
        errors, warnings = check_options(args.tool, opts, ranges)
        print(f"\n{conf}:")
        for e in errors:
            print(f"   ERROR   {e}")
        for w in warnings:
            print(f"   RESET   {w}")
        if not errors and not warnings:
            print("   OK — every value is accepted as-is by the platform and the validator.")
        status |= bool(errors or warnings)
    return status


def _load_user_conf(args, ranges, conf: str) -> Optional[Dict[str, Any]]:
    """A user .conf as the platform will see it (unknown keys are refused)."""
    opts = parse_conf(args.tool, Path(conf))
    errors, warnings = check_options(args.tool, opts, ranges)
    if errors:
        print(f"   {conf}: rejected by the validator ({'; '.join(errors)})")
        return None
    for w in warnings:
        print(f"   {conf}: {w}")
    # Score what validators would actually run: reset out-of-range values.
    full = base_options(args.tool, ranges)
    for k, v in opts.items():
        if k in ranges:
            full[k] = v if in_range(v, ranges[k]) else ranges[k]["default"]
    return full


def cmd_baseline(args) -> int:
    run_dir, ranges = _setup(args)
    samples = prepare_samples(args, run_dir)
    ev = Evaluator(args.tool, run_dir, samples, args)
    confs = args.conf or [str(REPO_ROOT / "configs" / f"{args.tool}.conf")]
    rows = [("defaults", base_options(args.tool, ranges))]
    for c in confs:
        opts = _load_user_conf(args, ranges, c)
        if opts is not None:
            rows.append((Path(c).name, opts))
    print()
    for label, opts in rows:
        recs = ev.evaluate_all(opts, samples, label)
        scores = [r["score"] for r in recs]
        print(f"=> {label:30s} mean {statistics.fmean(scores):6.2f}  min {min(scores):6.2f}  "
              f"max {max(scores):6.2f}  ({len(scores)} samples)\n")
    return 0


def cmd_search(args) -> int:
    try:
        import optuna
    except ImportError:
        raise SystemExit("optuna is required for search: pip install optuna")

    run_dir, ranges = _setup(args)
    samples = prepare_samples(args, run_dir)
    tune = [s for s in samples if s["split"] == "tune"]
    ev = Evaluator(args.tool, run_dir, samples, args)
    base = base_options(args.tool, ranges)
    names = search_params(args.tool, ranges, args.params)
    print(f"\nSearching {len(names)} parameters on {len(tune)} tuning samples "
          f"(+{len(samples) - len(tune)} held out): {', '.join(names)}\n", flush=True)

    study = optuna.create_study(
        study_name=f"minos-{args.tool}",
        storage=f"sqlite:///{run_dir / 'study.db'}",
        load_if_exists=True,
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=args.seed, multivariate=True,
                                           n_startup_trials=max(8, len(names) // 2)),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=8, n_warmup_steps=1)
        if args.prune else optuna.pruners.NopPruner(),
    )

    # Seeds: defaults, the operator's current conf(s), and known-good presets.
    if len(study.trials) == 0:
        seeds = [dict(base)]
        for c in args.conf or [str(REPO_ROOT / "configs" / f"{args.tool}.conf")]:
            if Path(c).exists():
                opts = _load_user_conf(args, ranges, c)
                if opts:
                    seeds.append(opts)
        for preset in SEED_CONFIGS.get(args.tool, []):
            seeds.append({**base, **{k: clamp(v, ranges[k]) for k, v in preset.items() if k in ranges}})
        for s in seeds:
            study.enqueue_trial({k: s[k] for k in names if k in s}, skip_if_exists=True)

    def objective(trial):
        opts = dict(base)
        for name in names:
            opts[name] = suggest(trial, name, ranges[name])
        consistent(args.tool, opts)
        trial.set_user_attr("config_hash", config_hash(args.tool, opts))
        trial.set_user_attr("options", opts)
        scores = []
        for step, sample in enumerate(tune):
            rec = ev.evaluate(opts, sample, f"trial {trial.number}")
            scores.append(rec["score"])
            trial.report(statistics.fmean(scores), step)
            if trial.should_prune():
                raise optuna.TrialPruned()
        return objective_value(scores, args.risk)

    remaining = args.trials - len([t for t in study.trials if t.state.is_finished()])
    if remaining > 0:
        optuna.logging.set_verbosity(optuna.logging.WARNING)
        study.optimize(objective, n_trials=remaining, catch=(Exception,),
                       callbacks=[_progress_callback])
    return _finalize(args, run_dir, ranges, samples, ev, study)


def _progress_callback(study, trial):
    try:
        best = study.best_trial
    except ValueError:
        return
    state = trial.state.name.lower()
    val = f"{trial.value:.3f}" if trial.value is not None else "-"
    print(f"   trial {trial.number}: {state} {val} | best #{best.number} = {best.value:.3f}\n",
          flush=True)


def _finalize(args, run_dir, ranges, samples, ev, study=None) -> int:
    """Re-score the top candidates on held-out samples and write the winner."""
    tune = [s for s in samples if s["split"] == "tune"]
    holdout = [s for s in samples if s["split"] == "holdout"]

    # Candidates: every config that has a score on the full tuning set.
    by_hash: Dict[str, Dict[str, Any]] = {}
    tune_ids = {s["sample_id"] for s in tune}
    for (h, sid), rec in ev.cache.data.items():
        if rec.get("tool") != args.tool:
            continue
        by_hash.setdefault(h, {"options": rec["options"], "recs": {}})["recs"][sid] = rec
    complete = []
    for h, c in by_hash.items():
        if tune_ids <= set(c["recs"]):
            scores = [c["recs"][sid]["score"] for sid in tune_ids]
            complete.append((objective_value(scores, args.risk), h, c))
    if not complete:
        print("No config has been scored on every tuning sample yet.")
        return 1
    complete.sort(key=lambda x: -x[0])

    finalists = complete[:args.top_k]
    defaults_h = config_hash(args.tool, base_options(args.tool, ranges))
    if defaults_h in by_hash and defaults_h not in [f[1] for f in finalists]:
        finalists.append(next(x for x in complete if x[1] == defaults_h))
    if holdout:
        print(f"\nValidating the top {len(finalists)} configs on {len(holdout)} held-out samples...\n")
        for _, h, c in finalists:
            for rec in ev.evaluate_all(c["options"], holdout, f"final {h[:8]}"):
                c["recs"][rec["sample_id"]] = rec

    all_ids = [s["sample_id"] for s in samples]
    hold_ids = [s["sample_id"] for s in holdout]
    table = []
    for _, h, c in finalists:
        t_scores = [c["recs"][s]["score"] for s in tune_ids]
        h_scores = [c["recs"][s]["score"] for s in hold_ids if s in c["recs"]]
        a_scores = [c["recs"][s]["score"] for s in all_ids if s in c["recs"]]
        table.append({
            "hash": h, "options": c["options"], "tune": statistics.fmean(t_scores),
            "holdout": statistics.fmean(h_scores) if h_scores else None,
            "all": objective_value(a_scores, args.risk),
            "min": min(a_scores), "recs": c["recs"],
            "runtime": max(r.get("runtime_s", 0) for r in c["recs"].values()),
        })
    table.sort(key=lambda r: -r["all"])

    print(f"\n{'rank':<5}{'config':<18}{'tune':>8}{'holdout':>9}{'all':>8}{'worst':>8}{'max t(s)':>10}")
    for i, r in enumerate(table, 1):
        tag = " (defaults)" if r["hash"] == defaults_h else ""
        ho = f"{r['holdout']:.2f}" if r["holdout"] is not None else "-"
        print(f"{i:<5}{r['hash']:<18}{r['tune']:>8.2f}{ho:>9}{r['all']:>8.2f}"
              f"{r['min']:>8.2f}{r['runtime']:>10.0f}{tag}")

    best = table[0]
    _print_breakdown(best, samples)
    if study is not None and len(study.trials) >= 10:
        _print_importances(study)

    out = Path(args.output or REPO_ROOT / "configs" / f"{args.tool}.optimized.conf")
    ho = f"{best['holdout']:.2f}" if best["holdout"] is not None else "n/a"
    write_conf(out, args.tool, best["options"], [
        f"Generated by scripts/optimize_config.py on {time.strftime('%Y-%m-%d %H:%M')}",
        f"v2 mean score: {best['all']:.2f} over {len(samples)} practice samples "
        f"(tune {best['tune']:.2f}, held-out {ho}, worst {best['min']:.2f})",
        "Every value is inside the live /scoring/parameter-ranges bounds.",
        f"To use it: cp {out.relative_to(REPO_ROOT) if out.is_relative_to(REPO_ROOT) else out} "
        f"configs/{args.tool}.conf",
    ])
    errors, warnings = check_options(args.tool, best["options"], ranges)
    assert not errors and not warnings, (errors, warnings)
    print(f"\nBest config written to {out}")
    if best["runtime"] > 0.5 * args.timeout:
        print(f"WARNING: its slowest run took {best['runtime']:.0f}s; a slower validator may "
              f"time out. Check the runtime column before shipping it.")
    if args.per_chrom:
        _write_per_chrom(args, table, samples)
    return 0


def _print_breakdown(best, samples) -> None:
    print("\nBest config, per sample:")
    print(f"   {'region':<28}{'diff':<8}{'split':<9}{'score':>7}{'core':>7}{'germ':>7}"
          f"{'fp/tgt':>8}  indel_1bp  snp_het  indel_2_3  vs top1  vs top20")
    for s in samples:
        r = best["recs"].get(s["sample_id"])
        if r is None:
            continue
        cf = r.get("class_f1") or {}
        fmt = lambda v: f"{v:.3f}" if isinstance(v, (int, float)) else "  -  "  # noqa: E731
        print(f"   {s['region']:<28}{(s.get('difficulty') or '-'):<8}{s['split']:<9}"
              f"{r['score']:>7.2f}{r.get('core', 0):>7.3f}{r.get('germline', 0):>7.3f}"
              f"{r.get('fp_per_target', 0):>8.2f}  {fmt(cf.get('indel_1bp')):>9}  "
              f"{fmt(cf.get('snp_het')):>7}  {fmt(cf.get('indel_2_3bp')):>9}"
              + (f"  {r['score'] - float(s['network_top']):+7.2f}" if s.get("network_top") is not None else "")
              + (f"  {r['score'] - float(s['network_r20']):+8.2f}" if s.get("network_r20") is not None else "")
              + (f"   ERROR {r['error'][:60]}" if r.get("error") else ""))


def _print_importances(study) -> None:
    try:
        import optuna
        imp = optuna.importance.get_param_importances(study)
    except Exception:  # noqa: BLE001 - informational only
        return
    print("\nParameter importance (which knobs actually moved the score):")
    for name, v in list(imp.items())[:10]:
        print(f"   {name:<48}{v:6.1%}")


def _write_per_chrom(args, table, samples) -> None:
    """Per-chromosome files, only where one clearly beats the global winner.

    The miner picks configs/<tool>_<chrom>.conf over configs/<tool>.conf when
    it exists. With few samples per chromosome a per-chromosome winner is
    usually noise, so it is written only with >= 3 samples and a margin.
    """
    best = table[0]
    by_chrom: Dict[str, List[str]] = defaultdict(list)
    for s in samples:
        by_chrom[s["chromosome"]].append(s["sample_id"])
    for chrom, ids in sorted(by_chrom.items()):
        if len(ids) < 3:
            print(f"   {chrom}: {len(ids)} sample(s) — too few for a per-chromosome config")
            continue
        mean = lambda r: statistics.fmean(r["recs"][i]["score"] for i in ids)  # noqa: E731
        cand = max((r for r in table if all(i in r["recs"] for i in ids)), key=mean)
        gain = mean(cand) - mean(best)
        if cand is best or gain < args.per_chrom_margin:
            print(f"   {chrom}: global config is best (or within {args.per_chrom_margin} pts)")
            continue
        out = REPO_ROOT / "configs" / f"{args.tool}_{chrom}.conf"
        write_conf(out, args.tool, cand["options"], [
            f"Per-chromosome config for {chrom} from scripts/optimize_config.py",
            f"+{gain:.2f} pts over the global winner on {len(ids)} {chrom} samples",
        ])
        print(f"   {chrom}: wrote {out.name} (+{gain:.2f} pts on {len(ids)} samples)")


def cmd_report(args) -> int:
    run_dir, ranges = _setup(args)
    samples_file = run_dir / "samples.json"
    if not args.manifest and not samples_file.exists():
        raise SystemExit(f"No run found in {run_dir}")
    samples = prepare_samples(args, run_dir)
    ev = Evaluator(args.tool, run_dir, samples, args)
    study = None
    if (run_dir / "study.db").exists():
        try:
            import optuna
            study = optuna.load_study(study_name=f"minos-{args.tool}",
                                      storage=f"sqlite:///{run_dir / 'study.db'}")
        except Exception:  # noqa: BLE001 - the report does not need the study
            study = None
    return _finalize(args, run_dir, ranges, samples, ev, study)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, samples=True):
        sp.add_argument("--tool", choices=["gatk", "deepvariant", "bcftools"],
                        default=os.getenv("MINER_TEMPLATE", "gatk").lower())
        sp.add_argument("--platform-url", default=os.getenv("PLATFORM_URL", DEFAULT_PLATFORM_URL))
        sp.add_argument("--run-dir", help="Where results live (default optimize_runs/<tool>)")
        if not samples:
            return
        sp.add_argument("--n-samples", type=int, default=8, help="Tuning samples (default 8)")
        sp.add_argument("--holdout", type=int, default=3,
                        help="Held-out samples for the final check (default 3)")
        sp.add_argument("--chroms", help="Restrict samples, e.g. chr18,chr19")
        sp.add_argument("--manifest", help="JSON list of local samples (e.g. real network "
                                           "rounds) to use instead of practice samples")
        sp.add_argument("--seed", type=int, default=107)
        sp.add_argument("--threads", type=int, default=os.cpu_count() or 4)
        sp.add_argument("--timeout", type=int, default=1800,
                        help="Per-run variant-calling timeout in seconds")
        sp.add_argument("--keep-outputs", action="store_true",
                        help="Keep VCFs and hap.py outputs (uses disk)")
        sp.add_argument("--risk", type=float, default=0.0,
                        help="Objective = mean - risk*stdev; >0 favours consistent configs")
        sp.add_argument("--top-k", type=int, default=5, help="Finalists re-scored on holdout")
        sp.add_argument("--output", help="Where to write the best config "
                                         "(default configs/<tool>.optimized.conf)")
        sp.add_argument("--per-chrom", action="store_true",
                        help="Also write configs/<tool>_<chrom>.conf where clearly better")
        sp.add_argument("--per-chrom-margin", type=float, default=1.0)

    sp = sub.add_parser("fetch-rounds", help="Download recent real rounds + write a manifest")
    common(sp, samples=False)
    sp.add_argument("--chroms", help="Chromosomes to take rounds from (default chr12,chr13)")
    sp.add_argument("--count", type=int, default=8, help="How many recent rounds (default 8)")
    sp.add_argument("--out", default=str(REPO_ROOT / "datasets" / "rounds"))
    sp.set_defaults(func=cmd_fetch_rounds)

    sp = sub.add_parser("check", help="Check a .conf against the live platform ranges")
    common(sp, samples=False)
    sp.add_argument("--conf", nargs="+", required=True)
    sp.set_defaults(func=cmd_check)

    sp = sub.add_parser("baseline", help="Score defaults and given .conf files")
    common(sp)
    sp.add_argument("--conf", nargs="*", help="Configs to score (default configs/<tool>.conf)")
    sp.set_defaults(func=cmd_baseline)

    sp = sub.add_parser("search", help="Bayesian search for the best config (resumable)")
    common(sp)
    sp.add_argument("--trials", type=int, default=60, help="Total trials in the study")
    sp.add_argument("--params", default="curated",
                    help="'curated' (default), 'all', or a comma-separated list")
    sp.add_argument("--conf", nargs="*", help="Extra configs to seed the search with")
    sp.add_argument("--no-prune", dest="prune", action="store_false",
                    help="Score every trial on every tuning sample")
    sp.set_defaults(func=cmd_search)

    sp = sub.add_parser("report", help="Leaderboard + write the best config")
    common(sp)
    sp.set_defaults(func=cmd_report)
    return p


def main(argv=None) -> int:
    from dotenv import load_dotenv

    load_dotenv(REPO_ROOT / ".env")
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
