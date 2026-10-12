#!/usr/bin/env python3
"""Write a 2^(7-3) fractional factorial design around soft_thr65.

Seven false-positive levers, each at two levels: soft_thr65's value (-1) or
the value probed in phase 3b (+1). Sixteen configs instead of 128, resolution
IV: every main effect is estimated free of two-factor interactions, which is
what one-at-a-time probing cannot do.

Generators: A, B, C, D form a full 2^4 design; E = ABC, F = BCD, G = ACD.

For active_probability_threshold the high level is whichever of 0.01 / 0.05
scored better in the phase 3b probe (read from the evaluation cache).
"""

import itertools
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "configs" / "candidates" / "factorial"


def probe_score(name: str) -> float:
    """Mean score of configs/candidates/fp_probe/<name>.conf over its rounds."""
    sys.path.insert(0, str(REPO / "scripts"))
    sys.path.insert(0, str(REPO))
    import optimize_config as oc

    live = json.loads((REPO / "optimize_runs" / "parameter_ranges.json").read_text())
    ranges = oc.effective_ranges("gatk", live)
    opts = oc.base_options("gatk", ranges)
    opts.update(oc.parse_conf("gatk", REPO / "configs" / "candidates" / "fp_probe" / f"{name}.conf"))
    h = oc.config_hash("gatk", opts)
    scores = []
    for line in (REPO / "optimize_runs" / "real" / "evals.jsonl").read_text().splitlines():
        rec = json.loads(line)
        if rec["config_hash"] == h:
            scores.append(rec["score"])
    return sum(scores) / len(scores) if scores else float("-inf")


def main() -> None:
    base = {}
    for line in (REPO / "configs" / "candidates" / "user" / "soft_thr65.conf").read_text().splitlines():
        if line and not line.startswith("#"):
            k, _, v = line.partition("=")
            base[k] = v
    # Live values: the platform resets these two.
    base["dont_use_soft_clipped_bases"] = "false"
    base["max_num_haplotypes_in_population"] = "128"

    apt_high = "0.01" if probe_score("fp_apt0.01") >= probe_score("fp_apt0.05") else "0.05"
    factors = [  # (letter, parameter, high level)
        ("A", "max_reads_per_alignment_start", "25"),
        ("B", "active_probability_threshold", apt_high),
        ("C", "assembly_region_padding", "25"),
        ("D", "max_alternate_alleles", "2"),
        ("E", "contamination_fraction_to_filter", "0.03"),
        ("F", "indel_heterozygosity", "0.00003"),
        ("G", "recover_all_dangling_branches", "true"),
    ]

    OUT.mkdir(parents=True, exist_ok=True)
    design = []
    for i, (a, b, c, d) in enumerate(itertools.product((-1, 1), repeat=4)):
        e, f, g = a * b * c, b * c * d, a * c * d
        signs = dict(zip("ABCDEFG", (a, b, c, d, e, f, g)))
        conf = dict(base)
        for letter, param, high in factors:
            if signs[letter] > 0:
                conf[param] = high
        name = f"fac{i:02d}"
        (OUT / f"{name}.conf").write_text(
            f"# Factorial run {name}: " + " ".join(f"{k}{'+' if v > 0 else '-'}" for k, v in signs.items())
            + "\n" + "\n".join(f"{k}={v}" for k, v in sorted(conf.items())) + "\n")
        design.append({"name": name, "signs": signs})
    (OUT / "design.json").write_text(json.dumps(
        {"factors": [{"letter": l, "param": p, "high": h} for l, p, h in factors], "runs": design}, indent=1))
    print(f"16 configs written to {OUT} (apt high level = {apt_high})")


if __name__ == "__main__":
    main()
