#!/usr/bin/env python3
"""Describe the false positives of a kept run, next to its true positives.

For every scored (config, round) kept with --keep-outputs, joins hap.py's
annotated VCF (TP / FP decision per record) with GATK's own VCF (DP, GQ, AD,
MQ, QD, FS, SOR, rank sums, genotype), and prints, per variant type, how FPs
and TPs are distributed over the features a GATK parameter could act on. A
feature where FPs and TPs separate cleanly points at the parameter to try.

    python scripts/fp_autopsy.py optimize_runs/diag/work/<config_hash>
"""

import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

FEATURES = ("DP", "GQ", "AB", "MQ", "QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum", "QUAL")


def gatk_records(vcf: Path):
    out = {}
    with pysam.VariantFile(str(vcf)) as f:
        for rec in f:
            s = rec.samples[0]
            ad = s.get("AD")
            dp = s.get("DP")
            alt_reads = sum(ad[1:]) if ad and None not in ad else None
            info = rec.info
            out[(rec.chrom, rec.pos)] = {
                "GT": "/".join(str(a) for a in (s.get("GT") or ())),
                "DP": dp, "GQ": s.get("GQ"),
                "AB": (alt_reads / sum(ad)) if ad and None not in ad and sum(ad) else None,
                "MQ": info.get("MQ"), "QD": info.get("QD"), "FS": info.get("FS"),
                "SOR": info.get("SOR"), "MQRankSum": info.get("MQRankSum"),
                "ReadPosRankSum": info.get("ReadPosRankSum"), "QUAL": rec.qual,
                "indel": any(len(a) != len(rec.ref) for a in rec.alts or ()),
            }
    return out


def decisions(happy_vcf: Path):
    """(chrom, pos, 'TP'|'FP') for every query-side decision."""
    with pysam.VariantFile(str(happy_vcf)) as f:
        for rec in f:
            q = rec.samples["QUERY"]
            bd = q.get("BD")
            if bd in ("TP", "FP"):
                yield rec.chrom, rec.pos, bd


def quantiles(values):
    vals = sorted(v for v in values if v is not None)
    if len(vals) < 4:
        return "n<4"
    q = statistics.quantiles(vals, n=4)
    return f"{q[0]:7.2f} {q[1]:7.2f} {q[2]:7.2f}"


def main(work: str) -> None:
    groups = defaultdict(list)  # (type, decision) -> feature dicts
    for run in sorted(Path(work).iterdir()):
        gvcf = run / "output.vcf.gz"
        hvcf = next(iter(sorted(run.glob("happy_*.vcf.gz"))), None)
        if not (gvcf.exists() and hvcf):
            continue
        calls = gatk_records(gvcf)
        for chrom, pos, bd in decisions(hvcf):
            c = calls.get((chrom, pos))
            if c is None:
                continue  # record re-represented by hap.py; rare, skipped
            groups[("INDEL" if c["indel"] else "SNP", bd)].append(c)

    for vtype in ("SNP", "INDEL"):
        fp, tp = groups[(vtype, "FP")], groups[(vtype, "TP")]
        if not fp:
            continue
        print(f"\n=== {vtype}: {len(fp)} FP vs {len(tp)} TP ===")
        print(f"{'feature':<16}{'FP  q1      q2      q3':>26}   {'TP  q1      q2      q3':>26}")
        for feat in FEATURES:
            print(f"{feat:<16}{quantiles([c[feat] for c in fp]):>26}   {quantiles([c[feat] for c in tp]):>26}")
        gts = Counter(c["GT"] for c in fp)
        print("FP genotypes:", dict(gts.most_common(4)),
              "| TP genotypes:", dict(Counter(c["GT"] for c in tp).most_common(4)))
        # How many FPs a single hard threshold would remove, against TPs lost.
        for feat, op, cuts in (("AB", "<", (0.2, 0.25, 0.3)), ("MQ", "<", (50, 55, 58)),
                               ("DP", "<", (8, 10, 12)), ("GQ", "<", (20, 30, 40)),
                               ("QD", "<", (2, 5, 10)), ("SOR", ">", (3, 4, 5)),
                               ("FS", ">", (20, 40, 60))):
            parts = []
            for cut in cuts:
                hit = (lambda v: v is not None and v < cut) if op == "<" else (lambda v: v is not None and v > cut)
                parts.append(f"{op}{cut}: -{sum(hit(c[feat]) for c in fp)} FP / -{sum(hit(c[feat]) for c in tp)} TP")
            print(f"  {feat:<4}", " | ".join(parts))


if __name__ == "__main__":
    main(sys.argv[1])
