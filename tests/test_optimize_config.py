"""Tests for scripts/optimize_config.py (the parts that need no Docker)."""

import importlib.util
from pathlib import Path

import pytest

from templates.tool_params import validate_and_build_flags

_SPEC = importlib.util.spec_from_file_location(
    "optimize_config", Path(__file__).parent.parent / "scripts" / "optimize_config.py"
)
oc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(oc)


# A trimmed copy of /scoring/parameter-ranges: narrower than the code on
# purpose, which is the case the optimizer has to get right.
LIVE = {
    "tools": {
        "gatk": {"parameters": [
            {"name": "standard_min_confidence_threshold_for_calling", "type": "float",
             "min": 30.0, "max": 100.0, "default": 30.0},
            {"name": "min_pruning", "type": "int", "min": 2, "max": 10, "default": 2},
            {"name": "pcr_indel_model", "type": "enum", "default": "CONSERVATIVE",
             "allowed_values": ["NONE", "HOSTILE", "AGGRESSIVE", "CONSERVATIVE"]},
            {"name": "sample_ploidy", "type": "int", "min": 1, "max": 10, "default": 2,
             "allowed_values": [2]},
            {"name": "dont_use_soft_clipped_bases", "type": "bool", "default": False,
             "allowed_values": [False]},
            {"name": "emit_ref_confidence", "type": "enum", "default": "NONE",
             "allowed_values": ["NONE", "GVCF", "BP_RESOLUTION"]},
        ]},
    }
}


@pytest.fixture
def gatk_ranges():
    return oc.effective_ranges("gatk", LIVE)


class TestEffectiveRanges:
    def test_live_bounds_narrow_the_code_bounds(self, gatk_ranges):
        d = gatk_ranges["standard_min_confidence_threshold_for_calling"]
        assert (d["min"], d["max"]) == (30.0, 100.0)
        assert gatk_ranges["min_pruning"]["min"] == 2

    def test_parameter_missing_from_live_is_dropped(self, gatk_ranges):
        assert "heterozygosity" not in gatk_ranges

    def test_single_allowed_value_is_not_tunable(self, gatk_ranges):
        assert not oc.is_tunable(gatk_ranges["sample_ploidy"])
        assert not oc.is_tunable(gatk_ranges["dont_use_soft_clipped_bases"])
        assert oc.is_tunable(gatk_ranges["min_pruning"])

    def test_without_live_ranges_code_ranges_are_used(self):
        ranges = oc.effective_ranges("gatk", None)
        assert ranges["standard_min_confidence_threshold_for_calling"]["min"] == 0.0
        assert "heterozygosity" in ranges


class TestCheckOptions:
    def test_out_of_live_range_is_reported_as_reset(self, gatk_ranges):
        errors, warnings = oc.check_options(
            "gatk", {"standard_min_confidence_threshold_for_calling": 20.0}, gatk_ranges)
        assert errors == []
        assert len(warnings) == 1 and "resets it to 30.0" in warnings[0]

    def test_unknown_key_is_an_error(self, gatk_ranges):
        errors, _ = oc.check_options("gatk", {"not_a_param": 1}, gatk_ranges)
        assert errors

    def test_value_that_crashes_the_tool_is_an_error(self, gatk_ranges):
        live = {"tools": {"gatk": {"parameters": LIVE["tools"]["gatk"]["parameters"] + [
            {"name": "base_quality_score_threshold", "type": "int", "min": 0, "max": 50,
             "default": 18}]}}}
        ranges = oc.effective_ranges("gatk", live)
        assert ranges["base_quality_score_threshold"]["min"] == 6
        errors, _ = oc.check_options("gatk", {"base_quality_score_threshold": 2}, ranges)
        assert errors and "crashes" in errors[0]

    def test_assembly_region_min_never_exceeds_max(self):
        opts = oc.consistent("gatk", {"min_assembly_region_size": 250,
                                      "max_assembly_region_size": 150})
        assert opts["min_assembly_region_size"] == 150

    def test_in_range_config_is_clean(self, gatk_ranges):
        assert oc.check_options("gatk", {"pcr_indel_model": "NONE", "min_pruning": 3},
                                gatk_ranges) == ([], [])


class TestSearchSpace:
    def test_base_options_are_valid_everywhere(self, gatk_ranges):
        base = oc.base_options("gatk", gatk_ranges)
        assert base["emit_ref_confidence"] == "NONE"
        assert oc.check_options("gatk", base, gatk_ranges) == ([], [])
        assert validate_and_build_flags("gatk", base)["valid"]

    def test_fixed_params_are_never_searched(self, gatk_ranges):
        names = oc.search_params("gatk", gatk_ranges, "all")
        assert "emit_ref_confidence" not in names
        assert "sample_ploidy" not in names
        assert "min_pruning" in names

    def test_unknown_param_name_is_refused(self, gatk_ranges):
        with pytest.raises(SystemExit):
            oc.search_params("gatk", gatk_ranges, "min_pruning,typo_param")

    @pytest.mark.parametrize("tool", ["gatk", "deepvariant", "bcftools"])
    def test_curated_lists_name_real_parameters(self, tool):
        known = oc.tool_params(tool)
        assert all(p in known for p in oc.CURATED_PARAMS[tool])

    def test_clamp(self, gatk_ranges):
        d = gatk_ranges["standard_min_confidence_threshold_for_calling"]
        assert oc.clamp(10.0, d) == 30.0
        assert oc.clamp(45.0, d) == 45.0
        assert oc.clamp(True, gatk_ranges["dont_use_soft_clipped_bases"]) is False

    def test_config_hash_is_order_independent(self):
        assert oc.config_hash("gatk", {"a": 1, "b": 2}) == oc.config_hash("gatk", {"b": 2, "a": 1})
        assert oc.config_hash("gatk", {"a": 1}) != oc.config_hash("deepvariant", {"a": 1})


class TestConfRoundTrip:
    def test_written_conf_parses_back_identically(self, tmp_path):
        opts = {"pcr_indel_model": "NONE", "min_pruning": 3,
                "standard_min_confidence_threshold_for_calling": 32.5,
                "heterozygosity": 0.0012345678, "recover_all_dangling_branches": True}
        path = tmp_path / "gatk.conf"
        oc.write_conf(path, "gatk", opts, ["header"])
        back = oc.parse_conf("gatk", path)
        assert back["heterozygosity"] == pytest.approx(0.001235)
        assert {k: v for k, v in back.items() if k != "heterozygosity"} == \
            {k: v for k, v in opts.items() if k != "heterozygosity"}
        assert validate_and_build_flags("gatk", back)["valid"]

    def test_integral_float_stays_a_float(self, tmp_path):
        # A float parameter written as "30" would parse back as int, which the
        # validator rejects as the wrong type and the whole round is lost.
        path = tmp_path / "gatk.conf"
        oc.write_conf(path, "gatk", {"standard_min_confidence_threshold_for_calling": 30.0}, [])
        back = oc.parse_conf("gatk", path)
        assert isinstance(back["standard_min_confidence_threshold_for_calling"], float)


class TestPickSamples:
    LISTING = (
        [{"sample_id": f"e{c}{i}", "chromosome": c, "difficulty": "easy"}
         for c in ("chr18", "chr19", "chr20") for i in range(4)]
        + [{"sample_id": f"h{c}", "chromosome": c, "difficulty": "hard"} for c in ("chr18", "chr19")]
    )

    def test_spreads_over_chromosomes_and_prefers_hard(self):
        picked = oc.pick_samples(self.LISTING, 3, None, seed=1)
        assert len({s["sample_id"] for s in picked}) == 3
        assert len({s["chromosome"] for s in picked}) == 3
        assert sum(s["difficulty"] == "hard" for s in picked) == 2

    def test_chromosome_filter(self):
        picked = oc.pick_samples(self.LISTING, 4, ["chr20"], seed=1)
        assert {s["chromosome"] for s in picked} == {"chr20"}

    def test_deterministic_for_a_seed(self):
        a = oc.pick_samples(self.LISTING, 5, None, seed=7)
        b = oc.pick_samples(self.LISTING, 5, None, seed=7)
        assert [s["sample_id"] for s in a] == [s["sample_id"] for s in b]


class TestObjectiveAndCache:
    def test_risk_penalises_spread(self):
        assert oc.objective_value([80, 80], 1.0) == 80
        assert oc.objective_value([70, 90], 1.0) == pytest.approx(70)
        assert oc.objective_value([70, 90], 0.0) == 80
        assert oc.objective_value([], 0.0) == 0.0

    def test_cache_survives_reload_and_truncated_line(self, tmp_path):
        path = tmp_path / "evals.jsonl"
        cache = oc.EvalCache(path)
        cache.put({"config_hash": "h", "sample_id": "s", "score": 81.5})
        with path.open("a") as f:
            f.write('{"config_hash": "h2", "sam')  # crash mid-write
        again = oc.EvalCache(path)
        assert again.get("h", "s")["score"] == 81.5
        assert again.get("h2", "s") is None


class TestManifest:
    def _files(self, tmp_path):
        for name in ("a.bam", "a.bam.bai", "t.vcf.gz", "m.vcf.gz", "ref.fa"):
            (tmp_path / name).write_text("x")
        (tmp_path / "ref.sdf").mkdir()

    def test_relative_paths_splits_and_defaults(self, tmp_path):
        import json
        self._files(tmp_path)
        entry = {"region": "chr13:22051001-27704000", "bam": "a.bam", "truth": "t.vcf.gz",
                 "mutations": "m.vcf.gz", "reference": "ref.fa", "sdf": "ref.sdf"}
        manifest = tmp_path / "rounds.json"
        manifest.write_text(json.dumps([dict(entry), dict(entry, sample_id="r2", network_top=85.07),
                                        dict(entry, sample_id="r3")]))
        samples = oc.load_manifest(manifest, holdout=1)
        assert [s["split"] for s in samples] == ["tune", "tune", "holdout"]
        assert samples[0]["sample_id"] == "chr13:22051001-27704000"
        assert samples[0]["chromosome"] == "chr13"
        assert samples[0]["bam"] == str(tmp_path / "a.bam")
        assert samples[1]["network_top"] == 85.07

    def test_missing_file_is_refused(self, tmp_path):
        import json
        self._files(tmp_path)
        manifest = tmp_path / "rounds.json"
        manifest.write_text(json.dumps([{"region": "chr12:1-100", "bam": "nope.bam",
                                         "truth": "t.vcf.gz", "mutations": "m.vcf.gz",
                                         "reference": "ref.fa", "sdf": "ref.sdf"}]))
        with pytest.raises(SystemExit):
            oc.load_manifest(manifest, holdout=0)


def test_sha256_matches_hashlib(tmp_path):
    import hashlib
    f = tmp_path / "x.bin"
    f.write_bytes(b"minos" * 100000)
    assert oc._sha256(f) == hashlib.sha256(f.read_bytes()).hexdigest()


def test_network_scores_takes_top1_and_rank20(monkeypatch):
    import httpx

    class _Resp:
        def raise_for_status(self):
            pass

        def json(self):
            return {"entries": [{"combined_final": 0.80 + i / 1000} for i in range(25)]
                    + [{"combined_final": None}]}

    monkeypatch.setattr(httpx, "get", lambda *a, **k: _Resp())
    top, r20 = oc._network_scores("https://x", "r")
    assert top == pytest.approx(82.4)
    assert r20 == pytest.approx(80.5)


def test_infrastructure_failure_is_not_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(oc, "score_sample", lambda *a, **k: {
        "score": 0.0, "error": "variant calling failed: Docker not running"})

    class _Args:
        threads, timeout, keep_outputs = 1, 10, False

    ev = oc.Evaluator("gatk", tmp_path, [], _Args())
    with pytest.raises(SystemExit):
        ev.evaluate({"min_pruning": 2}, {"sample_id": "s", "region": "chr12:1-2"})
    assert ev.cache.data == {}
    assert not (tmp_path / "evals.jsonl").exists()
