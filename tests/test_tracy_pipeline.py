import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from src.core.batch import Batch
from src.core.models import Sample, Tool, Variant
from src.core.sample import sample_to_dict
from src.tools.tracy.etl import process
from src.tools.tracy.pipeline import (
    BatchOptions,
    PreparedSample,
    _build_combined_samples,
    _process_prepared_sample,
    _write_full_region_jsons,
    _write_sample_qc_reports,
    process_batch,
)
from src.tools.tracy.preprocessing import DecomposeConfig, decompose_sample
from src.tools.tracy.quality_control import NoiseConfig, NoiseRange, TraceQCResult


def test_decompose_returns_only_current_json_outputs(tmp_path, monkeypatch):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    (input_dir / "S1_HV1F.ab1").write_bytes(b"trace")
    stale = output_dir / "S1" / "stale_HV1R.json"
    stale.parent.mkdir(parents=True)
    stale.write_text("{}", encoding="utf-8")

    monkeypatch.setattr("src.tools.tracy.preprocessing.shutil.which", lambda _binary: "/bin/tracy")

    def fake_run(command, **_kwargs: object) -> SimpleNamespace:
        output_prefix = command[command.index("--outprefix") + 1]
        output_path = Path(f"{output_prefix}.json")
        output_path.write_text("{}", encoding="utf-8")
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr("src.tools.tracy.preprocessing.subprocess.run", fake_run)

    outputs = decompose_sample(
        "S1",
        input_dir,
        output_dir,
        "ref/rCRS.fasta",
        config=DecomposeConfig(trim=7, pratio=0.3, maxindel=1000),
    )

    assert [path.name for path in outputs] == ["S1_HV1F.json"]
    assert stale not in outputs


def test_etl_uses_reference_coordinates_across_alignment_gap(tmp_path):
    ref_align = "AAAA-AAAAAAAAAA"
    trace = {
        "ref1pos": 16024,
        "ref1align": ref_align,
        "alt1align": "A" * len(ref_align),
        "basecallPos": list(range(len(ref_align))),
        "basecallQual": [50] * len(ref_align),
        "peakA": [1000] * len(ref_align),
        "peakC": [10] * len(ref_align),
        "peakG": [10] * len(ref_align),
        "peakT": [10] * len(ref_align),
    }
    trace_path = tmp_path / "S1_HV1F.json"
    trace_path.write_text(json.dumps(trace), encoding="utf-8")

    sample = process("S1", trace_path, ref_seq="A" * 16569)

    assert len(sample.variants) == 1
    assert str(sample.variants[0].pos) == "16027.1"


def test_per_trace_masking_updates_report_before_sample_merge(tmp_path, monkeypatch):
    forward_path = tmp_path / "S1_HV1F.json"
    reverse_path = tmp_path / "S1_HV1R.json"
    noise_range = NoiseRange(
        start=100,
        end=100,
        alignment_start=0,
        alignment_end=0,
        supporting_windows=2,
        median_snr=1.5,
        low_end_snr=1.0,
        median_purity=0.4,
    )
    prepared = PreparedSample(
        sample_id="S1",
        json_paths=(forward_path, reverse_path),
        qc_results=(
            TraceQCResult("S1", forward_path.name, "likely_noisy", "HV1F", "forward", 15, ranges=(noise_range,)),
            TraceQCResult("S1", reverse_path.name, "clean", "HV1R", "reverse", 15),
        ),
    )

    def fake_process(sample_id, input_path, ref_seq) -> Sample:
        del ref_seq
        return Sample(
            sample_id=sample_id,
            source_tool=Tool.TRACY,
            variants=[Variant(pos=100, ref="A", seq="G", files=[input_path.name])],
            intervals={"HV1": [[90, 110]]},
        )

    monkeypatch.setattr("src.tools.tracy.pipeline.process", fake_process)

    processed = _process_prepared_sample(prepared, "A" * 200, noise_mask_enabled=True)

    assert processed[0].sample.variants == []
    assert processed[0].sample.intervals == {"HV1": [[90, 99], [101, 110]]}
    assert [variant.pos for variant in processed[1].sample.variants] == [100]

    combined = _build_combined_samples(processed)
    _write_sample_qc_reports(tmp_path, [prepared], processed, NoiseConfig.from_settings())

    assert [variant.pos for variant in combined["S1"].variants] == [100]
    assert combined["S1"].intervals == {"HV1": [[90, 110]]}
    report = json.loads((tmp_path / "S1" / "qc_report.json").read_text(encoding="utf-8"))
    assert report["traces"][0]["excluded_variants"][0]["pos"] == 100
    assert report["traces"][1]["excluded_variants"] == []


def test_both_noisy_traces_are_absent_from_all_final_outputs(tmp_path, monkeypatch, ref_sequence):
    paths = (tmp_path / "S1_HV1F.json", tmp_path / "S1_HV1R.json")
    noise_range = NoiseRange(16100, 16100, 0, 0, 2, 1.5, 1.0, 0.4)
    prepared = PreparedSample(
        "S1",
        paths,
        tuple(
            TraceQCResult(
                "S1",
                path.name,
                "likely_noisy",
                "HV1F" if "HV1F" in path.name else "HV1R",
                "forward" if "HV1F" in path.name else "reverse",
                15,
                ranges=(noise_range,),
            )
            for path in paths
        ),
    )

    def fake_process(sample_id, input_path, ref_seq) -> Sample:
        del ref_seq
        return Sample(
            sample_id=sample_id,
            source_tool=Tool.TRACY,
            variants=[Variant(pos=16100, ref="A", seq="G", files=[input_path.name])],
            intervals={"HV1": [[16024, 16365]]},
        )

    monkeypatch.setattr("src.tools.tracy.pipeline.process", fake_process)

    processed = _process_prepared_sample(prepared, "A" * 200, noise_mask_enabled=True)
    combined = _build_combined_samples(processed)

    assert all(trace.sample.variants == [] for trace in processed)
    assert combined["S1"].variants == []
    assert combined["S1"].intervals == {"HV1": [[16024, 16099], [16101, 16365]]}
    output = sample_to_dict(combined["S1"], ref_sequence)
    assert all(not variants for variants in output["variants"].values())
    assert output["intervals"] == {"HV1": [[16024, 16099], [16101, 16365]]}

    batch_path = Batch([combined["S1"]], ref_sequence).write(
        tmp_path / "json",
        "TEST",
        nest_batch_id=False,
    )
    batch_output = json.loads(batch_path.read_text(encoding="utf-8"))
    assert all(not variants for variants in batch_output["S1"]["variants"].values())
    assert batch_output["S1"]["intervals"] == {"HV1": [[16024, 16099], [16101, 16365]]}

    regions_dir = tmp_path / "regions"
    _write_full_region_jsons(combined["S1"], "S1", regions_dir, ref_sequence)
    region_outputs = list(regions_dir.rglob("*.json"))
    assert region_outputs
    for region_path in region_outputs:
        region = json.loads(region_path.read_text(encoding="utf-8"))
        assert all(not variants for variants in region["variants"].values())
        assert region["intervals"] == {"HV1": [[16024, 16099], [16101, 16365]]}


def test_disabling_mask_preserves_variants_and_empty_exclusions(tmp_path, monkeypatch):
    trace_path = tmp_path / "S1_HV1F.json"
    noise_range = NoiseRange(100, 100, 0, 0, 2, 1.5, 1.0, 0.4)
    prepared = PreparedSample(
        "S1",
        (trace_path,),
        (TraceQCResult("S1", trace_path.name, "likely_noisy", "HV1F", "forward", 15, ranges=(noise_range,)),),
    )

    def fake_process(sample_id, input_path, ref_seq) -> Sample:
        del ref_seq
        return Sample(
            sample_id=sample_id,
            source_tool=Tool.TRACY,
            variants=[Variant(pos=100, ref="A", seq="G", files=[input_path.name])],
            intervals={"HV1": [[90, 110]]},
        )

    monkeypatch.setattr("src.tools.tracy.pipeline.process", fake_process)

    config = replace(NoiseConfig.from_settings(), mask_enabled=False)
    processed = _process_prepared_sample(prepared, "A" * 200, noise_mask_enabled=config.mask_enabled)
    _write_sample_qc_reports(tmp_path, [prepared], processed, config)

    assert [variant.pos for variant in processed[0].sample.variants] == [100]
    assert processed[0].sample.intervals == {"HV1": [[90, 110]]}
    report = json.loads((tmp_path / "S1" / "qc_report.json").read_text(encoding="utf-8"))
    assert report["settings"]["mask_enabled"] is False
    assert report["traces"][0]["excluded_variants"] == []


def test_sample_without_trace_writes_qc_report_in_sample_folder(tmp_path, ref_sequence):
    sample_ids = tmp_path / "samples.txt"
    sample_ids.write_text("S1\n", encoding="utf-8")
    output_dir = tmp_path / "output"
    legacy_report = output_dir / "preprocess" / "qc_report.json"
    legacy_report.parent.mkdir(parents=True)
    legacy_report.write_text("{}", encoding="utf-8")

    result = process_batch(
        input_dir=tmp_path,
        output_dir=output_dir,
        options=BatchOptions(ref_path=ref_sequence, samples_path=sample_ids),
    )

    assert result == {}
    report_path = output_dir / "preprocess" / "S1" / "qc_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["tool"] == "tracy"
    assert report["total"] == 1
    assert report["analyzed"] == 0
    assert report["unavailable"] == 1
    assert report["traces"][0]["sample_id"] == "S1"
    assert not legacy_report.exists()
