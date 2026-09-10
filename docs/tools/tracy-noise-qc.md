# Tracy Noise QC

> **Architecture direction**: This document describes the current implemented
> Tracy noise-QC behavior. The proposed modular evolution into artifact-specific
> assessors, generic trusted-coverage trimming, and sample-level QC promotion is
> documented in [tracy-qc-architecture.md](tracy-qc-architecture.md).

## Goal

The Tracy QC step finds trace areas that are probably noisy, removes variants in
`likely_noisy` ranges from final output, subtracts those ranges from the
originating trace's coverage intervals, and records each exclusion in the QC
report.

This is for Tracy only. Tracy JSON files contain the A, C, G, and T peak values
needed for this analysis.

## Current Tracy flow

The complete masking flow is:

```text
AB1 files
  -> Tracy decompose
  -> Tracy QC ranges
  -> Tracy ETL
  -> mask each trace Sample
  -> write per-sample QC report with exclusions
  -> combine forward/reverse traces
  -> final sample, region, and batch output
```

Tracy QC cannot run before decompose because the signal values are stored in the
Tracy JSON files.

## QC module and report

The Tracy package has a module similar to Sequencher's QC module:

```text
src/tools/tracy/quality_control.py
```

This module owns the QC calculations. The Tracy pipeline runs the module and
writes each report from the parent process, avoiding concurrent writes.

Each sample gets its own report:

```text
<output-dir>/preprocess/<sample-id>/qc_report.json
```

For example:

```text
results/tools/tracy/MS_260326_004/preprocess/LN_26_AA8550/qc_report.json
```

The report describes every Tracy JSON trace for that sample. It contains clean,
suspicious, likely-noisy, unavailable, and error results and remains separate
from the normal `Sample` JSON files. No batch-level `preprocess/qc_report.json`
is created.

## Metrics for each base

For each base, use the four Tracy peak values:

- **Signal**: the largest A, C, G, or T peak.
- **Noise**: the median of the other three peaks.
- **Peak SNR**: `signal / max(noise, 1.0)`, which avoids division by zero when the background channels are zero.
- **Purity**: `signal / sum of all four peaks`.
- **Background ratio**: `noise / signal`.
- **Second-peak ratio**: second-largest peak / largest peak.
- **Base quality**: Tracy's `basecallQual` value.

This is a peak-separation SNR. It is not a measurement of the complete raw
chromatogram background.

A trace is `unavailable` when required JSON fields are absent, no valid peak
measurements remain, or the trace has fewer valid measurements than the configured
window minimum. A missing or invalid measurement at one base is skipped; it does
not by itself classify the complete trace as noisy.

## Window analysis

One bad base should not hide a real variant. Therefore, the classifier uses
overlapping windows of **10 to 20 bases**.

Initial settings:

- Window size: **15 bases**
- Step size: **5 bases**
- Minimum valid bases: **10**
- Minimum supporting windows per noisy range: **2**

Example windows:

```text
1-15, 6-20, 11-25, 16-30, ...
```

Each serialized window reports:

- Median SNR
- Low-end SNR value
- Median purity
- Median background ratio
- Median second-peak ratio
- Median base quality
- Percentage of bases with low SNR
- Percentage of bases with low purity
- Percentage of bases with high background ratio
- Percentage of bases with high second-peak ratio
- Percentage of bases with low quality
- Percentage of bases with weak signal

## Initial noise rules

The default thresholds are configurable through `MTDNA_TRACY_NOISE_*` and should
be validated against representative traces:

| Measurement | Possible noise |
| --- | --- |
| SNR | below `3.0` |
| Purity | below `0.60` |
| Background ratio | above `0.35` |
| Second-peak ratio | above `0.50` |
| Base quality | below `20` |
| Signal | below `75` |

A window is **likely noisy** only when there is enough evidence:

- At least 40% of its bases have low SNR and the low-end SNR is below 2; or
- At least two different noise checks are positive, such as low SNR, low purity,
  low quality, or high background.

A window with some warning signs but not enough evidence is marked
**suspicious**. Suspicious windows are never used for masking.

## Building noisy ranges

Overlapping or directly adjacent likely-noisy windows are joined into one range
only when at least **two windows** support the range. One isolated
`likely_noisy` window remains in window details, makes the trace `suspicious`,
and does not create a range or mask variants.

```text
16120-16134
16125-16139
16130-16144
```

becomes:

```text
16120-16144
```

Each serialized range includes:

- Reference `start` and `end`
- Alignment `alignment_start` and `alignment_end`
- `supporting_windows`
- `median_snr`, `low_end_snr`, and `median_purity`
- The fixed `classification` value `likely_noisy`

The enclosing trace record supplies the filename, primer, and strand. Both
coordinate types are needed: reference coordinates drive canonical-variant
masking, while alignment coordinates retain the trace-local context for reverse
reads and indels.

## Report example

Simplified example (the real `settings` and `windows` sections contain all
configured thresholds and window metrics):

```json
{
  "tool": "tracy",
  "total": 2,
  "analyzed": 2,
  "clean": 1,
  "likely_noisy": 1,
  "suspicious": 0,
  "unavailable": 0,
  "error": 0,
  "settings": {
    "mask_enabled": true,
    "window_size": 15,
    "window_step": 5,
    "min_valid_bases": 10,
    "min_supporting_windows": 2
  },
  "traces": [
    {
      "sample_id": "2513827",
      "filename": "HV1F.json",
      "status": "likely_noisy",
      "primer": "HV1F",
      "strand": "forward",
      "ranges": [
        {
          "start": 16120,
          "end": 16144,
          "alignment_start": 80,
          "alignment_end": 104,
          "supporting_windows": 3,
          "median_snr": 2.5,
          "low_end_snr": 1.8,
          "median_purity": 0.56,
          "classification": "likely_noisy"
        }
      ],
      "excluded_variants": []
    },
    {
      "sample_id": "2513827",
      "filename": "HV1R.json",
      "status": "clean",
      "ranges": [],
      "excluded_variants": []
    }
  ]
}
```

## Variant masking behavior

The pipeline uses `likely_noisy` reference ranges to exclude variants and
per-trace coverage from final Tracy output.

The goal is:

```text
variant overlaps a likely-noisy range
  -> remove it from the per-trace Sample
  -> subtract the same reference span from that trace's intervals
  -> do not include the call or noisy interval in final output unless another clean trace covers it
  -> record it in that sample's QC report
```

This is a final-output mask. Tracy may still create an internal candidate before
the mask is applied, but the candidate must not reach the final output.

### Where the mask is applied

The current pipeline processes each trace separately and then combines forward
and reverse traces:

```text
Tracy JSON
  -> etl.process()
  -> apply noise mask to this trace's Sample
  -> collect excluded variants
  -> write this sample's QC report
  -> combine same-sample traces
  -> write final output
```

The mask must be applied before `_combine_same_tool_samples()`. This keeps trace
provenance and means a clean reverse trace can still support a variant when the
forward trace is noisy.

The existing `etl.process()` return type remains `Sample`. The pipeline applies
a Tracy-specific mask to each returned per-trace `Sample`, collects excluded
variants, recomputes variant and sample flags, and then writes the per-sample
report before combining samples. Core `Sample` and `Variant` models are
unchanged.

### Range overlap rules

Use canonical reference coordinates from the QC report:

- **SNP**: exclude when its position is inside the noisy range.
- **Insertion**: use its anchor/base position; exclude when the anchor is inside
  the range.
- **Deletion**: exclude when any reference base covered by the deletion overlaps
  the range.
- **Generated Tracy variants**: apply the same canonical-position check. Do not
  rely only on a raw alignment index because some Tracy transformations create
  variants after alignment processing.

The overlap helper must be one shared implementation. It should not be copied
into the QC module, pipeline, and ETL separately.

### Report changes

Each per-sample report keeps its existing trace entries and adds excluded
variants for the corresponding trace:

```json
{
  "filename": "LN_26_AA8550_HV1F.json",
  "status": "likely_noisy",
  "ranges": [
    {"start": 16120, "end": 16144}
  ],
  "excluded_variants": [
    {
      "pos": 16130,
      "ref": "C",
      "seq": "T",
      "reason": "overlaps likely-noisy range",
      "range": {"start": 16120, "end": 16144}
    }
  ]
}
```

Excluded variants are retained only in QC metadata. They must not be present in
normal `Sample` output, region JSONs, or batch JSONs.

### Mask configuration

The Tracy mask setting is:

```text
noise_mask_enabled = true
```

The mask is enabled by default so the normal pipeline excludes noisy-range
variants and coverage. It can be disabled through the Tracy environment
configuration when comparing masked and unmasked results:

```text
MTDNA_TRACY_NOISE_MASK_ENABLED=false
```

This is an explicit operational switch, not a compatibility path. The QC report
must still be generated when masking is disabled.

## File responsibilities and logic

### `src/config.py`

`TracySettings` defines:

```python
noise_mask_enabled: bool = True
```

The default makes the normal pipeline exclude `likely_noisy` variants and
coverage. The setting is read once and passed to per-trace processing. It never
disables QC report generation.

### `src/tools/tracy/quality_control.py`

Keep this module responsible for QC data, not `Sample` mutation or variant
calling.

Each `TraceQCResult` contains `excluded_variants`. The immutable
`with_excluded_variants()` helper lets the parent process attach mask results
without changing the QC calculations.

The module continues to own:

- Base metrics
- Window metrics
- `clean`, `suspicious`, and `likely_noisy` classification
- Merged `NoiseRange` values
- JSON report serialization

It must not import `Sample`, `Variant`, or `etl.py`.

### `src/tools/tracy/utils.py`

The canonical range-overlap helper uses primitive variant values:

```text
variant_overlaps_range(pos, ref, range_start, range_end) -> bool
```

The helper must use `pos_base()` for decimal insertion positions. It must cover
all reference bases for deletions and use the insertion anchor for insertions.
No pipeline, QC classifier, or ETL module may implement its own overlap logic.

### `src/tools/tracy/noise_mask.py`

This module owns post-ETL masking:

```text
apply_noise_mask(sample, ranges) -> NoiseMaskResult
```

Its steps are:

1. Iterate over variants from one completed trace `Sample`.
2. Check each variant against merged `likely_noisy` ranges with
   `utils.variant_overlaps_range()`.
3. Create an `ExcludedVariant` record containing position, ref, alt, reason,
   and matched range.
4. Keep only non-overlapping variants in the returned sample.
5. Subtract every inclusive noisy reference range from that trace's coverage
   intervals, splitting a span when the range is internal and removing a region
   key when no coverage remains.
6. Preserve other sample metadata.
7. Recalculate `variant_flags` with `flag_variants()`.
8. Recalculate sample flags with `recompute_sample_flags()` against the retained
   variants and reduced coverage.

This prevents a masked call from being treated as absent in coverage the final
profile no longer claims. During the same-LID merge, a clean trace can still
supply coverage for a span removed from a noisy mate.

Keeping this in `noise_mask.py` avoids adding masking responsibilities to the
already substantial ETL module. `etl.process()` still returns a `Sample` and
all existing quality, poly-C, position-specific, and region filters run before
masking. Because masking uses final canonical positions, it also covers variants
created by Tracy transformations.

### `src/tools/tracy/pipeline.py`

The two-phase pipeline works as follows:

1. `_prepare_sample` decomposes each sample and calculates QC.
2. `PreparedSample.qc_results` supplies ranges for each trace filename.
3. `_process_prepared_sample` calls `etl.process()` for one JSON file.
4. When `noise_mask_enabled` is true, it calls `apply_noise_mask()` with only
   the ranges belonging to that JSON file.
5. It returns a typed `ProcessedTrace` containing:
   - sample ID;
   - source filename;
   - masked `Sample`;
   - excluded variants.
6. The parent updates the matching `TraceQCResult` with exclusions through
   `with_excluded_variants()`.
7. The parent writes
   `preprocess/<sample-id>/qc_report.json`.
8. The parent combines only the masked `Sample` objects with
   `_combine_same_tool_samples()`.
9. Existing final sample JSON, region JSON, and batch writing continue using
   the masked samples.

If masking is disabled, the same QC report contains empty exclusion lists and
unmasked `Sample` objects continue through the pipeline.

Reports are written after per-trace masking so exclusions are included, but
before `_combine_same_tool_samples()` so source-trace information is still
available.

### Documentation and configuration files

These files describe the default-on behavior:

```text
docs/config.md       # noise_mask_enabled and environment variable
docs/tools/tracy.md  # masking position in the pipeline and report contents
CHANGELOG.md         # user-visible default masking behavior
```

### Tests

Masking coverage is split by responsibility:

```text
tests/test_tracy_noise_mask.py      # overlap and sample masking
tests/test_tracy_quality_control.py # configuration and report records
tests/test_tracy_pipeline.py        # trace matching, reports, and merging
```

No new dependency or core model change is needed.

Do not modify or add `src/core/models.py`, `src/core/sample.py`, or
`src/core/region.py`. Do not apply the mask after trace merging. Do not create a
second noise classifier or a second range-coordinate system.

## Test coverage

Shared-helper tests cover:

1. SNP outside a range returns `False`.
2. SNP inside a range returns `True`.
3. An insertion uses its anchor position.
4. A deletion overlaps when its first covered base is in the range.
5. A deletion overlaps when its last covered base is in the range.
6. Decimal positions such as `315.1` use their base anchor.

`apply_noise_mask()` tests cover:

7. A variant outside a noisy range remains in the masked sample.
8. A SNP inside a noisy range is absent from the masked sample.
9. An insertion is excluded when its anchor overlaps the range.
10. A deletion is excluded when any covered reference base overlaps the range.
11. A generated Tracy variant is excluded by canonical position.
12. Coverage intervals split around an internal noisy range and remove fully
    masked coverage even when the trace has no variants.
13. Variant flags and sample flags are recalculated after masking against the
    reduced coverage.
14. Excluded records contain the matched range and exclusion reason.

Pipeline/report tests cover:

15. Excluded variants appear only in the correct per-sample QC report.
16. Forward noisy plus reverse clean retains the clean-supported variant and
    clean-trace coverage in the merged profile.
17. Both noisy traces remove their shared noisy range from final intervals and
    produce no final variant there.
18. The final sample, region, and batch JSONs contain no masked variant or
    noisy-only coverage span.
19. Disabling `noise_mask_enabled` restores the original variant output and
    intervals while still writing QC reports with empty exclusions.
