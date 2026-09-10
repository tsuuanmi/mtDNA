# FIS Tool

## Purpose

`src/tools/fis/` converts an R1 DataProduction workbook and an R2 MT workbook into canonical mtDNA `Sample` JSON. It is the reusable FIS tool boundary; NGS-specific comparison, identity mapping, enrichment, and set splitting remain in `src/modules/NGS/`.

## Inputs and profile policy

The tool reads:

- `T1.BasicStatistics` in R1 for `QCFail` exclusion;
- `T1.MT.Summary.QC` in R2 for nomenclature and quality metadata;
- `T2.Genotype` in R2 for raw genotype audit data;
- `T3.Seq` in R2 for required NGS sequence artifacts.

Variants are built from the R2 `Nomenclature correction genotype` field. Blank values and the FIS sentinel values `-` and `/` fall back to the raw T2 profile. Raw calls remain `FIS_Variants` metadata; they are not the canonical comparison profile.

## Core integration

The ETL constructs reference-backed `Variant` objects and `Sample(source_tool=Tool.FIS)` values. `src/core/flagging.py` is the only owner of shared FIS variant and sample flag logic. R2 quality metadata remains an independent FIS QC input.

## CLI

```bash
python -m src.tools.fis.pipeline \
  --r1-file R1.DataProduction.xlsx \
  --r2-file R2.MT.xlsx \
  --output-dir results/tools/fis/BATCH \
  --batch-id BATCH
```

The standard output is `json/statistic_fullbatch.json`, per-sample JSON, and `regions/` JSON. Region JSON writes run concurrently with 75% of the available CPU cores (minimum one worker). Add `--write-ngs-artifacts` for the established NGS root-level `FIS.tsv`, `FIS_transformed.tsv`, `FIS.json`, `FIS.seq.tsv`, `FIS_QCFail_samples.txt`, and `FASTA/` artifacts.

## Pipeline position

```text
R1/R2 Excel -> src.tools.fis -> canonical Sample JSON
                         -> src.modules.NGS.prepare_fis_comparison
                         -> core comparison -> enriched NGS report
```
