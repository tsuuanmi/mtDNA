# FASTA -> Sample JSON Module

> Per-HV-region consensus FASTA -> standardized `Sample` JSON for the mtDNA pipeline.

## Overview

The FASTA module (`src/fasta_to_sample.py`) is the **reverse flow** of
`src/generate_fasta.py`.  It imports an externally
supplied consensus FASTA (records `HV1`/`HV2`/`HV3`, coordinate-aligned to rCRS)
and calls variants (SNPs / insertions / deletions) into the standard `Sample`
model used by the comparison, flagging, and batch aggregation layers.

It reuses the canonical serialization (`Sample` / `Batch` /
`sample_to_dict`) and does not reimplement JSON output.  Variants are emitted
with `source_tool = Tool.FASTA`.

## Algorithm

For each present HV region:

1. **Slice** the rCRS reference to the region span (`HV1=16024-16365`,
   `HV2=73-340`, `HV3=438-576`, from `src.config`).
2. **Global pairwise alignment** of the reference slice vs the sample record
   via `Bio.Align.PairwiseAligner` (mode `global`, match `+2`, mismatch `-1`,
   gap open `-2`, extend `-1`).
3. **Walk alignment columns** -> raw calls:
   - identical bases -> match (no variant);
   - ACGT mismatch -> SNP at the 1-based reference position;
   - reference gap + query base -> insertion;
   - reference base + query gap -> deletion;
   - query `N` (truly unknown) -> **unread** (no variant emitted); a real
     IUPAC ambiguity code (Y/R/S/W/K/M/B/D/H/V) -> SNP whose `seq` is the code.
4. **Canonical normalization** (see below).
5. Build `Variant` objects and a `Sample(sample_id=<fasta stem>,
   source_tool=Tool.FASTA, intervals={...}, information={"source_tool":"fasta"})`.

Missing HV records are omitted — partial samples are allowed (no crash).

## Canonical indel normalization

The PairwiseAligner places gaps at the **leftmost** position of a tandem
repeat.  The canonical mtDNA numbering (matching `generate_fasta.py` and the
other tool pipelines) uses the **rightmost equivalent position** within the
maximal tandem-repeat run.  This module translocates each indel to the
rightmost position within the maximal homopolymer (period 1) or dinucleotide
(period 2) repeat run where sliding the gap preserves the alignment score:

- **Deletions** are placed at the rightmost L positions of the run
  (e.g. `248` -> `249` in the `AA` run; `514,515` -> `523,524` in the poly-CA
  run).  Each deleted base is emitted as a separate `Variant` with `seq="-"`.
  Because the aligner places the gap leftmost, a SNP inside the same tandem run
  is initially placed relative to that leftmost gap; when the deletion is then
  right-normalized, such a SNP is re-anchored `L` bases left (the deletion
  length) to the canonical rightmost-gap coordinate, so the variant list keeps
  round-tripping the consensus (e.g. an `A>G` in the poly-CA run with a
  concurrent `CA` deletion lands at `519`, not `521`).
- **Insertions** are right-normalized to the end of the tandem-repeat run
  containing the insertion point.  Homopolymer (period-1) insertions are slid
  per base: in the polyC region (`HV2 303-315`, two C-runs split by `T@310`) a
  `+C` in the first run (303-309) lands at `309.x` and one in the second run
  (311-315) at `315.x` (e.g. 2557844 -> `309.1` + `315.1`).  Dinucleotide
  (period-2) insertions are slid as a block to the run end and re-ordered to
  the run phase, so a `+CA` in the poly-CA run (`HV3 514-524`) lands at
  `524.1` + `524.2` (e.g. 2557921).  These conventions match the regenerate,
  sequencher and tracy pipelines and round-trip the consensus FASTA faithfully
  (the lossy all-`315` variant used by the blastn pipeline is intentionally not
  used here).  Sequential insertions at the same base are numbered
  `<base>.<i>` (e.g. `309.1`, `315.1`, `315.2`, `524.1`, `524.2`).

`validate_position()` is applied at the construction boundary: `int` for base
positions, `str` `"<n>.<i>"` for insertions.

### Nomenclature conventions (reviewer knowledge)

Most sites are normalized purely by biological run detection above.  A few
sites are a naming *convention* that no alignment can infer — equivalent
placements that all round-trip the consensus but differ in nomenclature.  The
curated/regenerate convention is encoded as an explicit, data-driven
`REGION_CONVENTIONS` table in `src/fasta_to_sample.py` (a
`RegionConvention(name, apply)` per site; each `apply` sees the region's
collected SNPs/deletions/insertion anchors *after* biological normalization and
may rewrite them to the canonical naming, including adding an insertion).  This
is the "way to apply reviewer knowledge": add an entry as data when a new
convention is agreed, with no algorithm change.

- `16189 T>C`: the singleton `T` between the `16184-16188` and `16190-16193`
  polyC runs is the named SNP `16189 T>C`; the length change is reported at the
  rightmost `C` of the merged run (e.g. `16193`), never as a `T` deletion
  (e.g. 2557851 -> `16185 C>T` + `16189 T>C` + `del 16193`).  The aligner's
  optimal placement deletes the `T` (higher score); the convention overrides
  it to the canonical naming, and the consensus still round-trips.
- **polyC run shift** (`HV2 303-315`): an adjacent `309 C>T` + `310 T>C`
  (or `308 C>T` + `del 310`) is a C-run shift, named as a deletion of the run
  end plus an insertion at `315.1` — e.g. `309 C>T` + `310 T>C` ->
  `del 309` + `ins 315.1`; `308 C>T` + `del 310` -> `del 308` + `del 309` +
  `ins 315.1`.  Both forms round-trip; the convention uses the indel form
  matching the regenerate/sequencher/tracy pipelines.
- `513 G>A` + poly-CA length change (`HV3`): a deletion of `513 G` + `514 C`
  (the unique flank base + the first run base) is named `513 G>A` plus a
  deletion of the trailing `CA` unit (`523`, `524`), not a deletion of the
  unique `G`.  The flank base is kept as a SNP; the length change moves to the
  run end.  Both forms round-trip; the convention uses the curated naming.

### Ambiguity warning

When an indel **cannot be bounded to a single tandem-repeat run** (for example
a non-repeat indel in truly external FASTA), the caller emits a `loguru`
warning and still emits a best-effort call (leftmost placement).  This is a
best-effort import heuristic, not a universal lossless guarantee; correctness
on validated/curated samples is asserted by the targeted tests.

### N / IUPAC-ambiguity fidelity

A query `N` (truly unknown) emits **no variant** ("unread"), matching the
regenerate/curated pipelines.  A real IUPAC ambiguity code
(`Y`/`R`/`S`/`W`/`K`/`M`/`B`/`D`/`H`/`V`) is a callable base and is emitted as
a SNP whose `seq` is the ambiguity code (e.g. `152 T>Y`, `16042 G>R`).

### Intervals (callable range) and N preservation

Per-region `intervals` default to the config region span (`src/config.py`)
**minus any query `N` positions** -- i.e. the *callable* range.  Because
`generate_sequence` marks positions outside `intervals` as `N`, excluding the
`N` positions here preserves them as `N` in the reconstructed HV strings
(e.g. a trailing `NNN` low-coverage tail at `574-576` -> interval `HV3`
`[438,573]`, and the HV3 string still ends `NNN`).  This both round-trips the
consensus FASTA faithfully and matches the regenerate/curated interval
convention.  Internal `N` positions split the interval into sub-intervals
(`intervals[region]` is a list of `[start, end]` spans).  `no_unread` counts
the excluded positions.

### External-FASTA robustness

Record IDs are normalized with `strip().upper()` and the query sequence is
uppercased at parse time, so a mixed-case or lowercase consensus (e.g.
`a/c/g/t`) and headers with lowercase ids or trailing descriptions are treated
as real callable bases and mapped to the correct HV region rather than being
silently dropped.

## CLI

```bash
# Process every per-sample .fasta file in an input directory (one or more files)
python -m src.fasta_to_sample \
    --input-dir results/fasta/20260615_mtDNA_444 \
    --output-dir results/tools/fasta \
    --batch-id 20260615_mtDNA_444 \
    [--ref-path ref/rCRS.fasta]
```

## Output layout

Mirrors the other tool pipelines (`results/tools/<tool>/<batch_id>/json/`):

```
results/tools/fasta/<batch_id>/json/
├── statistic_fullbatch.json
└── <sample_id>/
    └── <sample_id>.json
```

`Batch(samples).write(json_dir, batch_id, nest_batch_id=False)` produces both
the aggregate `statistic_fullbatch.json` and the per-sample JSONs.  HV
sequences, counts (`no_snps`/`no_ins`/`no_dels`/`no_identicals`/`no_unread`),
intervals, and `variant_flags` are populated by `sample_to_dict` at
aggregation time, exactly as for the other tools.

## Public API

| Symbol | Location | Description |
|--------|----------|-------------|
| `FastaToSample(ref_path="ref/rCRS.fasta")` | `src/fasta_to_sample.py` | Single class: loads ref+HV regions once; `call_variants(fasta_path) -> Sample`, `process_batch(...)`, `process_sample(...)` |
| `FastaToSample.call_variants(fasta_path)` | `src/fasta_to_sample.py` | The only variant-calling entry point (per-sample FASTA -> `Sample`) |
| `Tool.FASTA` | `src/core/models.py` | Tool provenance enum value (`"fasta"`) |
