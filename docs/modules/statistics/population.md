# mtDNA population analysis

Analyze the **current upstream-approved export**, not raw reads or QC flags:

```bash
uv run python -m src.modules.statistics.population \
  --input /home/bca/mtDNA/science/results/merged/merged_statistics.json \
  --output-dir /home/bca/mtDNA/science/results/modules/statistics/population/test
```

The input path above is the CLI default. Without `--output-dir`, reports land in
`<results>/modules/statistics/population/current-cohort-<timestamp>` under the
configured results directory. The output directory must be new or empty; use a
separate directory for each execution. `--precision` controls executive-text
rounding only. Wilson intervals default to 95%; use `--confidence-level 0.99` or
`--no-confidence-intervals`.

## Eligibility and input contract

The operator supplies an upstream-approved JSON object keyed by the repository's
string sample IDs, in the `merged_statistics.json` format. Every record is eligible.
The file is reloaded once per execution and its exact bytes are SHA-256 hashed.
Decimal JSON position tokens are decoded losslessly as text before the shared
position normalizer runs: numeric `309.10` remains insertion index 10, not 1.
Python floating-point positions and boolean coordinates are not accepted; interval
endpoints must be integer JSON values.
The producer must refresh it to contain all current PASS samples, one biological
sample per ID. This module cannot certify completeness/freshness of that export.
The existing merger's file availability and last-batch overwrite behavior are **not**
independent proof of PASS or biological deduplication. Approval belongs upstream.

Entries must declare canonical `variants` (flat, type-grouped or HV-grouped) and
per-region callable `intervals`. The finalized `HV1`/`HV2`/`HV3` sequence strings are
deliberately **not consumed**: historical sequence generation is known to disagree
with the normalized variant lists for part of the cohort, so population statistics
trust the pipeline's declared intervals and normalized variant alleles instead.
Missing intervals are allowed but do not establish callable reference states.
Duplicate JSON keys, conflicting calls and unsupported schemas fail explicitly; they
are not silently repaired or aligned. Trace evidence and personal metadata are not
loaded into the analytical models.

Coordinates come solely from `get_settings().regions.REGIONS`, with all three HV
regions required. The configured `ref/rCRS.fasta` is used to validate representation,
not to call new variants. JSON and reference hashes identify the analyzed snapshot.

## Callability and exact profiles

- Coverage comes from the pipeline's declared intervals: a base is callable when the
  interval covers it and no unresolved variant call sits at that position. Positions
  outside every interval are not callable and are never treated as reference calls.
- Variant alleles carry the resolution evidence. Accepted uppercase IUPAC calls are
  exact upstream states: `Y` counts as `Y`, not separate C and T carriers, and matches
  only `Y`. `N` and lowercase alleles are unresolved: such a position (or insertion
  anchor) is not callable for that sample. Upstream must mask unresolved ambiguity.
  No heteroplasmy threshold or read-fraction averaging is introduced.
- A normalized entry that cannot be verified against the configured reference — a
  mis-referenced base, conflicting distinct calls at one position, a reference-equal
  entry or an unsupported allele representation — never counts as a carrier and
  leaves that position (or insertion anchor) non-callable for that sample, with an
  explicit `unverifiable_variant_representation` reason. It is not silently repaired
  or realigned. Structurally unlocatable entries (malformed or zero insertion
  positions) fail the run naming the sample; the export producer must fix those.
- A called deletion is a resolved state. An insertion is callable when its anchor and
  the next reference base are both callable, a single interval spans them, and the
  insertion's own allele is resolved. Insertions immediately after interval endpoints
  are therefore not callable. A represented insertion at the configured right boundary
  makes that region incomplete.
- Carrier counts include each callable normalized variant once per sample. The
  denominator is the number callable for that base/insertion anchor, not total N.
  Non-callable observations do not enter variant counts.
- Complete region profiles require the full configured span and every internal
  insertion anchor to be callable. Adjacent intervals that do not span their shared
  insertion boundary do not certify a complete profile. Combined
  profiles require all three complete regions. Incomplete profiles have null keys
  and explicit reasons; they never match complete reference-only profiles.
- Keys are sorted normalized `(position, reference, alternate)` tuples. Insertion
  indices sort numerically (`.9` before `.10`); duplicates count once. Space-separated
  display strings use the shared formatter (`249DEL`, `315.1C`). The existing `None`
  display sentinel means a **complete reference-only** profile, not missing data.
- Stable SHA-256 haplotype IDs bind canonical variants, comparison scope, reference
  hash and region-configuration hash. No poly-C or special-position exclusions,
  raw evidence reinterpretation, or approximate IUPAC matching is performed.
  The approved export must use one upstream normalized nomenclature across samples;
  sequence-equivalent but differently named indels are not renamed here.
- Identical plain profiles — the same sorted normalized variant keys plus the same
  per-region declared intervals — are computed once as a group and aggregated with
  the group's member count. Callable bases and insertion anchors derive from interval
  range arithmetic exactly equivalent to per-position membership, and validated
  export entries are model-checked in one bulk pass, so runtime scales with the
  number of distinct profiles instead of cohort size. No reported value changes.

## Outputs

| File | Contents |
|---|---|
| `mtDNA_population_summary.json` | Timestamp, eligibility, counts, variant summaries by region/all, complete combined and regional diversity, configuration/reference/input hashes, methods/options and checksums of the other six outputs |
| `mtDNA_variant_population_frequencies.tsv` | Variant/region/position/reference/alternate, carrier and callable counts, population frequency, rare flags, interval fields, optional per-database annotations |
| `mtDNA_haplotype_frequencies.tsv` | Complete combined scope, stable ID, canonical haplotype, count, frequency, singleton and interval fields |
| `mtDNA_region_haplotype_frequencies.tsv` | The same fields for each independently complete regional cohort |
| `mtDNA_sample_haplotypes.tsv` | Sample ID, per-region/combined complete profiles, combined ID, completeness map and reasons |
| `mtDNA_population_data_issues.tsv` | One row per sample, region and issue for fixable upstream representation problems: `missing_intervals` (no declared intervals), `duplicated_or_overlapping_intervals` (the offending pairs) and `unverifiable_variant_representation` (each entry with the reference base it contradicts) |
| `mtDNA_population_executive_summary.txt` | Human-readable current-cohort interpretation and identification limitations |

Tab-separated nested fields are deterministic JSON. Missing scalar values are blank
in TSV, null in JSON. Tables retain headers when empty. Only sample IDs, not names
or other personal metadata, are exported. Machine-readable values retain full
floating-point precision. Variants sort genomically, sample mappings by ID,
haplotypes by descending count then canonical genomic order.

Rare flags always include singleton, doubleton, carrier count ≤2 and strict
frequency `<0.001`, `<0.005`, `<0.01`. Additional explicit thresholds can be requested
through `AnalysisOptions.additional_frequency_thresholds`; mandatory flags do not
change. Summary keys use `variants.all`, `variants.HV1`, `variants.HV2`, `variants.HV3`.

## Statistics and interpretation

The combined complete cohort size is the haplotype denominator N. Given group counts n:

- Empirical random match probability: `sum(n*(n-1)) / (N*(N-1))`.
- Discrimination power: `1 - empirical_random_match_probability`.
- Frequency-squared sum (separate with-replacement statistic): `sum((n/N)**2)`.
- Sample-corrected haplotype diversity: `N/(N-1) * (1 - frequency_squared_sum)`.

Pairwise probability, DP and corrected diversity are null for N<2. For N=0,
frequency-squared sum and fractions without denominators are null; counts are zero.
For N=1, frequency-squared sum is 1. Group counts and frequencies are invariant-checked.
These statistics are computable directly from the exported haplotype counts; no
quadratic sample-to-sample comparison is performed.

Wilson binomial intervals are isolated in `statistics.frequency_interval`; no
pseudo-counts are added to reported frequencies. They are marginal sampling intervals,
not forensic match-probability confidence bounds. Cohort ascertainment, relatedness and
population sampling affect inference; these intervals do not adjust for them.
A singleton or an unobserved profile does not have a precisely known true population
frequency. mtDNA is maternally inherited: a complete profile match contributes evidence
but does not uniquely identify an individual or prove common specimen origin.

## Python API and optional reference annotation

```python
from src.modules.statistics.population import (
    AnalysisOptions, analyze_population, load_cohort, write_reports,
)

cohort = load_cohort("current-approved-export.json")
result = analyze_population(cohort, options=AnalysisOptions(confidence_level=0.95))
write_reports(result, "new-report-directory")
```

A typed `ReferenceProvider` may be passed as `reference_provider` to
`analyze_population`. Its `annotate` method receives a canonical `Variant`, reference
hash and actual region definitions, and returns a dictionary keyed by database name
of `ReferenceVariantAnnotation` objects. Each records `observed: bool | None`, optional
frequency/population, database version and provenance. Providers must validate reference,
nomenclature and scope compatibility themselves and perform a successful lookup before
returning `False`. Missing annotations and `None` mean unknown, never absent.

No live reference databases, network calls or external dependencies are required.
Only explicit `False` annotations support “not observed in selected reference database”.
Internal rarity is independent of reference absence; neither proves a globally novel mutation.
Provider exceptions propagate rather than being converted into false negative annotations.
