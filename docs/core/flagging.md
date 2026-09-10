# mtDNA Sample Flagging

This document explains how the mitochondrial DNA (mtDNA) sample flagging system
works. It is written for bioinformaticians and reviewers who need to understand
**what gets flagged, why, and what to do about it** — without reading the source
code.

Flagging is the quality-control layer that sits between variant calling and the
two-tool comparison. After two analysis methods (the pipeline and Sequencher)
each produce a variant list for a sample, the flagger inspects those variants and
marks anything that looks suspicious, unusual, or worth a human's attention. The
output is a small, structured set of *flag reasons* that the comparison report and
the downstream review workflow consume.

---

## At a glance

- **Where it lives:** `src/core/flagging.py`
- **Two kinds of output:** per-variant flags (`variant_flags`) and sample-level
  flags (`sample_flags`).
- **Severity:** every flag has a level from **0** (clean) to **4** (critical); a
  sample's level is the highest level among its flags.
- **Most flags are *presence* flags** — they fire when a particular variant is
  *present*. One flag, `No 315.1 variant`, is an *absence* flag — it fires when an
  expected variant is *missing*.
- **Exceptions are deliberate:** common, benign events (the expected `315.1C`
  insertion, the expected C-insertion at the noisy 309/573 polyC positions, the
  common 249 deletion, the paired 523-524 AC deletion) are **not** flagged, so the
  report focuses on genuinely unusual findings.

---

## Why we flag (the short version)

mtDNA has a few famously unreliable stretches — the polyC tracts around positions
309, 315, and 573, the AC repeat at 513-525, and the HV1 flank at 16180-16193.
These regions throw up artefacts and length-variation that the sequencers handle
inconsistently. Rather than treat every variant equally, the flagger:

1. **Suppresses** the known-benign events so they don't generate noise.
2. **Flags** the genuinely unusual events — a non-C insertion in a polyC tract, a
   deletion in a sensitive region, mixed bases (heteroplasmy), a run of consecutive
   indels that suggests an assembly error.
3. **Ranks** them by severity so a reviewer can triage.

The guiding principle throughout is: *only the expected, normal form is
suppressed; anything else is surfaced.* So `315.1C` (the rCRS reference insertion)
is suppressed, but `315.1T`, `315.2C`, and `315.3C` are all flagged. Likewise
`309.1C` and `573.1C` are suppressed, but `309.1T` and `573.1T` are flagged.

---

## How it works

The flagger has three parts, each with a clear job:

| Part | Job |
|---|---|
| **`VariantFlag`** | A tiny data container for one variant's flag results: its position, ref, alt, and the list of flag reasons that apply to it. |
| **`VariantAnalyzer`** | The "looking at individual variants" layer. It answers questions like "are there 5+ consecutive indels?", "is any base a low-confidence lowercase call?", "is the 459-deletion-plus-16192-issue complex pattern present?", and "does a variant exist at this position?" |
| **`SampleFlagger`** | The main entry point. You hand it one sample's variants and it returns three things: whether the sample should be flagged, the list of specific flag reasons, and a severity level. |

You usually don't call these directly. The two convenience entry points are:

- `analyze()` → `(should_flag, reasons, level)` for a whole sample.
- `flag_variants()` → a dict keyed by `"pos|ref|alt"` of per-variant flag reasons.

---

## Public API

### `VariantFlag`

```python
@dataclass
class VariantFlag:
    pos: int | str          # Variant position (normalized)
    ref: str                # Reference base
    seq: str                # Alternate sequence base
    flags: list[str]        # Flag reason strings for this variant

    def key() -> str         # Stable identity key: "pos|ref|alt"
```

### `VariantAnalyzer`

```python
VariantAnalyzer(variants: Variant | dict | list[dict[str, Any]])
```

Accepts a single `Variant`, a single dict, or a list of either, and normalizes
them internally. Its methods are the building blocks the sample flagger uses:

| Method | What it tells you |
|--------|-------------------|
| `has_consecutive_indels(count=5)` | Are there `count` or more *truly consecutive* same-type indels? (Consecutive deletion base positions, or consecutive insertion indices at one base.) |
| `has_lowercase_variants()` | Does any variant carry a lowercase base — a low-confidence call? |
| `check_combined_conditions()` | Is the 459-deletion + 16192-issue complex pattern present? |
| `has_variant_at(position, seq=None)` | Does a variant exist at `position` (optionally with a specific `seq`)? |
| `has_523_524_ac_deletion()` | Are both 523A→del and 524C→del present (the paired AC deletion)? |

### `SampleFlagger`

```python
SampleFlagger(variants: Variant | dict | list, intervals: dict | None = None)
```

The `intervals` argument is optional and only affects the `No 315.1 variant`
absence check (see [that section](#the-no-3151-variant-absence-flag)). When
omitted, full coverage is assumed.

| Method | Returns |
|--------|---------|
| `level_from_reasons(flag_reasons)` (classmethod) | The max severity across a list of flag reasons |
| `analyze()` | `(should_flag: bool, reasons: list[str], level: int)` |
| `flag_variants()` | `list[VariantFlag]` — one per input variant |
| `classify(fis_flag, fis_flag_info, flag_reasons, genotypes_str)` | A severity level 0-4, used by the FIS comparison path |

### `flag_variants()` — the unified entry point

```python
flag_variants(variants: Variant | dict | list) -> dict[str, list[str]]
```

Returns a dict keyed by `"pos|ref|alt"` mapping each flagged variant to its list of
flag reasons. This is the canonical way to get per-variant flags and is used across
all tools.

---

## Severity levels (0-4)

A sample's severity is the **maximum** level across all its flag reasons. The
levels are defined in `SampleFlagger.FLAG_REASON_MAP` (pattern → level); the
`classify()` method additionally consults the FIS QC info and genotype string. The FIS tool supplies nomenclature-corrected, reference-backed variants and the same selected profile string; R2 QC remains an independent quality input.

| Level | Name | What it means | Example flags |
|---|---|---|---|
| **0** | None | Nothing flagged | A clean sample |
| **1** | Low | Benign or common; unlikely to affect interpretation | A variant in the 16180-16193 polyC flank; a lone `Deletion at 523/524`; `Has 309T variant` |
| **2** | Medium | Technical or nomenclature worth a look, not necessarily a problem | An AC-repeat SNP (`Has 519A-G variant`); an insertion at a named polyC-flank position (290.1, 16182.1, 16193.1, 16261.1, 16296.1); `Has 309Y variant`; `Has 460 variant` |
| **3** | High | Complex or difficult; may affect interpretation | Any general `Insertion X at Y` or `Deletion at X`; `Heteroplasmy at X`; `No 315.1 variant`; `Has 310 variant`; the transversion flags (`16293A-C`, `16258A-C`, `16263T-C`); `Deletion at 315/459/513/514/16182/16193/309`; `Consecutive indels (5+)`; `Lowercase variant detected`; `Complex: 459 deletion with 16192 issue` |
| **4** | Critical | Quality issue that likely compromises the result | The QC keywords `Low depth`, `Failed position`, `Possible wrong consensus` (these come from the FIS flag-info string, not from the variants themselves) |

> **Why the asymmetry between 310 and 460?** Both are position-presence flags, but
> 310 sits in a region historically associated with review-worthy complexity, so it
> is Level 3, while 460 is Level 2. The levels are not a judgement of "how bad is
> this base" but of "how much does this deserve a reviewer's attention given where it
> is."

---

## The complete flag catalog

Every flag string the flagger can produce is listed below, grouped by kind. The
**Output** column says where the flag lives: `V` = per-variant (`variant_flags`),
`S` = sample-level (`sample_flags`). The **Level** is the severity from
`FLAG_REASON_MAP`.

### Insertions and deletions

These are the workhorse flags — one per indel that isn't suppressed by an exception.

| Flag string | Level | Output | When it fires |
|---|---|---|---|
| `Insertion {seq} at {pos}` | 3 — High | V + S | Any insertion that is *not* suppressed and *not* at one of the named Level-2 positions below |
| `Insertion {seq} at 290.1` / `290.2` | 2 — Medium | V + S | Insertion at these named polyC-flank positions |
| `Insertion {seq} at 16182.1` / `16192.1` / `16193.1` / `16261.1` / `16296.1` | 2 — Medium | V + S | Insertion at these named HV-region positions |
| `Insertion {seq} at 16188.1` | 3 — High | V + S | Insertion at 16188.1 (named High) |
| `Deletion at {pos}` | 3 — High | V + S | Any deletion that isn't suppressed and isn't a Level-1 named position |
| `Deletion at 523` / `Deletion at 524` | 1 — Low | V + S | A *lone* AC-base deletion (only one of the 523/524 pair is present) |
| `Deletion at 315` / `459` / `513` / `514` / `16182` / `16193` / `309` | 3 — High | V + S | Deletion at these named positions |

**Example.** A single `567.1C` insertion produces:
`['Insertion C at 567.1', 'No 315.1 variant']`. Five consecutive deletions at
100-104 produce five `Deletion at …` flags *plus* `Consecutive indels (5+)`.

### Heteroplasmy and low-confidence calls

| Flag string | Level | Output | When it fires |
|---|---|---|---|
| `Heteroplasmy at {pos}` | 3 — High | V + S | `seq` contains an IUPAC ambiguity code (`R Y S W K M B D H V N`) — a mixed-base call |
| `Lowercase variant detected` | 3 — High | V + S | `seq` contains a lowercase letter — a low-confidence base call |

> Heteroplasmy flags fire wherever the ambiguity code appears; they are not gated
> to any particular region. The Mutation Surveyor pipeline normalizes its `_het`
> tokens to IUPAC codes (e.g. `309C_het` → `seq="Y"`) before flagging, so core
> detects `Heteroplasmy at 309` automatically.

### Region and absence flags

| Flag string | Level | Output | When it fires |
|---|---|---|---|
| `16180-16193 region` | 1 — Low | V | Per-variant: the base position falls in 16180-16193 (the HV1 polyC flank) |
| `16180-16193 region (p1, p2, …)` | 1 — Low | S | Aggregated summary of all the region hits, sorted and deduplicated |

The `No 315.1 variant` flag is special enough to get its own section — see
[below](#the-no-3151-variant-absence-flag).

### Position-specific variant flags

These fire when a variant lands on a position of particular interest.

| Flag string | Level | Output | When it fires |
|---|---|---|---|
| `Has 309T variant` | 1 — Low | V + S | A SNP at 309 with `seq=T` |
| `Has 309Y variant` | 2 — Medium | V + S | A SNP at 309 with `seq=Y` (heteroplasmy at 309) |
| `Has 310 variant` | 3 — High | V + S | Any variant at position 310 |
| `Has 460 variant` | 2 — Medium | V + S | Any variant at position 460 |
| `Has 16293A-C variant` | 3 — High | V + S | Transversion at 16293: ref=A, seq=C |
| `Has 16258A-C variant` | 3 — High | V + S | Transversion at 16258: ref=A, seq=C |
| `Has 16263T-C variant` | 3 — High | V + S | Transition at 16263: ref=T, seq=C |
| `Has {pos}{ref}-{seq} variant` | 2 — Medium | V + S | A SNP in the AC repeat region 513-525 (e.g. `Has 519A-G variant`) |

### Pattern and complex flags

These look across multiple variants rather than at a single one.

| Flag string | Level | Output | When it fires |
|---|---|---|---|
| `Consecutive indels (5+)` | 3 — High | S | 5+ *truly consecutive* same-type indels — consecutive deletion base positions (e.g. 100,101,102,103,104) or consecutive insertion indices at one base (e.g. 100.1,100.2,100.3,100.4,100.5). Indels that are merely adjacent in the sorted list but don't meet these criteria do **not** trigger it. |
| `Complex: 459 deletion with 16192 issue` | 3 — High | S | A 459 deletion paired with a 16192 issue — either heteroplasmy at 16192, or `seq=T`/an ambiguity code there. This combination is a known complex pattern that needs expert review. |

### Quality flags (from `classify()`, not `analyze()`)

These come from the FIS comparison path, read out of the `fis_flag_info` string —
they are *not* produced by inspecting variants.

| Flag string | Level | When it fires |
|---|---|---|
| `Low depth` | 4 — Critical | The keyword appears in `fis_flag_info` |
| `Failed position` | 4 — Critical | The keyword appears in `fis_flag_info` |
| `Possible wrong consensus` | 4 — Critical | The keyword appears in `fis_flag_info` |

---

## The `No 315.1 variant` absence flag

This is the only *absence* flag — it reports a variant that **should** be there but
isn't. The rCRS reference carries a `315.1C` insertion in the HV2 polyC stretch
(303-315); almost every real sample has it. When it's missing, something is off and
the flag fires.

It is **region-aware**, which is the subtle part:

- It fires only when `315.1` is absent **and** the HV2 region (the region containing
  position 315) was actually analyzed.
- If a sample only sequenced HV1 (HV2 was never run), `315.1` can't be expected, so
  the flag does **not** fire. A sample that sequenced HV2 but left a gap over 303-315
  (a data problem) *still* gets flagged — the gap doesn't excuse the missing call.
- `intervals=None` means "assume full coverage" (the legacy/FIS path). An empty
  intervals dict `{}` (no region analyzed) does **not** flag — it is not the same
  as full coverage.

Because the flag depends on the actual variant set and coverage, it is recomputed
at JSON-write time by `recompute_no_315_1_flag()` rather than trusted from stale
per-file or per-region output. This keeps each output JSON (the combined 3-region
file and each per-region file) consistent with its own contents.

---

## Sample-level vs per-variant flags

Flag reasons are split across two output fields so the comparison report can show
both a position-level view and a sample-level view without duplicating:

- **`variant_flags`** — per-variant, keyed by `"pos|ref|alt"`. Holds
  `Insertion {seq} at {pos}`, `Deletion at {pos}`, `Heteroplasmy at {pos}`, the
  AC-repeat `Has {pos}{ref}-{seq} variant`, `Has 309T/Y variant`, `Has 310/460
  variant`, the transversion flags, the per-variant `16180-16193 region`, and
  `Lowercase variant detected`.
- **`sample_flags`** — sample-level only. Holds the absence check
  `No 315.1 variant` (region-aware), the pattern flags
  (`Consecutive indels (5+)`, `Complex: 459 deletion …`), and tool-level QC flags.
  Per-variant flag strings and the aggregated `16180-16193 region (…)` summary are
  **not** kept here.

To keep the two from overlapping, `deduplicate_sample_flags()` drops any
sample-level flag whose exact string already appears as a per-variant flag, and
also drops the aggregated `16180-16193 region (…)` form whenever the per-variant
`16180-16193 region` flag is present (the comparison report consolidates the
per-variant entries into the same aggregated form for its Variant Flags column, so
the summary isn't duplicated in `sample_flags`).

---

## Exception rules — what is *not* flagged, and why

Most indels are flagged. The exceptions below are the **only** cases where an indel
is suppressed, and each has a biological reason.

| Position | Type | Rule | Why |
|----------|------|------|-----|
| 315.1 | **C**-insertion only | **Not flagged** | `315.1C` is the expected rCRS reference insertion; nearly every sample has it. |
| 315.1 (non-C) | Insertion | **Flagged** | Only the C form is normal; another base at 315.1 is unusual. |
| 315.2 / 315.3 / … | Insertion (any base) | **Flagged** | Extra insertions beyond the expected 315.1C are real length variants. |
| 309, 573 | **C**-insertion only | **Not flagged** | The expected C-insertion at these noisy polyC positions is a known artefact. |
| 309.1 / 573.1 (non-C) | Insertion | **Flagged** | Only the C form is the expected artefact; another base warrants review. |
| 455.1 / 463.1 | Insertion (any base) | **Flagged** | Special positions, but not C-exception sites, so insertions are flagged. |
| 16193.1 | Insertion (any base, incl. C) | **Flagged** | 16193 is deliberately *not* in the C-exception set; `16193.1C` is flagged. |
| 249 | Deletion | **Not flagged** | A common, well-known deletion. |
| 523-524 (both present) | AC deletion pair | **Not flagged** | The common AC-repeat deletion is only benign as a *pair*. |
| 523 or 524 alone | Single AC-base deletion | **Flagged** (Level 1) | Only the paired deletion is benign; a lone one is flagged. |
| 513A + 523-524 deletion | Combo | **Not flagged** (the 513A SNP is suppressed) | Linked pattern; the 513A→X SNP is suppressed when the AC deletion is present. |
| 309 | Deletion (base position, `309DEL`) | **Flagged** | A real variant event (e.g. the Tracy `309C>T` + `310T>C` conversion). |
| 455 / 463 / 573 | Deletion (base position) | **Not flagged** | Special-position deletions are skipped; only `309DEL` is kept. |

The recurring theme: **only the expected normal form is suppressed.** `315.1C`,
`309.1C`, `573.1C`, the 249 deletion, and the paired 523-524 AC deletion are the
"expected" events; everything else is surfaced.

---

## Special positions — flagging vs concordance

Positions 309, 455, 463, and 573 (and their sub-positions 309.1, 309.2, 573.1, …)
are "special positions" **for concordance**: variants there are ignored when
deciding whether the two tools agree, because the polyC stretches are too
unreliable for a fair comparison.

**Flagging is different from concordance.** The special-position guard only skips
special-position *deletions*, not insertions:

```python
if is_special_position(pos, ref, seq) and ref != "-" and not self._is_309_deletion(pos, seq):
    return ...   # skip: special-position deletion (and non-indel SNP)
```

What this means in practice:

- **Insertions at special positions are checked.** The guard's `ref != "-"` lets
  them through to the insertion check, where the C-exception suppresses only the
  expected C form:
  - `309.1C` / `573.1C` → **not flagged** (expected polyC artefact)
  - `309.1T` / `573.1T` / `573.1A` … → **flagged** (`Insertion {seq} at 309.1/573.1`)
  - `455.1C` / `463.1C` → **flagged** (special positions, but not C-exception sites)
- **Deletions at special positions are skipped**, with one exception:
  - `309DEL` (a base-position 309 deletion) → **flagged** as `Deletion at 309` —
    a real variant event (e.g. the Tracy `309C>T` + `310T>C` conversion).
  - Deletions at 455 / 463 / 573 (base) and at sub-positions like 309.1 → **skipped**.
- **Non-indel SNPs at special positions** are not "special" for flagging at all
  (`is_special_position` returns `False` for non-indels), so `309T` / `309Y` are
  still flagged as `Has 309T variant` / `Has 309Y variant`.

> This is intentionally narrower than the legacy comparator, which blanket-skipped
> *all* special-position indels — including non-C insertions, which masked genuine
> events. Suppressing only the expected C form surfaces the non-C insertions and
> the `309DEL` event. (Per TASKS.md #2.)

---

## Worked examples

A few real inputs and the flags they produce (the `No 315.1 variant` absence flag
appears in most because these toy inputs don't include `315.1C` and assume full
coverage):

**A single HV1-region SNP** (`16185 A>G`, `16189 T>C`):
```
['16180-16193 region (16185, 16189)', 'No 315.1 variant']
```
Two variants in the 16180-16193 flank are consolidated into one aggregated region
flag.

**A 459 deletion + heteroplasmy at 16192** — the complex pattern:
```
['Deletion at 459',
 'Heteroplasmy at 16192', '16180-16193 region (16192)',
 'No 315.1 variant', 'Complex: 459 deletion with 16192 issue']
```
The 16192 variant is in the flank *and* heteroplasmic, and together with the 459
deletion it triggers the complex-pattern flag.

**Five consecutive deletions** (100-104):
```
['Deletion at 100', 'Deletion at 101', 'Deletion at 102',
 'Deletion at 103', 'Deletion at 104', 'No 315.1 variant',
 'Consecutive indels (5+)']
```
Each deletion gets its own flag *and* the run triggers the pattern flag.

**A real sample — `LN_26_AA0050`** (Sequencher side):
```
variants: 73G 263G 315.1C 523DEL 524DEL 573.1C 16181C 16182C 16183C 16189C 16213A 16217C 16261T 16292T
flags:    ['16180-16193 region (16181, 16182, 16183, 16189)']
```
Note what is **not** here: `573.1C` produces **no** insertion flag (the expected C
insertion is suppressed), `315.1C` produces nothing (expected, and present so no
absence flag), and the paired `523DEL 524DEL` produces no deletion flag. Only the
genuinely review-worthy 16180-16193 flank variants surface.

---

## Configuration constants

| Constant | Value | Purpose |
|----------|-------|---------|
| `IUPAC_AMBIGUITY_CODES` | `{R, Y, S, W, K, M, B, D, H, V, N}` | IUPAC codes that indicate heteroplasmy |
| `FLAGGING_REGION_START` | `16180` | Start of the HV1 flagging region |
| `FLAGGING_REGION_END` | `16193` | End of the HV1 flagging region |
| `AC_REPEAT_REGION_START` | `513` | Start of the AC repeat region |
| `AC_REPEAT_REGION_END` | `525` | End of the AC repeat region |
| `INSERTION_C_EXCEPTION_POSITIONS` | `{309, 573}` | Positions where the expected **C**-insertion is suppressed. 16193 is deliberately *excluded* so `16193.1C` is flagged. `315.1C` is handled separately by the `normalize_position(pos) == "315.1" and seq == "C"` check, not this set. |
| `_POSITION_VARIANT_FLAGS` | `[("309","T","Has 309T variant"), ("309","Y","Has 309Y variant")]` | Position+seq-specific flags |
| `_POSITION_FLAGS` | `[("310","Has 310 variant"), ("460","Has 460 variant")]` | Position-presence flags |
| `_TRANSVERSION_CHECKS` | `[(16293,A,C), (16258,A,C), (16263,T,C)]` | Transversion/transition flags |
| `POS_*` | `249, 309, 459, 513, 523, 524, 16192, 16258, 16263, 16293` | Well-known mtDNA positions used in the logic |

> `POS_315` was removed when the 315 suppression was narrowed from "all base-315
> insertions" to "only `315.1C`".

---

## Current system vs legacy

The current flagger catches every variant condition the original inline logic did,
and adds several flags beyond it (`16258A-C`, `16263T-C`, AC-repeat SNPs,
consecutive-indel detection, lowercase detection, severity levels, and per-variant
reasons). No regressions.

Key improvements over the legacy inline logic:

- Het deletions (`seq="DEL"`) and het insertions (lowercase `seq`) now produce
  `Heteroplasmy at {pos}` flags.
- `16258A-C` and `16263T-C` variants are flagged.
- SNPs in the AC repeat region (513-525) are flagged.
- Consecutive indels (5+) are flagged.
- Lowercase (low-confidence) calls are flagged.
- Severity levels (0-4) let reviewers triage.
- Per-variant flag reasons enable position-level comparison reports.

The one deliberate divergence from the legacy *comparator*: the comparator
blanket-skipped all special-position indels (including non-C insertions). The
current flagger suppresses only the expected C form, so non-C insertions at
309/573 and the `309DEL` event are surfaced. See
[Special positions](#special-positions--flagging-vs-concordance).

---

## Cross-tool variant flagging

`flag_variants()` is the unified entry point used across all tools:

- **Mutation Surveyor** normalizes `_het` tokens to IUPAC codes in `Variant.seq`
  (e.g. `309C_het` → `seq="Y"`), so core detects `Heteroplasmy at 309` automatically.
- Het deletions (`seq="DEL"`) and het insertions (lowercase `seq`) also produce
  `Heteroplasmy at {pos}`.
- An explicit `"het"` flag (review-triggering, sample-level) is added by
  `build_hv1_flags()` / `build_hv23_flags()` when het variants are in the merged
  profile; source het flags (`HV1 source het`, `HV1 source polyC het`) come from
  `source_het_flags()` when source traces contain het variants.
- Insertions and deletions are flagged by core with the exception rules above
  (only `315.1C`, the expected C-insertion at 309/573, the 249 deletion, and the
  paired AC 523-524 deletion are suppressed; everything else is flagged).
- Tool-specific sample-level flags (coverage, pair status, region status) stay in
  each module's flagging code.

Tool-specific flags that remain in the modules:

- **PHP**, **PolyC informational**, **HV1/HV2-3 region flags** →
  `src/tools/mutation_surveyor/flagging.py`
- Sample-level routing (`sample_flags` vs `variant_flags`) →
  `src/tools/mutation_surveyor/etl.py`

> **Deferred:** PHP flag routing and PolyC informational flags may move to core in a
> future update, pending lab/clinical input.

---

## Future extensions

The flagger's per-variant and per-sample signals feed the downstream routing
workflow (see [`docs/tools/mutation_surveyor.md`](../tools/mutation_surveyor.md)):

- **Region-based auto-pass:** `AUTOPASS_HV1` and `AUTOPASS_HV2_AND_HV3` flags (in
  `src/tools/mutation_surveyor/flagging.py`) combined with coverage logic determine
  per-region pass/review routing.
- **Het filtering in final output:** het variants identified by
  `Heteroplasmy at {pos}` are kept in comparison/audit outputs but excluded from the
  final Excel/JSON/FASTA via a filter at export time.
- **Tool priority for multi-source merge:** when several tools cover the same
  region, configurable priority rules owned by `src/core/mtdna_merger.py` pick the
  best source.

See [`docs/tools/mutation_surveyor.md`](../tools/mutation_surveyor.md) for the full
MS AutoRun + MS1 target workflow.
