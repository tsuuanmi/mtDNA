# Tracy Modular QC and Trusted-Coverage Architecture

> **Status**: Design proposal. This document defines the intended evolution of
> Tracy quality control; it does not describe behavior that has already been
> fully implemented.
>
> **Current implementation**: See [tracy-noise-qc.md](tracy-noise-qc.md).
>
> **Scope**: Tracy/Sanger trace quality assessment, per-trace coverage masking,
> forward/reverse evidence merging, QC findings, and sample-level QC flags.

## 1. Goal

Evolve the existing Tracy noise QC into a modular, optional trace-evidence QC
layer that can distinguish different Sanger failure modes while preserving the
current strengths of per-trace masking.

The central rule is:

```text
aligned coverage != trusted coverage
```

A trace may align through a position without providing evidence that is safe to
interpret. Conversely, one unreliable trace must not invalidate a position when
another independent trace provides trusted evidence.

The module must therefore keep these concepts separate:

- **Raw coverage**: positions reached by the Tracy-derived alignment.
- **Trusted coverage**: positions whose evidence is accepted for final
  interpretation after enabled QC policies are applied.
- **Observed candidates**: variants or other events detected from the trace,
  including candidates that are later excluded or sent for review.
- **Accepted variants**: candidates supported by trusted evidence and allowed
  into normal final output.
- **QC findings**: structured explanations of trace-local uncertainty.
- **Sample flags**: final unresolved issues after trace evidence has been merged.

## 2. Why this is needed

### 2.1 The existing noise mask solves an important coverage problem

The current Tracy noise QC already does more than classify a trace as noisy.
For a `likely_noisy` range it:

1. removes variants overlapping that range from the per-trace `Sample`;
2. subtracts the same range from that trace's coverage intervals;
3. allows another clean trace to restore the removed coverage during the
   same-sample merge.

This is the correct semantic pattern for untrusted evidence.

For example:

```text
Raw Tracy coverage:
73 -------------------------------------------------- 340

Problematic internal range:
                    180 -------- 190

Trusted coverage after masking:
73 --------------- 179     191 ---------------------- 340
```

Tracy-derived coverage is fundamentally bounded by a start and end, while the
post-ETL mask can split one aligned span into multiple trusted spans. The mask is
therefore best understood as an **advanced coverage trim** rather than only a
"noise filter".

### 2.2 Variant removal without coverage removal is unsafe

Some existing Tracy transformations can reject a candidate variant while the
trace interval still covers that position.

That creates an ambiguity:

```text
covered position + no retained variant
```

can later be interpreted as an implicit reference call even when the real
meaning was:

```text
trace evidence at this position is not trustworthy
```

The QC architecture must ensure that an exclusion caused by unreliable evidence
can also remove that evidence from the originating trace's trusted coverage.

### 2.3 Different artifacts need different evidence

Generic signal noise, read-edge uncertainty, homopolymer/polyC phase problems,
and future artifacts such as dye blobs do not have identical signatures.

A single global classifier should not be forced to explain all of them.

The system should instead support focused assessors that share common trace
metrics but apply artifact-specific reasoning.

## 3. Design principles

### 3.1 Per-trace first, merge second

QC is evaluated independently for each Tracy trace before forward/reverse traces
are combined.

```text
HV2F -> QC -> trusted HV2F evidence --\
                                      +-> merge -> final Sample
HV3R -> QC -> trusted HV3R evidence --/
```

This preserves provenance and allows one trusted mate trace to rescue evidence
that was rejected in another trace.

### 3.2 Detection is separate from action

An assessor detects and describes a problem. A policy decides what to do with
that finding.

```text
trace evidence
    -> assessor
    -> QCFinding
    -> policy
       -> ACCEPT
       -> REVIEW
       -> EXCLUDE
```

Assessors must not directly mutate `Sample` coverage or delete variants.

This separation allows a rule to begin as `REVIEW`, be validated independently,
and later become `EXCLUDE` without rewriting the detector.

### 3.3 Masking is a generic operation

Coverage subtraction must not belong specifically to the signal-noise
classifier.

Any finding that resolves to `EXCLUDE` may contribute an exclusion range:

```text
SIGNAL_NOISE                    -> trim
READ_EDGE_UNCERTAINTY           -> trim
POLYC_DIRECTIONAL_UNCERTAINTY   -> trim
DYE_BLOB                        -> trim   # future
```

The generic operation is:

```text
trusted_coverage = raw_coverage - exclusion_ranges
```

It must support ranges at the beginning, end, or middle of a trace.

### 3.4 Raw observations are preserved

QC changes interpretation, not history.

A candidate that is excluded from normal output should remain available in QC
metadata with its reason and source trace so that it can be audited or
re-evaluated later.

For example:

```text
310 C>T
status: EXCLUDED
reason: POLYC_DIRECTIONAL_UNCERTAINTY
trace: HV2F
```

This is intentionally different from silently deleting the observation.

### 3.5 The QC layer remains optional

The QC layer is an add-on to Tracy processing, not a dependency of the core
canonical model or downstream consumers.

The dependency direction should remain:

```text
Tracy decompose / ETL
        -> trace evidence
        -> optional Tracy QC
        -> per-trace Sample
        -> same-LID merge
        -> canonical Sample
```

Tracy ETL must not import individual QC assessors.

## 4. Proposed QC architecture

```text
AB1
 |
 v
Tracy decompose
 |
 v
Tracy ETL / raw per-trace Sample
 |
 +----------------------------------------------------+
 |                                                    |
 | QC disabled                                        | QC enabled
 |                                                    v
 |                                                QC Engine
 |                                                    |
 |                        +---------------------------+-------------------+
 |                        |                           |                   |
 |                        v                           v                   v
 |                 SignalNoiseAssessor         PolyCAssessor      ReadEdgeAssessor
 |                                                                      ...
 |                                                    |
 |                                                    v
 |                                               QC Findings
 |                                                    |
 |                                                    v
 |                                                QC Policy
 |                                  +-----------------+-----------------+
 |                                  |                 |                 |
 |                                  v                 v                 v
 |                               ACCEPT             REVIEW            EXCLUDE
 |                                                                        |
 |                                                                        v
 |                                                        generic coverage masking
 |                                                                        |
 +------------------------------------------------------------------------+
                                  |
                                  v
                         interpreted per-trace Sample
                                  |
                                  v
                      merge trusted trace evidence
                                  |
                                  v
                              final Sample
                                  |
                                  +-- trusted variants
                                  +-- trusted coverage
                                  +-- existing sample_flags
```

The exact internal class names may change during implementation. The boundaries
and responsibilities above are the contract.

## 5. QC finding model

A trace-local QC finding should be structured rather than encoded only as a
free-form string.

Conceptually:

```text
QCFinding
{
    code: READ_EDGE_UNCERTAINTY,
    region: HV1,
    interval: [16024, 16031],
    traces: ["HV1F"],
    decision: REVIEW | EXCLUDE,
    metrics: {...}
}
```

Required semantics:

- `code`: stable machine-readable reason.
- `region`: canonical mtDNA region when known.
- `interval`: inclusive canonical reference interval affected by the finding.
- `traces`: one or more source trace identifiers or primer names.
- `decision`: policy result for this finding.
- `metrics`: optional diagnostic evidence used to create the finding.

The initial implementation does not need to change the canonical `Sample` model
to store this structured object. Trace-local findings may remain in the Tracy QC
report while unresolved final conditions are rendered into the existing
`Sample.sample_flags` system.

## 6. Assessor responsibilities

### 6.1 Signal-noise assessor

This assessor should preserve the behavior that is already working in
`quality_control.py`:

- signal;
- peak-separation SNR;
- purity;
- background ratio;
- secondary-peak ratio;
- Tracy base quality;
- weak-signal fraction;
- overlapping window analysis;
- sustained support before producing a maskable range.

The first implementation should be behavior-preserving: refactoring the
classifier into a modular assessor must not silently change existing masking
thresholds or final output.

### 6.2 Read-edge assessor

Read-edge uncertainty is not the same as generic sustained noise.

A read-edge assessor may use:

- distance from the usable Tracy alignment edge;
- local per-base metrics;
- shorter edge-specific windows;
- primer/direction context;
- validation-derived boundaries.

Example finding:

```text
code: READ_EDGE_UNCERTAINTY
region: HV1
interval: [16024, 16031]
traces: ["HV1F"]
```

The assessor must be able to describe only the affected edge span rather than
masking a larger generic 15-base window when the evidence does not support it.

### 6.3 PolyC/homopolymer assessor

Homopolymer uncertainty must be treated separately from weak-signal noise.
Strong peaks or high base quality do not by themselves prove that homopolymer
length or downstream sequence phase is correct.

The polyC assessor may use:

- known polyC context;
- primer direction;
- local purity changes;
- secondary-peak behavior;
- downstream phase deterioration;
- forward/reverse consistency;
- validated position-specific policies.

Existing directional polyC exclusions are candidates for migration into this
model because they describe evidence that is intentionally not trusted from one
trace direction.

A directional exclusion should affect the originating trace's trusted coverage,
not automatically the whole sample.

### 6.4 Future artifact assessors

The architecture must allow additional assessors without changing Tracy ETL,
merge logic, or the canonical downstream API.

Examples include:

```text
DYE_BLOB
SIGNAL_SATURATION
BASELINE_DISTURBANCE
```

Each new assessor should only need to:

1. inspect trace evidence;
2. emit structured findings;
3. register its policy/configuration;
4. reuse the common masking and reporting machinery.

## 7. Generic coverage trimming

The existing ability to subtract internal noise ranges should become a shared QC
primitive.

Input:

```text
raw intervals + exclusion ranges
```

Output:

```text
trusted intervals
```

Examples:

### Edge trim

```text
raw:       16024 ------------------------------ 16365
exclude:   16024 -- 16031
trusted:           16032 ---------------------- 16365
```

### Internal trim

```text
raw:       73 ------------------------------------ 340
exclude:                    180 -- 190
trusted:   73 ------------ 179     191 ----------- 340
```

### Multiple exclusions

```text
raw:       73 ------------------------------------ 340
exclude:       90-95       180-190       300-305
trusted:   73-89  96-179    191-299       306-340
```

This operation should remain independent of why each range was excluded.

## 8. Review versus exclusion

A finding does not automatically remove coverage.

### REVIEW

`REVIEW` means the system found evidence worth surfacing but the current
validation does not justify automatic exclusion.

Default behavior:

- retain the existing trusted coverage;
- retain normal output unless another rule excludes it;
- record the finding in Tracy QC metadata;
- promote a final unresolved review condition to a sample flag when appropriate.

### EXCLUDE

`EXCLUDE` means the affected trace evidence must not contribute to final
interpretation.

Behavior:

- remove overlapping accepted variants from that trace's normal output;
- retain excluded observations in QC metadata;
- subtract the affected interval from that trace's trusted coverage;
- recompute affected trace/sample flags before merge;
- allow another trusted trace to rescue the position.

This distinction lets exploratory polyC/read-edge rules start in review-only
mode and move to automatic exclusion only after independent validation.

## 9. Merge semantics

Merge must operate on trusted per-trace evidence.

### 9.1 Trusted plus untrusted

```text
HV2F @310: untrusted G
HV3R @310: trusted C

final: C
```

This is not a conflict. Only the trusted observation contributes to final
interpretation.

### 9.2 Both untrusted

```text
HV2F @310: untrusted
HV3R @310: untrusted

final: NO_CALL at 310
```

The final sample must not infer the reference allele merely because both
untrusted candidate calls were removed.

### 9.3 Both trusted and concordant

```text
HV2F @263: trusted G
HV3R @263: trusted G

final: G
```

No conflict flag is required.

### 9.4 Both trusted and discordant

```text
HV2F @263: trusted G
HV3R @263: trusted A

final: unresolved trusted-evidence disagreement
```

Do not create a parallel conflict subsystem. Reuse the existing sample flagging
mechanism and add an appropriate sample-level flag for the disagreement.

The final flag should contain enough context in its rendered reason or linked QC
metadata to identify the affected region/interval/traces.

Conceptually the source finding is:

```text
code: TRACE_CONFLICT
region: HV2
interval: [263, 263]
traces: ["HV2F", "HV3R"]
```

The canonical `Sample.sample_flags` representation remains the existing flag
system unless a separate schema change is explicitly approved later.

## 10. Trace findings versus sample flags

Not every trace warning should become a final sample warning.

Trace-level findings describe what happened to an individual source trace.
Sample-level flags should describe unresolved conditions that still matter after
all trusted trace evidence has been merged.

Example:

```text
HV2F @310: EXCLUDE because of polyC uncertainty
HV3R @310: trusted and clean
```

The HV2F finding remains in the QC report for audit, but the final sample does
not need a conflict/no-coverage flag at 310 because HV3R rescued it.

By contrast:

```text
HV2F @310: EXCLUDE
HV3R @310: EXCLUDE
```

leaves a real final coverage problem and should produce the appropriate
sample-level quality flag.

Likewise:

```text
HV2F @263: trusted G
HV3R @263: trusted A
```

leaves a final disagreement and should produce an existing-system sample flag.

## 11. FULL REGION semantics

`FULL REGION` should mean that the merged trusted evidence covers the complete
configured region, not merely that Tracy aligned through it.

Example:

```text
HV2F trusted: 73 -------- 303
HV3R trusted:       300 ---------------- 340

merged trusted coverage: 73 ------------ 340
=> FULL REGION
```

But:

```text
HV2F trusted: 73 -------- 303
HV3R trusted:                  316 ------ 340

304-315 has no trusted evidence
=> NOT FULL REGION
```

This prevents excluded/untrusted positions from being silently represented as
trusted reference sequence.

## 12. Optional configuration

The QC layer must be switchable as a module and, eventually, by assessor.

Conceptual configuration:

```yaml
qc:
  enabled: true

  signal_noise:
    enabled: true

  polyc:
    enabled: true

  read_edge:
    enabled: true

  dye_blob:
    enabled: false
```

The exact configuration format should follow existing project settings when
implemented.

Required behavior:

- `qc.enabled = false`: Tracy continues through its normal non-QC path.
- disabling one assessor does not disable the others;
- disabling an assessor does not destroy or mutate raw trace evidence;
- reporting may still expose available raw diagnostics when masking is disabled,
  following the existing noise-QC pattern;
- the optional layer must not create a reverse dependency from Tracy ETL into
  assessor implementations.

## 13. Proposed module boundaries

A future refactor may use a layout similar to:

```text
src/tools/tracy/qc/
+-- models.py
+-- engine.py
+-- policy.py
+-- coverage.py
+-- signal_noise.py
+-- polyc.py
+-- read_edge.py
+-- reporting.py
```

Responsibilities:

### `models.py`

Shared Tracy-QC data contracts such as base/window evidence, `QCFinding`, and
policy decisions.

### `engine.py`

Runs the enabled assessors and collects findings. It contains orchestration, not
artifact-specific biological rules.

### `policy.py`

Maps validated findings to `ACCEPT`, `REVIEW`, or `EXCLUDE`.

### `coverage.py`

Owns the generic inclusive-range subtraction used to derive trusted intervals
from raw intervals. It must not know whether the source finding was noise,
polyC, a read edge, or another artifact.

### `signal_noise.py`

Owns the existing generic signal-noise metrics and window classifier after a
behavior-preserving extraction from `quality_control.py`.

### `polyc.py`

Owns polyC/homopolymer-specific detection and directional uncertainty rules.

### `read_edge.py`

Owns edge-specific evaluation and edge range findings.

### `reporting.py`

Serializes trace findings, exclusions, metrics, and provenance into the Tracy QC
report without adding dependencies to canonical downstream models.

This layout is illustrative. Implementation should prefer the smallest set of
files that creates clear boundaries and should not introduce empty abstractions
before they are needed.

## 14. Relationship to the current noise QC

The current implementation should be treated as the first working instance of
the broader model:

```text
current NoiseRange
    -> current masking
    -> split per-trace intervals
```

becomes conceptually:

```text
SIGNAL_NOISE finding
    -> EXCLUDE policy
    -> generic QC coverage mask
```

The initial refactor should preserve existing noise behavior and tests. The goal
is to generalize ownership of the masking operation, not to retune the current
noise classifier at the same time.

Current behavior that must remain valid:

- masking happens before same-LID trace merge;
- only the originating trace loses the excluded range;
- a trusted mate trace can restore final coverage;
- variants overlapping excluded ranges do not enter normal final output;
- excluded observations remain visible in QC metadata;
- disabling masking restores unmasked output while QC reporting remains
  available according to current configuration behavior.

## 15. Candidate migration order

Implementation should be incremental.

### Phase 1: establish generic QC contracts

- Introduce structured finding/policy concepts.
- Extract generic coverage subtraction from noise-specific ownership.
- Preserve current signal-noise output and thresholds.
- Keep final canonical APIs unchanged.

### Phase 2: migrate current directional evidence exclusions

- Represent validated directional polyC exclusions as QC findings/policy.
- Remove the corresponding untrusted per-trace coverage as well as the accepted
  candidate call.
- Preserve excluded candidates in diagnostic metadata.
- Keep trusted reverse/forward rescue behavior.

### Phase 3: add review-only specialized assessors

- Add read-edge findings.
- Add homopolymer/polyC quality findings that are not yet validated for automatic
  masking.
- Keep exploratory rules as `REVIEW` until independently validated.

### Phase 4: add additional artifacts

- Add dye-blob and other assessors only when detection criteria and validation
  data are available.

## 16. Validation requirements

Do not promote a new rule from `REVIEW` to `EXCLUDE` only because it fits the
same data used to design it.

For each specialized assessor, validate on an independent batch and measure at
least:

- sensitivity;
- specificity;
- false-positive rate;
- false-negative rate;
- positive/negative predictive value where prevalence makes them meaningful;
- effect on final trusted coverage;
- effect on final variant concordance;
- **false trusted-reference rate**: positions represented as trusted reference
  despite reference evidence being unreliable.

The last metric is especially important because the architecture is designed to
prevent an excluded candidate from accidentally becoming an implicit trusted
reference call.

Manual Sequencher ranges are useful operational reference annotations, but
important disagreements should be adjudicated with original trace evidence and
forward/reverse consistency rather than assuming every manual boundary is an
absolute biological truth.

## 17. Invariants

Implementation must preserve these invariants:

1. **Raw evidence is not destroyed by QC.**
2. **Aligned coverage is not automatically trusted coverage.**
3. **Removing an untrusted candidate must not silently create a trusted
   reference call.**
4. **Masking occurs per trace before trace merge.**
5. **A trusted mate trace may rescue a range excluded from another trace.**
6. **Trusted-vs-trusted disagreement uses the existing sample flag system.**
7. **Trace-local warnings are not blindly copied into final sample flags.**
8. **Coverage trimming is generic and can remove edge or internal ranges.**
9. **Assessors detect; policy decides; coverage code masks.**
10. **The complete Tracy QC layer can be disabled without making Tracy ETL
    depend on it.**
11. **Individual assessors can be enabled or disabled independently once their
    configuration is introduced.**
12. **Current signal-noise behavior must remain stable during the first
    modularization refactor.**

## 18. Non-goals

This design does not require the first implementation to:

- replace the canonical `Sample` or `Variant` models;
- introduce a new downstream output format;
- create a second sample-flagging framework;
- globally raise Tracy quality or SNR thresholds;
- automatically mask exploratory read-edge or homopolymer rules before
  validation;
- require both forward and reverse traces to be trusted at every position;
- discard excluded candidate observations from diagnostic output;
- implement dye-blob detection before validated criteria exist.

The first goal is a clean evidence boundary: **raw Tracy coverage and calls are
preserved, QC determines trusted per-trace evidence, generic masking can trim any
problematic span, and the existing merge/flagging system receives only the
interpretation state it needs.**
