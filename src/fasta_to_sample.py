"""FASTA -> Sample JSON variant caller and batch writer.

A single ``FastaToSample`` class imports an externally supplied per-HV-region
consensus FASTA (``HV1``/``HV2``/``HV3``, coordinate-aligned to rCRS), calls
variants with canonical rightmost-in-tandem-repeat normalization, and writes
the standard Sample JSON batch via :class:`src.core.batch.Batch`.

Reverse flow of ``src/generate_fasta.py``; reuses ``Sample`` / ``Batch`` /
``sample_to_dict`` and does not reimplement JSON output.  See
``docs/tools/fasta.md`` for the full algorithm, indel normalization, IUPAC/N
handling, and the reviewer-knowledge nomenclature conventions.
"""

from argparse import ArgumentParser, Namespace
from collections.abc import Callable
from pathlib import Path
from typing import NamedTuple

import dotenv
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from loguru import logger

from src.config import get_settings
from src.core.batch import Batch
from src.core.models import Sample, Tool, Variant, validate_position
from src.core.variants import pos_sort_key

dotenv.load_dotenv()

# Minimum tandem-repeat run lengths required for an indel to have an
# equivalent (score-preserving) rightmost placement.
_MIN_HOMOPOLYMER_RUN = 2
_MIN_DINUCLEOTIDE_RUN = 4
_DINUCLEOTIDE_PERIOD = 2
# Minimum 1-based position with a usable left flank for a singleton split.
_SINGLETON_MIN_POS = 2


# ---------------------------------------------------------------------------
# Alignment setup (the aligner is stateless and reused across calls/regions)
# ---------------------------------------------------------------------------


def _build_aligner() -> PairwiseAligner:
    """Build the global pairwise aligner used to place indels gap-aware.

    Scoring: match +2, mismatch -1, gap open -2, extend -1.  Global mode keeps
    the region slice and the consensus record end-to-end aligned so
    trailing/leading indels are recovered too.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    return aligner


_ALIGNER = _build_aligner()


# ---------------------------------------------------------------------------
# Tandem-repeat run detection for canonical indel normalization
# ---------------------------------------------------------------------------


def _extend_run(ref: str, pos: int, base: str, *, forward: bool) -> int:
    """Extend a homopolymer run of ``base`` from ``pos`` (1-based) in one direction.

    Returns the 1-based position of the last run base in that direction.
    """
    n = len(ref)
    cur = pos
    if forward:
        while cur + 1 <= n and ref[cur] == base:
            cur += 1
    else:
        while cur - 1 >= 1 and ref[cur - 2] == base:
            cur -= 1
    return cur


def _detect_homopolymer_run(ref: str, p_start: int, p_end: int) -> tuple[int, int] | None:
    """Detect a homopolymer (period-1) tandem run containing ``[p_start, p_end]``."""
    bases = ref[p_start - 1 : p_end]
    if len(set(bases)) != 1:
        return None
    base = bases[0]
    run_start = _extend_run(ref, p_start, base, forward=False)
    run_end = _extend_run(ref, p_end, base, forward=True)
    if run_end - run_start + 1 >= _MIN_HOMOPOLYMER_RUN:
        return run_start, run_end
    return None


def _detect_dinucleotide_run(ref: str, p_start: int, p_end: int) -> tuple[int, int] | None:
    """Detect a dinucleotide (period-2) tandem run containing ``[p_start, p_end]``."""
    length = p_end - p_start + 1
    n = len(ref)
    if length < _DINUCLEOTIDE_PERIOD or length % _DINUCLEOTIDE_PERIOD != 0:
        return None
    unit = ref[p_start - 1 : p_start + 1]

    def matches(idx: int) -> bool:
        return ref[idx - 1] == unit[(idx - p_start) % _DINUCLEOTIDE_PERIOD]

    if not all(ref[p_start - 1 + k] == unit[k % _DINUCLEOTIDE_PERIOD] for k in range(length)):
        return None
    run_start = p_start
    while run_start - 1 >= 1 and matches(run_start - 1):
        run_start -= 1
    run_end = p_end
    while run_end + 1 <= n and matches(run_end + 1):
        run_end += 1
    if run_end - run_start + 1 >= _MIN_DINUCLEOTIDE_RUN:
        return run_start, run_end
    return None


def _detect_tandem_run(ref: str, p_start: int, p_end: int) -> tuple[int, int, int]:
    """Return the maximal tandem-repeat run containing ``[p_start, p_end]``.

    Tries period 1 (homopolymer) then period 2 (dinucleotide).  Returns
    ``(run_start, run_end, period)`` as 1-based inclusive bounds; ``period`` is
    0 when no tandem run bounds the block.
    """
    homo = _detect_homopolymer_run(ref, p_start, p_end)
    if homo is not None:
        return homo[0], homo[1], 1
    di = _detect_dinucleotide_run(ref, p_start, p_end)
    if di is not None:
        return di[0], di[1], _DINUCLEOTIDE_PERIOD
    return p_start, p_end, 0


# ---------------------------------------------------------------------------
# Raw alignment-column walking
# ---------------------------------------------------------------------------


def _walk_alignment(ref_slice: str, query: str, start: int) -> list[tuple]:
    """Walk one region alignment and emit raw (un-normalized) events, in order:

    - ``("snp", pos, ref, seq)``
    - ``("del", p_start, p_end, deleted_ref)``  (1-based inclusive)
    - ``("ins", after_pos, inserted)``
    - ``("unread", pos)``  (query N; no variant)

    Gaps are read from the aligner's block coordinates, so the arbitrary column
    order inside a gap does not matter -- only deleted ref bases and inserted
    query bases are recorded, each anchored to its true rCRS position.
    """
    events: list[tuple] = []
    alignment = _ALIGNER.align(ref_slice, query)[0]
    ref_blocks, query_blocks = alignment.aligned

    prev_ref_end = 0
    prev_query_end = 0
    for (rs, re_), (qs, qe) in zip(ref_blocks, query_blocks, strict=True):
        # Gap before this aligned block.
        deleted = ref_slice[prev_ref_end:rs]
        inserted = query[prev_query_end:qs]
        if deleted:
            events.append(("del", int(start + prev_ref_end), int(start + rs - 1), deleted))
        if inserted:
            # Anchored after ref position ``start + rs - 1``; a leading gap
            # (rs == 0) anchors just before the region.
            events.append(("ins", int(start + rs - 1), inserted))

        # Aligned columns: matches / SNPs / unread.
        for i in range(int(re_ - rs)):
            rbase = ref_slice[rs + i]
            qbase = query[qs + i]
            gpos = int(start + rs + i)
            if rbase == qbase:
                continue
            if qbase == "N":
                events.append(("unread", gpos))  # N is unread; curated pipelines drop it too
            else:
                events.append(("snp", gpos, rbase, qbase))  # IUPAC ambiguity codes are callable

        prev_ref_end = re_
        prev_query_end = qe

    # Trailing gap after the last aligned block.
    deleted = ref_slice[prev_ref_end:]
    inserted = query[prev_query_end:]
    if deleted:
        events.append(("del", int(start + prev_ref_end), int(start + len(ref_slice) - 1), deleted))
    if inserted:
        events.append(("ins", int(start + len(ref_slice) - 1), inserted))

    return events


# ---------------------------------------------------------------------------
# Canonical normalization (biological, rightmost-in-tandem-repeat)
# ---------------------------------------------------------------------------


def _split_singleton_deletion(ref: str, p: int) -> list[tuple[int, str, str]] | None:
    """Singleton deletion -> SNP (p, X->Y) + del at the rightmost Y of the merged
    run.  Canonical convention for a base flanked by two identical homopolymer
    runs; see docs/tools/fasta.md (Nomenclature conventions).  Returns the two
    ``(pos, ref, seq)`` tuples, or None if the pattern does not apply.
    """
    if p < _SINGLETON_MIN_POS or p >= len(ref):
        return None
    x = ref[p - 1]
    left = _detect_homopolymer_run(ref, p - 1, p - 1)
    right = _detect_homopolymer_run(ref, p + 1, p + 1)
    if left is None or right is None:
        return None
    y = ref[p - 2]
    if y != ref[p] or x == y:
        return None
    return [(p, x, y), (right[1], y, "-")]


# A region convention sees the region's collected SNPs/deletions/insertion
# anchors (after biological normalization) and may rewrite them to the canonical
# reviewer-naming form.  Each returns the new ``(snps, dels, anchors)`` or None
# when the pattern does not apply.  Add an entry to ``REGION_CONVENTIONS`` as
# data when a convention is agreed; see docs/tools/fasta.md.

ConvResult = (
    tuple[
        list[tuple[int, str, str]],
        list[tuple[int, str, str]],
        list[tuple[int, str]],
    ]
    | None
)
ConventionFn = Callable[
    [list[tuple[int, str, str]], list[tuple[int, str, str]], list[tuple[int, str]], str],
    ConvResult,
]


class RegionConvention(NamedTuple):
    """One reviewer-knowledge naming convention (see docs/tools/fasta.md)."""

    name: str
    apply: ConventionFn


def _conv_16189_singleton(
    snps: list[tuple[int, str, str]],
    dels: list[tuple[int, str, str]],
    anchors: list[tuple[int, str]],
    ref: str,
) -> ConvResult:
    """del 16189 T -> SNP 16189 T>C + del at the rightmost C of the merged run.

    Only a true singleton deletion (neighbour 16188/16190 not also deleted)
    matches; a multi-base deletion spanning 16189 is left to the biological
    normalizer."""
    if (16189, "T", "-") not in dels or any(d[0] in (16188, 16190) for d in dels):
        return None
    res = _split_singleton_deletion(ref, 16189)
    if res is None:
        return None
    dels = [d for d in dels if d != (16189, "T", "-")]
    snps.append(res[0])  # (16189, T, C)
    dels.append(res[1])  # (16193, C, "-")
    return snps, dels, anchors


def _conv_polyc_run_shift(
    snps: list[tuple[int, str, str]],
    dels: list[tuple[int, str, str]],
    anchors: list[tuple[int, str]],
    ref: str,
) -> ConvResult:
    """polyC HV2 303-315: a C-run shift is named del + ins 315.1, not SNPs."""
    # 309 C>T + 310 T>C  ->  del 309 + ins 315.1
    if (309, "C", "T") in snps and (310, "T", "C") in snps:
        snps = [s for s in snps if s not in ((309, "C", "T"), (310, "T", "C"))]
        dels.append((309, ref[308], "-"))
        anchors.append((315, ref[314]))
        return snps, dels, anchors
    # 308 C>T + del 310 T  ->  del 308 + del 309 + ins 315.1
    if (308, "C", "T") in snps and (310, "T", "-") in dels and not any(d[0] in (309, 311) for d in dels):
        snps = [s for s in snps if s != (308, "C", "T")]
        dels = [d for d in dels if d != (310, "T", "-")]
        dels.extend([(308, ref[307], "-"), (309, ref[308], "-")])
        anchors.append((315, ref[314]))
        return snps, dels, anchors
    return None


def _conv_513_514_flank(
    snps: list[tuple[int, str, str]],
    dels: list[tuple[int, str, str]],
    anchors: list[tuple[int, str]],
    ref: str,
) -> ConvResult:
    """del 513 G + del 514 C -> SNP 513 G>A + del 523 A + del 524 C.

    The unique flank base (513 G) is not deleted; the length change moves to the
    poly-CA run end (523, 524), the canonical reviewer naming (round-trips).
    Only an isolated 513-514 block matches (512/515 not also deleted)."""
    if (513, "G", "-") in dels and (514, "C", "-") in dels and not any(d[0] in (512, 515) for d in dels):
        dels = [d for d in dels if d not in ((513, "G", "-"), (514, "C", "-"))]
        dels.extend([(523, ref[522], "-"), (524, ref[523], "-")])
        snps.append((513, "G", ref[514]))
        return snps, dels, anchors
    return None


# Reviewer-knowledge naming conventions (see docs/tools/fasta.md).  Disjoint
# sites; at most one fires per region.  Add an entry as data when agreed.
REGION_CONVENTIONS: tuple[RegionConvention, ...] = (
    RegionConvention("16189_singleton", _conv_16189_singleton),
    RegionConvention("polyc_run_shift", _conv_polyc_run_shift),
    RegionConvention("513_514_flank", _conv_513_514_flank),
)


def _apply_region_conventions(
    snps: list[tuple[int, str, str]],
    dels: list[tuple[int, str, str]],
    anchors: list[tuple[int, str]],
    ref: str,
) -> tuple[list[tuple[int, str, str]], list[tuple[int, str, str]], list[tuple[int, str]]]:
    """Apply reviewer-naming conventions in order (disjoint sites)."""
    for conv in REGION_CONVENTIONS:
        res = conv.apply(snps, dels, anchors, ref)
        if res is not None:
            snps, dels, anchors = res
    return snps, dels, anchors


def _normalize_deletion(
    ref: str, p_start: int, p_end: int
) -> tuple[list[tuple[int, str, str]], tuple[int, int, int] | None]:
    """Right-normalize a deletion block to the tandem-repeat run end.  Returns
    the ``(pos, ref, seq)`` tuples plus an optional SNP *slid zone*
    ``(lo, hi, shift)``: when the deletion is right-shifted within a tandem
    repeat, a SNP the leftmost alignment placed at ``lo..hi`` must re-anchor
    ``shift`` bases left to stay consistent with the canonical (rightmost) gap
    and round-trip the consensus (see docs/tools/fasta.md).  Reviewer-naming
    conventions are applied afterwards at the region level."""
    _run_start, run_end, period = _detect_tandem_run(ref, p_start, p_end)
    length = p_end - p_start + 1

    if period == 0:
        if length > 1:
            logger.warning(
                "Deletion at {}-{} cannot be bounded to a single tandem-repeat run; "
                "emitting best-effort call at leftmost position.",
                p_start,
                p_end,
            )
        return [(p, ref[p - 1], "-") for p in range(p_start, p_end + 1)], None

    new_start = run_end - length + 1
    # Right-shifted within the run -> SNPs in (old gap, new gap) re-anchor left.
    zone = (p_end + 1, new_start - 1, length) if new_start > p_start else None
    return [(p, ref[p - 1], "-") for p in range(new_start, new_start + length)], zone


def _normalize_dinucleotide_insertion(ref: str, after_pos: int, bases: str) -> list[tuple[int, str, str]] | None:
    """Right-normalize a period-2 repeat insertion as a block: slide to the run
    end and re-order bases to the run phase (see docs/tools/fasta.md).  Returns
    None if not a dinucleotide repeat insertion (caller falls back)."""
    if len(bases) < _DINUCLEOTIDE_PERIOD or len(bases) % _DINUCLEOTIDE_PERIOD != 0:
        return None
    n = len(ref)
    # Try the run just after the insertion point, then the run containing it.
    for p_start, p_end in ((after_pos + 1, after_pos + 2), (after_pos - 1, after_pos)):
        if not (p_start >= 1 and p_end <= n):
            continue
        run = _detect_dinucleotide_run(ref, p_start, p_end)
        if run is None:
            continue
        run_start, run_end = run
        x = ref[run_start - 1]
        y = ref[run_start]  # unit XY (ref[run_start], ref[run_start + 1])
        if x == y:
            continue  # homopolymer — handled by the period-1 path
        if set(bases) - {x, y}:
            continue
        if bases.count(x) != len(bases) // _DINUCLEOTIDE_PERIOD:
            continue
        # Phase of the run at run_end; re-ordered bases continue the pattern.
        phase_end = (run_end - run_start) % _DINUCLEOTIDE_PERIOD
        unit_at = (x, y)
        out = [unit_at[(phase_end + 1 + i) % _DINUCLEOTIDE_PERIOD] for i in range(len(bases))]
        if sorted(out) != sorted(bases):
            continue
        return [(run_end, "-", b) for b in out]
    return None


def _normalize_homopolymer_insertion(ref: str, after_pos: int, bases: str) -> list[tuple[int, str, str]]:
    """Right-normalize a period-1 insertion, per base: each base slides to the
    end of its flanking homopolymer run (polyC splits at T@310 -> 309.x / 315.x)."""
    results: list[tuple[int, str, str]] = []
    for b in bases:
        prev_in_run = ref[after_pos - 1] == b  # rCRS[after_pos] == b
        run_end = _extend_run(ref, after_pos, b, forward=True)
        next_in_run = run_end > after_pos  # rCRS[after_pos + 1] == b
        if prev_in_run or next_in_run:
            anchor = run_end if next_in_run else after_pos
        else:
            anchor = after_pos
            logger.warning(
                "Insertion of {} after {} is not in a tandem-repeat run; emitting best-effort call.",
                b,
                after_pos,
            )
        results.append((anchor, "-", b))
    return results


def _normalize_insertion(ref: str, after_pos: int, bases: str) -> list[tuple[int, str, str]]:
    """Normalize an insertion to the rightmost tandem-repeat position.  Tries
    period-2 (block) then period-1 (per base).  Returns ``(after_pos, "-", base)``
    tuples; the caller numbers sequential insertions at the same anchor."""
    di = _normalize_dinucleotide_insertion(ref, after_pos, bases)
    if di is not None:
        return di
    return _normalize_homopolymer_insertion(ref, after_pos, bases)


def _reanchor_snps(snps: list[tuple[int, str, str]], zones: list[tuple[int, int, int]]) -> list[tuple[int, str, str]]:
    """Re-anchor SNPs that fall in a deletion slid zone ``(lo, hi, shift)`` to the
    canonical rightmost-gap coordinate (``pos - shift``); see docs/tools/fasta.md."""
    out: list[tuple[int, str, str]] = []
    for pos, ref_base, seq_base in snps:
        adj = pos
        for lo, hi, shift in zones:
            if lo <= adj <= hi:
                adj -= shift
                break
        out.append((adj, ref_base, seq_base))
    return out


def _callable_intervals(start: int, end: int, unread: set[int]) -> list[list[int]]:
    """Callable sub-intervals of ``[start, end]`` excluding unread (N) positions.

    ``generate_sequence`` marks positions outside the intervals as ``N``, so
    excluding unread positions here preserves the FASTA's ``N`` bases in the
    reconstructed HV string (and matches the curated pipelines' "callable
    range" interval).  Returns ``[[start, end]]`` when nothing is unread.
    """
    if not unread:
        return [[start, end]]
    out: list[list[int]] = []
    cur = start
    for pos in range(start, end + 1):
        if pos in unread:
            if pos > cur:
                out.append([cur, pos - 1])
            cur = pos + 1
    if cur <= end:
        out.append([cur, end])
    return out


# ---------------------------------------------------------------------------
# Public class
# ---------------------------------------------------------------------------


class FastaToSample:
    """Call variants from per-HV-region consensus FASTA and write Sample JSON
    batches (see docs/tools/fasta.md).  Reference and regions load once per
    instance.  Use ``call_variants`` for one sample, ``process_batch`` /
    ``process_sample`` to also write the batch JSON layout."""

    def __init__(self, ref_path: str = "ref/rCRS.fasta") -> None:
        """Initialize the caller, loading the reference and region config.

        Args:
            ref_path: Path to the rCRS reference FASTA.
        """
        self.ref_path = ref_path
        self._ref_seq = str(next(SeqIO.parse(ref_path, "fasta")).seq).upper()
        self._regions = get_settings().regions.REGIONS

    def _process_region(
        self, start: int, end: int, query: str, stem: str
    ) -> tuple[list[Variant], list[tuple[int, str]], set[int]]:
        """Walk one region: return (variants, insertion anchors, unread positions).

        Insertion anchors are ``(after_pos, base)`` in alignment order; the caller
        numbers them with a global discovery index.  ``unread`` holds the rCRS
        positions of query ``N`` bases (no variant emitted).
        """
        ref_slice = self._ref_seq[start - 1 : end]
        snps: list[tuple[int, str, str]] = []  # (pos, ref, seq) before re-anchoring
        dels: list[tuple[int, str, str]] = []  # normalized (pos, ref, seq)
        anchors: list[tuple[int, str]] = []
        unread: set[int] = set()
        # Slid zones (lo, hi, shift) from right-normalized tandem-run deletions.
        zones: list[tuple[int, int, int]] = []

        for event in _walk_alignment(ref_slice, query, start):
            kind = event[0]
            if kind == "snp":
                _, pos, ref_base, seq_base = event
                snps.append((pos, ref_base, seq_base))
            elif kind == "del":
                _, p_start, p_end, _bases = event
                normalized, zone = _normalize_deletion(self._ref_seq, p_start, p_end)
                dels.extend(normalized)
                if zone is not None:
                    zones.append(zone)
            elif kind == "ins":
                _, after_pos, bases = event
                for anchor_after, _ref_base, seq_base in _normalize_insertion(self._ref_seq, after_pos, bases):
                    anchors.append((anchor_after, seq_base))
            elif kind == "unread":
                unread.add(event[1])  # N positions: excluded from the interval

        # Re-anchor SNPs in a deletion slid zone, then apply reviewer-naming
        # conventions (which may rewrite SNPs/deletions and add insertions).
        snps = _reanchor_snps(snps, zones)
        snps, dels, anchors = _apply_region_conventions(snps, dels, anchors, self._ref_seq)
        variants = [
            Variant(pos=validate_position(p), ref=r, seq=s, quality=[], files=[stem]) for p, r, s in dels + snps
        ]
        return variants, anchors, unread

    def call_variants(self, fasta_path: str | Path) -> Sample:
        """Call variants from one per-sample FASTA -> ``Sample`` (source_tool=FASTA).

        Records matched by id (HV1/HV2/HV3); missing regions omitted.  Query
        uppercased and ids normalized so mixed/lowercase external consensus imports.
        """
        fasta_path = Path(fasta_path)
        records = {rec.id.strip().upper(): str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")}

        stem = fasta_path.stem
        raw_variants: list[Variant] = []
        # (discovery_index, after_pos, base) for insertions; numbered sequentially
        # per anchor base after all regions are walked.
        insertion_anchors: list[tuple[int, int, str]] = []
        intervals: dict[str, list[list[int]]] = {}

        discovery = 0
        for region, (start, end) in self._regions.items():
            query = records.get(region)
            if query is None:
                continue
            region_variants, anchors, unread = self._process_region(start, end, query, stem)
            raw_variants.extend(region_variants)
            for anchor_after, seq_base in anchors:
                insertion_anchors.append((discovery, anchor_after, seq_base))
                discovery += 1
            intervals[region] = _callable_intervals(start, end, unread)

        # Number sequential insertions at the same anchor as "<anchor>.<i>",
        # sorted by (anchor, discovery order) so insertions sharing an anchor
        # receive .1, .2, ... in alignment order.
        insertion_anchors.sort(key=lambda item: (item[1], item[0]))
        insertions: list[Variant] = []
        suffix_counts: dict[int, int] = {}
        for _discovery, anchor_after, seq_base in insertion_anchors:
            idx = suffix_counts.get(anchor_after, 0) + 1
            suffix_counts[anchor_after] = idx
            insertions.append(
                Variant(pos=validate_position(f"{anchor_after}.{idx}"), ref="-", seq=seq_base, quality=[], files=[stem])
            )

        variants = raw_variants + insertions
        variants.sort(key=lambda v: pos_sort_key(v.pos))

        return Sample(
            sample_id=stem,
            variants=variants,
            source_tool=Tool.FASTA,
            intervals=intervals or None,
            sample_flags=[],
            information={"source_tool": "fasta"},
        )

    def process_batch(self, input_dir: str | Path, output_dir: str | Path, batch_id: str) -> Path | None:
        """Process every ``.fasta`` file in ``input_dir`` and write batch JSON.

        Returns the path to the written ``statistic_fullbatch.json``, or
        ``None`` when no samples were produced.
        """
        input_dir = Path(input_dir)
        files = [p for p in sorted(input_dir.iterdir()) if p.is_file() and p.suffix == ".fasta"]
        samples = [self._call_and_tag(f, batch_id) for f in files]
        samples = [s for s in samples if s is not None]
        if not samples:
            logger.warning("No samples produced from {}", input_dir)
            return None
        json_dir = Path(output_dir) / batch_id / "json"
        return Batch(samples, self.ref_path).write(json_dir, batch_id, nest_batch_id=False)

    def process_sample(self, fasta_path: str | Path, output_dir: str | Path, batch_id: str | None = None) -> Path:
        """Process a single FASTA file and write batch JSON.

        ``batch_id`` defaults to the FASTA filename stem.  Returns the path to
        the written ``statistic_fullbatch.json``.
        """
        path = Path(fasta_path)
        batch_id = batch_id or path.stem
        sample = self._call_and_tag(path, batch_id)
        if sample is None:
            msg = f"No variants produced for {path}"
            raise ValueError(msg)
        json_dir = Path(output_dir) / batch_id / "json"
        return Batch([sample], self.ref_path).write(json_dir, batch_id, nest_batch_id=False)

    def _call_and_tag(self, fasta_path: Path, batch_id: str) -> Sample | None:
        """Call variants for one FASTA, tag it with ``batch_id``, log the result."""
        try:
            sample = self.call_variants(fasta_path)
            sample.batch_id = batch_id
            logger.info("Called variants for {} ({} variants)", fasta_path.name, len(sample.variants))
        except (ValueError, KeyError, OSError) as exc:
            logger.error("Processing FASTA {}: {}", fasta_path, exc)
            return None
        else:
            return sample


def parse_args() -> Namespace:
    """Parse command line arguments."""
    parser = ArgumentParser(description="Call variants from per-sample FASTA and write Sample JSON")
    parser.add_argument("--input-dir", default=None, help="Directory of per-sample .fasta files")
    parser.add_argument("--output-dir", required=True, help="Base output directory (e.g. results/tools/fasta)")
    parser.add_argument("--batch-id", default=None, help="Batch identifier")
    parser.add_argument("--ref-path", default="ref/rCRS.fasta", help="Path to reference FASTA")
    return parser.parse_args()


def main() -> None:
    """Main entry point: process every per-sample FASTA in an input directory."""
    args = parse_args()
    converter = FastaToSample(ref_path=args.ref_path)

    if not args.input_dir:
        missing_input = "--input-dir is required"
        raise SystemExit(missing_input)
    if not args.batch_id:
        missing_batch = "--batch-id is required"
        raise SystemExit(missing_batch)

    out_path = converter.process_batch(args.input_dir, args.output_dir, args.batch_id)
    if out_path is not None:
        logger.success("Wrote batch JSON to {}", out_path)


if __name__ == "__main__":
    main()
