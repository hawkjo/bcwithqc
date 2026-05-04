import logging
import os
import pysam
import itertools

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import misc
from .bc_aligner import CustomBCAligner

log = logging.getLogger(__name__)

def _candidate_boundary_shifts(blocks, raw_bounds, seq_len):
    """
    Generate nearby boundary candidates for barcode rescue.

    The initial aligner can place block boundaries one base too early/late when an
    indel occurs in a wildcard block.  We therefore try small shifts around each
    inferred boundary, using the maxerrors of the adjacent blocks to limit the
    search space.  The unshifted boundaries are yielded first.
    """
    raw_bounds = [int(b) for b in raw_bounds]

    if not raw_bounds:
        return

    # Always try the original segmentation first.
    seen = {tuple(raw_bounds)}
    yield raw_bounds

    windows = []
    for i in range(len(raw_bounds)):
        left_err = blocks[i].get("maxerrors", 0) if i < len(blocks) else 0
        right_err = blocks[i + 1].get("maxerrors", 0) if i + 1 < len(blocks) else 0
        # One-base shifts catch the common case even for random/UMI blocks where
        # maxerrors can be 0.  The cap prevents combinatorial explosions for more
        # complex barcode structures.
        window = max(1, left_err, right_err)
        window = min(window, 2)
        windows.append(range(-window, window + 1))

    candidates = []
    for shifts in itertools.product(*windows):
        if all(shift == 0 for shift in shifts):
            continue

        bounds = [bound + shift for bound, shift in zip(raw_bounds, shifts)]
        bounds_tup = tuple(bounds)
        if bounds_tup in seen:
            continue
        seen.add(bounds_tup)

        if bounds[0] < 0 or bounds[-1] > seq_len:
            continue
        if any(left > right for left, right in zip(bounds, bounds[1:])):
            continue

        # Prefer minimal rescue shifts, then deterministic left-to-right order.
        candidates.append((sum(abs(s) for s in shifts), shifts, bounds))

    for _, _, bounds in sorted(candidates):
        yield bounds


def _pieces_from_bounds(seq, bounds):
    starts = [0] + list(bounds[:-1])
    return [str(seq[start:end]) for start, end in zip(starts, bounds)]


def _parse_pieces_against_config(blocks, raw_pieces, decoders, commonseqs):
    """
    Decode barcode blocks and validate constant regions for one segmentation.

    Returns:
      - parsed: parsed fields if the segmentation is valid, otherwise None
      - qc_fields: barcode decoding state for QC reporting
      - failure_debug: structured information for dropped-read logging
    """
    raw_bcs = [
        raw_piece.upper().replace('N', 'A')
        for raw_piece, block in zip(raw_pieces, blocks)
        if block["blocktype"] == "barcodeList"
    ]
    decode_results = [
        decoder.decode_with_status(raw_bc)
        for raw_bc, decoder in zip(raw_bcs, decoders)
    ]

    if decode_results:
        bcs, statuses, overlapping_bcs = map(list, zip(*decode_results))
    else:
        bcs, statuses, overlapping_bcs = [], [], []

    qc_fields = {
        "raw_bcs": raw_bcs,
        "decoded_bcs": bcs,
        "statuses": statuses,
        "conflict_bcs": overlapping_bcs,
    }

    bc_idx = commonseq_idx = 0
    corrected_pieces = []

    for i, block in enumerate(blocks):
        if block["blocktype"] == "barcodeList":
            if bcs[bc_idx] is None:
                failure_debug = {
                    "reason": "barcode",
                    "raw_pieces": raw_pieces,
                    "raw_bcs": raw_bcs,
                    "statuses": statuses,
                    "decoded_bcs": bcs,
                }
                return None, qc_fields, failure_debug
            corrected_pieces.append(bcs[bc_idx])
            bc_idx += 1

        elif block["blocktype"] == "constantRegion":
            if commonseqs[commonseq_idx] is None:
                failure_debug = {
                    "reason": "barcode",
                    "raw_pieces": raw_pieces,
                    "raw_bcs": raw_bcs,
                    "statuses": statuses,
                    "decoded_bcs": bcs,
                }
                return None, qc_fields, failure_debug

            observed_const = raw_pieces[i].upper().replace("N", "A")
            expected_const = commonseqs[commonseq_idx]

            distfun = misc.DistanceThresh("levenshtein", block["maxerrors"])
            dist = distfun(observed_const, expected_const)
            if dist is False:
                failure_debug = {
                    "reason": "constant",
                    "block_index": i,
                    "observed_const": observed_const,
                    "expected_const": expected_const,
                    "raw_pieces": raw_pieces,
                    "raw_bcs": raw_bcs,
                    "statuses": statuses,
                    "decoded_bcs": bcs,
                }
                return None, qc_fields, failure_debug

            corrected_pieces.append(expected_const)
            commonseq_idx += 1

        else:
            corrected_pieces.append('N' * block["length"])

    parsed = {
        **qc_fields,
        "corrected_pieces": corrected_pieces,
    }
    return parsed, qc_fields, None


def process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders):
    """
    Find barcodes etc in bc_rec.

    The fast path uses the boundaries from CustomBCAligner.  If that segmentation
    fails barcode decoding or constant-region validation, try nearby boundary
    shifts before dropping the read.  This rescues valid reads where an indel in a
    wildcard block shifted the true barcode/constant boundary.
    """

    scores_pieces_end_bounds = [al.find_norm_score_pieces_and_boundaries(bc_rec.seq, return_seq=True) for al in aligners]
    scores_pieces_end_pos = [(s, p, bounds[-1], sq) for s, p, bounds, sq in scores_pieces_end_bounds]

    raw_score, raw_pieces, raw_end_pos, raw_seq = max(scores_pieces_end_pos)
    raw_bounds = next(
        bounds for s, p, bounds, sq in scores_pieces_end_bounds
        if s == raw_score and p == raw_pieces
    )
    raw_bounds = [int(bound) for bound in raw_bounds]
    raw_end_pos = int(raw_end_pos)

    best_aligner = next(al for al, (s, p, e, sq) in zip(aligners, scores_pieces_end_pos) if s == raw_score)
    commonseqs = [best_aligner.prefixes[i] for i, block in enumerate(blocks) if block["blocktype"] == "constantRegion"]

    parsed = None
    failure_qc_fields = None
    failure_debug = None
    failure_bounds = raw_bounds
    selected_bounds = raw_bounds
    selected_pieces = raw_pieces

    # First try the aligner's original segmentation, then nearby rescue shifts.
    for candidate_bounds in _candidate_boundary_shifts(blocks, raw_bounds, len(raw_seq)):
        candidate_pieces = raw_pieces if candidate_bounds == raw_bounds else _pieces_from_bounds(raw_seq, candidate_bounds)
        candidate_parsed, candidate_qc_fields, candidate_failure_debug = _parse_pieces_against_config(
            blocks,
            candidate_pieces,
            decoders,
            commonseqs,
        )

        if failure_qc_fields is None:
            failure_qc_fields = candidate_qc_fields
            failure_debug = candidate_failure_debug
            failure_bounds = candidate_bounds

        if candidate_parsed is not None:
            parsed = candidate_parsed
            selected_bounds = candidate_bounds
            selected_pieces = candidate_pieces
            break

    if parsed is None:
        bc_qc = {
            "read_name": str(bc_rec.id),
            **failure_qc_fields,
        }

        # Keep the original dropped-read diagnostics, but emit them only after
        # all boundary-rescue candidates have failed.  This prevents rescued
        # reads from being logged as dropped reads.
        if failure_debug and failure_debug.get("reason") == "constant":
            log.warning(
                "Dropped read %s: constant failed block=%s observed=%s expected=%s raw_bounds=%s raw_pieces=%s",
                bc_rec.id,
                failure_debug["block_index"],
                failure_debug["observed_const"],
                failure_debug["expected_const"],
                failure_bounds,
                failure_debug["raw_pieces"],
            )
        else:
            log.warning(
                "Dropped read %s: score=%s raw_bounds=%s raw_pieces=%s raw_bcs=%s statuses=%s decoded=%s",
                bc_rec.id,
                raw_score,
                failure_bounds,
                failure_debug["raw_pieces"] if failure_debug else raw_pieces,
                failure_qc_fields.get("raw_bcs", []),
                failure_qc_fields.get("statuses", []),
                failure_qc_fields.get("decoded_bcs", []),
            )

        return raw_score, None, None, bc_qc

    raw_pieces = selected_pieces
    raw_bounds = selected_bounds
    raw_bcs = parsed["raw_bcs"]
    bcs = parsed["decoded_bcs"]
    statuses = parsed["statuses"]
    overlapping_bcs = parsed["conflict_bcs"]
    corrected_pieces = parsed["corrected_pieces"]
    raw_end_pos = raw_bounds[-1]

    # QC Metrics
    bc_qc = {
        "read_name": str(bc_rec.id),
        "raw_bcs": raw_bcs,
        "decoded_bcs": bcs,
        "statuses": statuses,
        "conflict_bcs": overlapping_bcs,
    }

    new_aligner = CustomBCAligner(*corrected_pieces)
    new_score, new_pieces = new_aligner.find_norm_score_and_pieces(raw_seq)

    bam_bcs = []
    bam_raw_bcs = []
    bam_umis = []
    bam_commonseqs = []

    bc_idx = commonseq_idx = 0
    for i, block in enumerate(blocks):
        if block["blockfunction"] != "discard":
            if block["blocktype"] == "barcodeList":
                bam_bcs.append(bcs[bc_idx])
                bam_raw_bcs.append(raw_pieces[i])
                bc_idx += 1
            elif block["blocktype"] == "constantRegion":
                bam_commonseqs.append(commonseqs[commonseq_idx])
                commonseq_idx += 1
            else:
                if block["blockfunction"] == "UMI":
                    bam_umis.append(new_pieces[i])
                else:
                    bam_bcs.append(new_pieces[i])
                    bam_raw_bcs.append(new_pieces[i])

    tags = []
    # Cell barcode
    tags.append(('CB', '.'.join(bam_bcs)))
    tags.append(('CR', '.'.join(bam_raw_bcs)))
    # Filler sequences
    tags.append(('FL', sum(len(seq) for seq in bam_commonseqs)))
    # And raw UMI
    tags.append(("UR", ".".join(bam_umis)))

    starts = [0] + list(raw_bounds[:-1])
    quals = bc_rec.letter_annotations["phred_quality"]

    keep_seq_parts = []
    keep_qual_parts = []

    # Keep selected defined blocks
    for block, piece, start, end in zip(blocks, raw_pieces, starts, raw_bounds):
        if block.get("keepblock", False):
            keep_seq_parts.append(piece)
            keep_qual_parts.extend(quals[start:end])

    # Keep undefined / trailing sequence if requested
    if keep_nonbarcode:
        trailing_rec = bc_rec[raw_end_pos:]
        keep_seq_parts.append(str(trailing_rec.seq))
        keep_qual_parts.extend(trailing_rec.letter_annotations["phred_quality"])

    new_seq = "".join(keep_seq_parts)

    new_rec = SeqRecord(
        Seq(new_seq),
        id=bc_rec.id,
        description=bc_rec.description,
    )
    new_rec.letter_annotations["phred_quality"] = keep_qual_parts

    return raw_score, new_rec, tags, bc_qc
