import csv
import logging
import os
from collections import Counter, defaultdict
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
import pysam


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


def sanitize_name(name):
    """Replace spaces so labels are safe for downstream filenames and plotting."""
    return str(name).replace(" ", "_")


def get_qc_paths(arguments):
    """Return the QC output directory and standard output file paths."""
    qc_dir = os.path.join(arguments.output_dir, "QC_metrics")
    os.makedirs(qc_dir, exist_ok=True)

    return {
        "qc_dir": qc_dir,
        "bcs_tsv": os.path.join(qc_dir, "bcs.tsv"),
        "bcs_summary_tsv": os.path.join(qc_dir, "bcs_summary.tsv"),
        "reads_tsv": os.path.join(qc_dir, "reads.tsv"),
        "reads_summary_tsv": os.path.join(qc_dir, "reads_summary.tsv"),
        "reads_blocks_plot": os.path.join(qc_dir, "reads_and_blocks_summary.png"),
        "no_match_reads_fq": os.path.join(qc_dir, "no_match_reads.fq"),
        "ambiguous_reads_fq": os.path.join(qc_dir, "ambiguous_reads.fq"),
    }

def get_config_metadata(arguments, read_key):
    """
    Get metadata for barcodeList blocks from barcode_struct_r1 or barcode_struct_r2.
    Returns an empty list if the requested barcode_struct is not present.
    """
    if read_key not in arguments.config:
        return []

    metas = []
    for block in arguments.config[read_key]["blocks"]:
        if block["blocktype"] == "barcodeList":
            metas.append({
                "read_label": "R1" if read_key == "barcode_struct_r1" else "R2",
                "blockname": block.get("blockname", f"{read_key}_barcode"),
                "whitelist": block["sequence"],
            })
    return metas

def get_sanitized_metadata(arguments, config_key):
    """Load metadata and sanitize read/block labels once at the source."""
    metadata = get_config_metadata(arguments, config_key)
    sanitized = []

    for meta in metadata:
        meta_copy = dict(meta)
        meta_copy["read_label"] = sanitize_name(meta_copy["read_label"])
        meta_copy["blockname"] = sanitize_name(meta_copy["blockname"])
        sanitized.append(meta_copy)

    return sanitized


def handle_read_qc(read_qc, meta_list, counts, conflict_counts, block_summary_counts):
    """
    Update counts and conflict_counts for one read_qc entry.
    Also update one unique summary count per block.
    """
    for block_idx, (decoded_bc, status, conflicts) in enumerate(
        zip(read_qc["decoded_bcs"], read_qc["statuses"], read_qc["conflict_bcs"])
    ):
        meta = meta_list[block_idx]
        read_label = meta["read_label"]
        blockname = meta["blockname"]

        summary_key = (read_label, blockname)

        if status in ("exact", "corrected") and decoded_bc is not None:
            # If the read failed threshold, change corrected to below_threshold
            if status == "corrected" and read_qc.get('status') == '__BELOW_THRESHOLD__':
                status = "below_threshold"
            key = (read_label, blockname, decoded_bc)
            counts[key][status] += 1
            block_summary_counts[summary_key][status] += 1
            continue

        if status == "no_match":
            key = (read_label, blockname, "__NO_MATCH__")
            counts[key]["no_match"] += 1
            block_summary_counts[summary_key]["no_match"] += 1
            continue

        if status == "ambiguous":
            block_summary_counts[summary_key]["ambiguous"] += 1

            candidate_bcs, conflict_meta = parse_conflicts(conflicts)

            if candidate_bcs:
                for candidate_bc in candidate_bcs:
                    key = (read_label, blockname, candidate_bc)
                    counts[key]["ambiguous"] += 1

                    other_candidates = [bc for bc in candidate_bcs if bc != candidate_bc]
                    for other_bc in other_candidates:
                        conflict_counts[key][other_bc] += 1
            else:
                key = (read_label, blockname, "__AMBIGUOUS_UNIDENTIFIED__")
                counts[key]["ambiguous"] += 1
                for entry in conflict_meta:
                    conflict_counts[key][entry] += 1


def get_reads_header(arguments):
    """
    Build the TSV header for read-level QC output.

    The header includes the read-level status columns plus one block-name /
    block-status pair for each barcode block in the active read layout.
    """
    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_config_metadata(arguments, "barcode_struct_r2")

    if arguments.single_end_reads:
        max_blocks = len(r1_meta)
    else:
        max_blocks = len(r1_meta) + len(r2_meta)

    header = ["read_index", "read_name", "read_status"]
    for i in range(1, max_blocks + 1):
        header.extend([f"block_name_{i}", f"block_status_{i}"])
    return header


def get_read_row_single(arguments, read_idx, read_qc):
    """
    Convert one single-end read QC record into a TSV row.

    Returns a tuple of (row, collapsed_status), where collapsed_status is the
    final read-level status used for summary counting.
    """
    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    block_entries = []
    statuses = []

    for block_idx, status in enumerate(read_qc["statuses"]):
        blockname = r1_meta[block_idx]["blockname"]
        block_entries.extend([blockname, status])
        statuses.append(status)

    collapsed_status = collapse_read_status(statuses)
    if read_qc.get("status") == "__BELOW_THRESHOLD__" and collapsed_status in ["exact", "corrected"]:
        collapsed_status = "below_threshold"

    while len(block_entries) < 2 * len(r1_meta):
        block_entries.extend(["", ""])

    return [
        read_idx,
        read_qc.get("read_name"),
        collapsed_status,
        *block_entries,
    ], collapsed_status


def get_read_row_pair(arguments, read_idx, read_qc_pair):
    """
    Convert one paired-end read QC pair into a TSV row.

    Returns a tuple of (row, collapsed_status), where collapsed_status is the
    final read-level status used for summary counting.
    """
    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_config_metadata(arguments, "barcode_struct_r2")
    block_entries = []
    statuses = []

    for read_member_idx, read_qc in enumerate(read_qc_pair):
        if read_qc is None:
            continue

        meta_list = r1_meta if read_member_idx == 0 else r2_meta
        for block_idx, status in enumerate(read_qc["statuses"]):
            read_label = meta_list[block_idx]["read_label"]
            blockname = meta_list[block_idx]["blockname"]
            full_blockname = f"{read_label}_{blockname}"
            block_entries.extend([full_blockname, status])
            statuses.append(status)

    collapsed_status = collapse_read_status(statuses)
    if (
        any(
            read_qc is not None and read_qc.get("status") == "__BELOW_THRESHOLD__"
            for read_qc in read_qc_pair
        )
        and collapsed_status in ["exact", "corrected"]
    ):
        collapsed_status = "below_threshold"

    while len(block_entries) < 2 * (len(r1_meta) + len(r2_meta)):
        block_entries.extend(["", ""])

    return [
        read_idx,
        common_read_name_from_qc_pair(read_qc_pair),
        collapsed_status,
        *block_entries,
    ], collapsed_status


def write_read_summary(arguments, status_counts):
    """
        Write the read-level summary TSV in output_dir/QC_metrics.

        The output contains one row per collapsed read-level status: status, count.
    """

    paths = get_qc_paths(arguments)
    output_summary = paths["reads_summary_tsv"]
    with open(output_summary, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["status", "count"])
        for status in ["exact", "corrected", "below_threshold", "ambiguous", "no_match"]:
            writer.writerow([status, status_counts.get(status, 0)])


def write_qc_metrics_from_counts(arguments, counts, conflict_counts, block_summary_counts):
    """
        Write barcode-level QC TSV outputs in output_dir/QC_metrics.

    Outputs:
        - bcs.tsv
            One row per barcode label, with per-status counts and conflict summaries.
        - bcs_summary.tsv
            One row per barcode block, with aggregated status counts.

        The counts are expected to be accumulated incrementally during streaming
        preprocessing.
    """

    paths = get_qc_paths(arguments)
    output_file_bcs = paths["bcs_tsv"]
    output_file_summary = paths["bcs_summary_tsv"]
    status_columns = ["exact", "corrected", "below_threshold", "ambiguous", "no_match"]

    with open(output_file_bcs, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["read", "blockname", "barcode", *status_columns, "total", "conflict_bcs"])

        for meta in get_sanitized_metadata(arguments, "barcode_struct_r1") + get_sanitized_metadata(arguments, "barcode_struct_r2"):
            read_label = meta["read_label"]
            blockname = meta["blockname"]

            barcode_iter = chain(
                meta["whitelist"],
                ["__NO_MATCH__", "__AMBIGUOUS_UNIDENTIFIED__"],
            )

            for barcode_label in barcode_iter:
                key = (read_label, blockname, barcode_label)
                status_counter = counts.get(key)
                if status_counter is None:
                    status_values = [0 for _ in status_columns]
                    total = 0
                else:
                    status_values = [
                        status_counter.get(status, 0)
                        for status in status_columns
                    ]
                    total = sum(status_values)

                conflict_counter = conflict_counts.get(key)
                if conflict_counter:
                    conflict_summary = ";".join(
                        f"{entry}:{count}"
                        for entry, count in sorted(conflict_counter.items())
                    )
                else:
                    conflict_summary = ""

                writer.writerow([
                    read_label,
                    blockname,
                    barcode_label,
                    *status_values,
                    total,
                    conflict_summary,
                ])

    with open(output_file_summary, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["read", "blockname", *status_columns, "total"])

        for read_label, blockname in sorted(block_summary_counts.keys()):
            summary = block_summary_counts[(read_label, blockname)]
            status_values = [summary.get(status, 0) for status in status_columns]
            total = sum(status_values)
            writer.writerow([
                read_label,
                blockname,
                *status_values,
                total,
            ])

    log.info(f"Wrote barcode-level QC metrics to: {output_file_bcs}")
    log.info(f"Wrote block summary QC metrics to: {output_file_summary}")

def parse_conflicts(conflicts):
    """
    Split conflict information into parsable whitelist candidates and metadata entries.
    Metadata entries are things like "conflict_level:...".
    """
    if not conflicts:
        return [], []

    whitelist_candidates = []
    metadata_entries = []

    for entry in conflicts:
        # Since I implemented correct conflicting bcs returning, this should always jumpt to else. 
        if isinstance(entry, str) and entry.startswith("conflict_level:"):
            metadata_entries.append(entry)
        else:
            whitelist_candidates.append(entry)

    return whitelist_candidates, metadata_entries

def collapse_read_status(statuses):
    if "no_match" in statuses:
        return "no_match"
    if "ambiguous" in statuses:
        return "ambiguous"
    if "corrected" in statuses:
        return "corrected"
    return "exact"

def common_read_name(name1, name2):
    """
    Return the shared part of two paired-end read names.

    Handles simple cases like:
    read_1 / read_2
    sample_R1_001 / sample_R2_001
    """
    if not name1:
        return name2
    if not name2:
        return name1

    if name1 == name2:
        return name1

    # Keep characters that are identical at the same positions
    # Useful for names differing only by 1/2.
    if len(name1) == len(name2):
        shared = "".join(c1 for c1, c2 in zip(name1, name2) if c1 == c2)
        return shared.rstrip("_/-:. ")

    # Fallback: common prefix
    prefix = os.path.commonprefix([name1, name2])
    return prefix.rstrip("_/-:. ")


def common_read_name_from_qc_pair(read_qc_pair):
    read_names = [
        read_qc.get("read_name")
        for read_qc in read_qc_pair
        if read_qc is not None and read_qc.get("read_name")
    ]

    if not read_names:
        return ""

    if len(read_names) == 1:
        return read_names[0]

    return common_read_name(read_names[0], read_names[1])

STATUS_ORDER = ["exact", "corrected", "below_threshold", "ambiguous", "no_match"]
STATUS_COLOR_MAP = {
    "exact": "darkgreen",
    "corrected": "lightgreen",
    "below_threshold": "lightgrey",
    "ambiguous": "orange",
    "no_match": "grey",
}


def _draw_status_stack(ax, x, values_by_status, width=0.8, show_legend=True):
    """Draw a stacked barplot on an existing axis."""
    bottom = np.zeros(len(x), dtype=float)

    for status in STATUS_ORDER:
        values = np.asarray(values_by_status[status], dtype=float)
        ax.bar(
            x,
            values,
            width=width,
            bottom=bottom,
            color=STATUS_COLOR_MAP[status],
            label=status,
        )
        bottom += values

    if show_legend:
        ax.legend()


def _load_read_status_percentages(arguments):
    paths = get_qc_paths(arguments)
    summary_fpath = paths["reads_summary_tsv"]

    counts = {status: 0 for status in STATUS_ORDER}

    with open(summary_fpath, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status = row["status"]
            count = int(row["count"])
            if status in counts:
                counts[status] = count

    total = sum(counts.values())
    if total > 0:
        percentages = {status: 100 * counts[status] / total for status in STATUS_ORDER}
    else:
        percentages = {status: 0 for status in STATUS_ORDER}

    return counts, percentages


def _load_block_status_percentages(arguments):
    paths = get_qc_paths(arguments)
    input_file = paths["bcs_summary_tsv"]

    labels = []
    values_by_status = {status: [] for status in STATUS_ORDER}

    with open(input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            read_label = row["read"]
            blockname = row["blockname"]
            total = int(row["total"])

            labels.append(f"{read_label}_{blockname}")

            for status in STATUS_ORDER:
                count = int(row[status])
                if total > 0:
                    values_by_status[status].append(100 * count / total)
                else:
                    values_by_status[status].append(0)

    return labels, values_by_status


def _plot_reads_stacked_barplot_on_ax(ax, arguments, show_legend=True):
    _, percentages = _load_read_status_percentages(arguments)
    values_by_status = {status: [percentages[status]] for status in STATUS_ORDER}

    _draw_status_stack(
        ax=ax,
        x=["all_reads"],
        values_by_status=values_by_status,
        width=0.8,
        show_legend=show_legend,
    )

    ax.set_ylabel("Reads (%)")
    ax.set_ylim(0, 100)
    ax.set_title("QC metrics across all reads")


def _plot_blocks_stacked_barplot_on_ax(ax, arguments, show_legend=True):
    labels, values_by_status = _load_block_status_percentages(arguments)
    x = np.arange(len(labels))

    _draw_status_stack(
        ax=ax,
        x=x,
        values_by_status=values_by_status,
        width=0.8,
        show_legend=show_legend,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Reads (%)")
    ax.set_ylim(0, 100)
    ax.set_title("QC metrics per barcode block")


def _maybe_make_combined_reads_blocks_stacked_barplot(arguments):
    """
    Create a combined read/block stacked barplot once both summary files exist.

    The read-level plot is placed above the block-level plot.
    """
    paths = get_qc_paths(arguments)

    if not (
        os.path.exists(paths["reads_summary_tsv"])
        and os.path.exists(paths["bcs_summary_tsv"])
    ):
        return

    labels, _ = _load_block_status_percentages(arguments)
    fig_width = min(max(8, len(labels) * 1.2), 30)

    fig, (ax_reads, ax_blocks) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(fig_width, 10),
        gridspec_kw={"height_ratios": [1, 1.2]},
    )

    _plot_reads_stacked_barplot_on_ax(ax_reads, arguments, show_legend=True)
    _plot_blocks_stacked_barplot_on_ax(ax_blocks, arguments, show_legend=False)

    plt.tight_layout()
    plt.savefig(paths["reads_blocks_plot"], dpi=300)
    plt.close(fig)

    log.info(
        f"Wrote combined read/block stacked barplot to: {paths['reads_blocks_plot']}"
    )


def _barcode_rank_step(n_bcs):
    if n_bcs <= 20:
        return 1
    if n_bcs <= 100:
        return 5
    if n_bcs <= 500:
        return 25
    return max(1, n_bcs // 20)


def _barcode_values_by_status(rows, normalize=False):
    labels = [row["barcode"] for row in rows]

    values_by_status = {
        status: np.array([row[status] for row in rows], dtype=float)
        for status in STATUS_ORDER
    }

    if normalize:
        totals = np.zeros(len(labels), dtype=float)
        for values in values_by_status.values():
            totals += values

        # Avoid division by zero for completely empty barcode rows.
        totals_safe = np.where(totals == 0, 1, totals)

        for status in STATUS_ORDER:
            values_by_status[status] = values_by_status[status] / totals_safe * 100

    return labels, values_by_status


def _plot_stacked_barcodes_on_ax(
    ax,
    rows,
    read_label,
    blockname,
    normalize=False,
    show_legend=True,
    show_xlabel=True,
    show_rank_labels=True,
    show_barcode_axis=False,
):
    labels, values_by_status = _barcode_values_by_status(rows, normalize=normalize)
    n_bcs = len(labels)
    x = np.arange(n_bcs)

    _draw_status_stack(
        ax=ax,
        x=x,
        values_by_status=values_by_status,
        width=1.0,
        show_legend=show_legend,
    )

    ax.set_xlim(-0.5, n_bcs - 0.5)

    if normalize:
        ax.set_ylabel("Fraction of reads (%)")
        ax.set_ylim(0, 100)
        ax.set_title(f"QC metrics for {read_label}_{blockname} normalized to 100%")
    else:
        ax.set_ylabel("Read count")
        ax.set_title(f"QC metrics for {read_label}_{blockname}")

    if show_xlabel:
        ax.set_xlabel("Barcode rank (sorted by exact count)")

    ax.margins(x=0)

    rank_step = _barcode_rank_step(n_bcs)
    rank_ticks = np.arange(0, n_bcs, rank_step)
    rank_ticklabels = [str(i + 1) for i in rank_ticks]
    ax.set_xticks(rank_ticks)

    if show_rank_labels:
        ax.set_xticklabels(rank_ticklabels)
    else:
        ax.set_xticklabels([])

    # Optional barcode names as a top axis only when readable.
    # In combined plots this is only shown on the upper, non-normalized panel.
    if show_barcode_axis and n_bcs <= 50:
        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        ax_top.set_xticks(x)
        fontsize = 10 if n_bcs <= 20 else 8
        ax_top.set_xticklabels(labels, rotation=90, fontsize=fontsize)
        ax_top.set_xlabel("Barcode")


def make_stacked_barplots_bcs(arguments):
    """
    Create one combined stacked barplot per barcode block.

    Each output PNG contains:
    - the absolute read-count plot on top
    - the normalized 100% plot below

    Each bar corresponds to one barcode from the 'barcode' column in bcs.tsv.
    Bars are sorted by descending exact count.
    """
    paths = get_qc_paths(arguments)
    qc_dir = paths["qc_dir"]
    input_file = paths["bcs_tsv"]

    grouped_rows = defaultdict(list)

    with open(input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            barcode = row["barcode"]

            # Skip special rows.
            if barcode in {"__NO_MATCH__", "__AMBIGUOUS_UNIDENTIFIED__"}:
                continue

            read_label = row["read"]
            blockname = row["blockname"]

            grouped_rows[(read_label, blockname)].append({
                "barcode": barcode,
                "exact": int(row["exact"]),
                "corrected": int(row["corrected"]),
                "below_threshold": int(row["below_threshold"]),
                "ambiguous": int(row["ambiguous"]),
                "no_match": int(row["no_match"]),
            })

    for (read_label, blockname), rows in grouped_rows.items():
        if not rows:
            continue

        rows = sorted(rows, key=lambda row: row["exact"], reverse=True)
        n_bcs = len(rows)
        fig_width = min(max(8, n_bcs * 0.12), 30)

        fig, (ax_counts, ax_norm) = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(fig_width, 10),
            sharex=True,
            gridspec_kw={"height_ratios": [1, 1]},
        )

        _plot_stacked_barcodes_on_ax(
            ax=ax_counts,
            rows=rows,
            read_label=read_label,
            blockname=blockname,
            normalize=False,
            show_legend=True,
            show_xlabel=False,
            show_rank_labels=False,
            show_barcode_axis=True,
        )

        _plot_stacked_barcodes_on_ax(
            ax=ax_norm,
            rows=rows,
            read_label=read_label,
            blockname=blockname,
            normalize=True,
            show_legend=False,
            show_xlabel=True,
            show_rank_labels=True,
            show_barcode_axis=False,
        )

        plt.tight_layout()

        output_file = os.path.join(
            qc_dir,
            f"barcodes_{read_label}_{blockname}.png",
        )

        plt.savefig(output_file, dpi=300)
        plt.close(fig)

        log.info(f"Wrote combined barcode-level stacked barplot to: {output_file}")


def make_stacked_barplot_blocks(arguments):
    """
    Backward-compatible wrapper.

    The standalone block plot is no longer written.
    Instead, create the combined reads + blocks summary plot once both
    reads_summary.tsv and bcs_summary.tsv exist.
    """
    _maybe_make_combined_reads_blocks_stacked_barplot(arguments)


def make_stacked_barplot_reads(arguments):
    """
    Backward-compatible wrapper.

    The standalone read plot is no longer written.
    Instead, create the combined reads + blocks summary plot once both
    reads_summary.tsv and bcs_summary.tsv exist.
    """
    _maybe_make_combined_reads_blocks_stacked_barplot(arguments)
