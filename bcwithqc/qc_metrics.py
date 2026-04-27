import csv
import logging
import os
from collections import Counter, defaultdict

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
        "bcs_tsv": os.path.join(qc_dir, "QC_metrics_bcs.tsv"),
        "bcs_summary_tsv": os.path.join(qc_dir, "QC_metrics_bcs_summary.tsv"),
        "reads_tsv": os.path.join(qc_dir, "QC_metrics_reads.tsv"),
        "reads_summary_tsv": os.path.join(qc_dir, "QC_metrics_reads_summary.tsv"),
        "blocks_plot": os.path.join(qc_dir, "QC_stacked_barplot_blocks.png"),
        "reads_plot": os.path.join(qc_dir, "QC_metrics_reads_stacked_barplot.png"),
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

def iter_read_qc_rows(arguments, read_qcs):
    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_config_metadata(arguments, "barcode_struct_r2")

    if arguments.single_end_reads:
        for read_idx, read_qc in enumerate(read_qcs, start=1):
            for block_idx, (raw_bc, decoded_bc, status, conflict_bcs) in enumerate(
                zip(
                    read_qc["raw_bcs"],
                    read_qc["decoded_bcs"],
                    read_qc["statuses"],
                    read_qc["conflict_bcs"],
                )
            ):
                meta = r1_meta[block_idx]
                read_name = read_qc.get("read_name", meta["read_label"])
                yield {
                    "read_index": read_idx,
                    "read_label": read_name,
                    "block_index": block_idx,
                    "blockname": meta["blockname"],
                    "raw_bc": raw_bc,
                    "decoded_bc": decoded_bc if decoded_bc is not None else "",
                    "status": status,
                    "conflict_bcs": ";".join(conflict_bcs) if conflict_bcs else "",
                }
    else:
        for read_idx, read_qc_pair in enumerate(read_qcs, start=1):
            for read_member_idx, read_qc in enumerate(read_qc_pair):
                if read_qc is None:
                    continue
                meta_list = r1_meta if read_member_idx == 0 else r2_meta
                for block_idx, (raw_bc, decoded_bc, status, conflict_bcs) in enumerate(
                    zip(
                        read_qc["raw_bcs"],
                        read_qc["decoded_bcs"],
                        read_qc["statuses"],
                        read_qc["conflict_bcs"],
                    )
                ):
                    meta = meta_list[block_idx]
                    read_name = read_qc.get("read_name", meta["read_label"])
                    yield {
                        "read_index": read_idx,
                        "read_label": read_name,
                        "block_index": block_idx,
                        "blockname": meta["blockname"],
                        "raw_bc": raw_bc,
                        "decoded_bc": decoded_bc if decoded_bc is not None else "",
                        "status": status,
                        "conflict_bcs": ";".join(conflict_bcs) if conflict_bcs else "",
                    }
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
    
def generate_qc_metrics(arguments, read_qcs):
    """
    Generate barcode-level and block-level QC TSV files in output_dir/QC_metrics.

    Outputs:
    - QC_metrics_bcs.tsv
    - QC_metrics_bcs_summary.tsv
    """
    log.info(f"Collected {len(read_qcs):,d} read QC entries")

    paths = get_qc_paths(arguments)
    output_file_bcs = paths["bcs_tsv"]
    output_file_summary = paths["bcs_summary_tsv"]
    status_columns = ["exact", "corrected", "ambiguous", "no_match"]

    r1_meta = get_sanitized_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_sanitized_metadata(arguments, "barcode_struct_r2")

    counts = defaultdict(Counter)
    conflict_counts = defaultdict(Counter)
    block_summary_counts = defaultdict(Counter)
    row_keys = []

    # Pre-initialize whitelist rows plus special rows.
    for meta in r1_meta + r2_meta:
        for whitelist_bc in meta["whitelist"]:
            key = (meta["read_label"], meta["blockname"], whitelist_bc)
            row_keys.append(key)
            counts[key]

        for special in ["__NO_MATCH__", "__AMBIGUOUS_UNIDENTIFIED__"]:
            key = (meta["read_label"], meta["blockname"], special)
            row_keys.append(key)
            counts[key]

        block_summary_counts[(meta["read_label"], meta["blockname"])]

    if arguments.single_end_reads:
        for read_qc in read_qcs:
            handle_read_qc(read_qc, r1_meta, counts, conflict_counts, block_summary_counts)
    else:
        for read_qc_pair in read_qcs:
            for read_idx, read_qc in enumerate(read_qc_pair):
                if read_qc is None:
                    continue
                meta_list = r1_meta if read_idx == 0 else r2_meta
                handle_read_qc(read_qc, meta_list, counts, conflict_counts, block_summary_counts)

    with open(output_file_bcs, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["read", "blockname", "barcode", *status_columns, "total", "conflict_bcs"])

        for key in row_keys:
            read_label, blockname, barcode_label = key
            status_counter = counts[key]
            total = sum(status_counter.values())

            conflict_counter = conflict_counts.get(key, Counter())
            if conflict_counter:
                conflict_summary = ";".join(
                    f"{entry}:{count}" for entry, count in sorted(conflict_counter.items())
                )
            else:
                conflict_summary = ""

            writer.writerow([
                read_label,
                blockname,
                barcode_label,
                *[status_counter.get(status, 0) for status in status_columns],
                total,
                conflict_summary,
            ])

    with open(output_file_summary, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["read", "blockname", *status_columns, "total"])

        for read_label, blockname in sorted(block_summary_counts.keys()):
            summary = block_summary_counts[(read_label, blockname)]
            exact = summary.get("exact", 0)
            corrected = summary.get("corrected", 0)
            ambiguous = summary.get("ambiguous", 0)
            no_match = summary.get("no_match", 0)
            total = exact + corrected + ambiguous + no_match

            writer.writerow([
                read_label,
                blockname,
                exact,
                corrected,
                ambiguous,
                no_match,
                total,
            ])

    log.info(f"Wrote barcode-level QC metrics to: {output_file_bcs}")
    log.info(f"Wrote block summary QC metrics to: {output_file_summary}")

def generate_qc_metrics_reads(arguments, read_qcs):
    """
    Generate read-level QC outputs in output_dir/QC_metrics.

    Outputs:
    - QC_metrics_reads.tsv
      One row per read (single-end) or read pair (paired-end), with:
        read_index, read_label, status, blockname_1, blockstatus_1, blockname_2, blockstatus_2, ...

    - QC_metrics_reads_summary.tsv
      One row per collapsed read-level status:
        status, count

    Uses collapse_read_status(statuses) to assign one final status per read/read-pair.
    """
    paths = get_qc_paths(arguments)
    output_reads = paths["reads_tsv"]
    output_summary = paths["reads_summary_tsv"]

    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_config_metadata(arguments, "barcode_struct_r2")

    read_rows = []
    summary_statuses = []

    if arguments.single_end_reads:
        max_blocks = len(r1_meta)

        for read_idx, read_qc in enumerate(read_qcs, start=1):
            block_entries = []
            statuses = []

            for block_idx, status in enumerate(read_qc["statuses"]):
                blockname = r1_meta[block_idx]["blockname"]
                block_entries.extend([blockname, status])
                statuses.append(status)

            collapsed_status = collapse_read_status(statuses)
            summary_statuses.append(collapsed_status)

            read_rows.append({
                "read_index": read_idx,
                "read_label": read_qc.get("read_name"),
                "status": collapsed_status,
                "block_entries": block_entries,
            })

    else:
        max_blocks = len(r1_meta) + len(r2_meta)

        for read_idx, read_qc_pair in enumerate(read_qcs, start=1):
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
            summary_statuses.append(collapsed_status)

            read_rows.append({
                "read_index": read_idx,
                "read_label": common_read_name_from_qc_pair(read_qc_pair),
                "status": collapsed_status,
                "block_entries": block_entries,
            })

    header = ["read_index", "read_name", "read_status"]
    for i in range(1, max_blocks + 1):
        header.extend([f"block_name_{i}", f"block_status_{i}"])

    with open(output_reads, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(header)

        for row in read_rows:
            padded_entries = list(row["block_entries"])
            while len(padded_entries) < 2 * max_blocks:
                padded_entries.extend(["", ""])

            writer.writerow([
                row["read_index"],
                row["read_label"],
                row["status"],
                *padded_entries,
            ])

    status_counts = Counter(summary_statuses)
    with open(output_summary, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t")
        writer.writerow(["status", "count"])
        for status in ["exact", "corrected", "ambiguous", "no_match"]:
            writer.writerow([status, status_counts.get(status, 0)])

    log.info(f"Wrote read-level QC metrics to: {output_reads}")
    log.info(f"Wrote read-level QC summary to: {output_summary}")

def make_stacked_barplot_blocks(arguments):
    """
    Read QC_metrics_bcs_summary.tsv and create one stacked barplot with one bar per block.

    Dark green  = exact
    Light green = corrected
    Orange      = ambiguous
    Grey        = no_match

    Y-axis is shown as percentage.
    """
    paths = get_qc_paths(arguments)
    input_file = paths["bcs_summary_tsv"]
    output_file = paths["blocks_plot"]

    labels = []
    exact_pct = []
    corrected_pct = []
    ambiguous_pct = []
    no_match_pct = []

    with open(input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            read_label = row["read"]
            blockname = row["blockname"]
            exact = int(row["exact"])
            corrected = int(row["corrected"])
            ambiguous = int(row["ambiguous"])
            no_match = int(row["no_match"])
            total = int(row["total"])

            labels.append(f"{read_label}_{blockname}")

            if total > 0:
                exact_pct.append(100 * exact / total)
                corrected_pct.append(100 * corrected / total)
                ambiguous_pct.append(100 * ambiguous / total)
                no_match_pct.append(100 * no_match / total)
            else:
                exact_pct.append(0)
                corrected_pct.append(0)
                ambiguous_pct.append(0)
                no_match_pct.append(0)

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(labels))

    ax.bar(x, exact_pct, label="exact", color="darkgreen")
    ax.bar(x, corrected_pct, bottom=exact_pct, label="corrected", color="lightgreen")
    ax.bar(
        x,
        ambiguous_pct,
        bottom=[e + c for e, c in zip(exact_pct, corrected_pct)],
        label="ambiguous",
        color="orange",
    )
    ax.bar(
        x,
        no_match_pct,
        bottom=[e + c + a for e, c, a in zip(exact_pct, corrected_pct, ambiguous_pct)],
        label="no_match",
        color="grey",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Reads (%)")
    ax.set_ylim(0, 100)
    ax.set_title("QC metrics per barcode block")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    log.info(f"Wrote block-level stacked barplot to: {output_file}")


def make_stacked_barplots_bcs(arguments):
    """
    Create one stacked barplot per barcode block, with one bar per whitelist barcode.

    Each bar corresponds to one barcode from the 'barcode' column in QC_metrics_bcs.tsv.
    Bars are sorted by descending exact count.

    Colors:
    - dark green: exact
    - light green: corrected
    - orange: ambiguous
    - grey: no_match

    The x-axis always shows barcode rank (1..n).
    Barcode names are only shown when there are few enough to remain readable.

    Plots are written into the QC_metrics subdirectory.
    """
    qc_dir = os.path.join(arguments.output_dir, "QC_metrics")
    os.makedirs(qc_dir, exist_ok=True)

    input_file = os.path.join(qc_dir, "QC_metrics_bcs.tsv")

    grouped_rows = defaultdict(list)

    with open(input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            barcode = row["barcode"]

            # Skip special rows
            if barcode in {"__NO_MATCH__", "__AMBIGUOUS_UNIDENTIFIED__"}:
                continue

            read_label = row["read"]
            blockname = row["blockname"]

            grouped_rows[(read_label, blockname)].append({
                "barcode": barcode,
                "exact": int(row["exact"]),
                "corrected": int(row["corrected"]),
                "ambiguous": int(row["ambiguous"]),
                "no_match": int(row["no_match"]),
            })

    for (read_label, blockname), rows in grouped_rows.items():
        if not rows:
            continue

        # Sort by highest exact count to lowest
        rows = sorted(rows, key=lambda x: x["exact"], reverse=True)

        labels = [r["barcode"] for r in rows]
        exact_vals = [r["exact"] for r in rows]
        corrected_vals = [r["corrected"] for r in rows]
        ambiguous_vals = [r["ambiguous"] for r in rows]
        no_match_vals = [r["no_match"] for r in rows]

        n_bcs = len(labels)

        fig_width = min(max(8, n_bcs * 0.12), 30)
        fig, ax = plt.subplots(figsize=(fig_width, 5))

        x = np.arange(n_bcs)

        ax.bar(x, exact_vals, width=1.0, label="exact", color="darkgreen")
        ax.bar(
            x,
            corrected_vals,
            width=1.0,
            bottom=exact_vals,
            label="corrected",
            color="lightgreen",
        )
        ax.bar(
            x,
            ambiguous_vals,
            width=1.0,
            bottom=[e + c for e, c in zip(exact_vals, corrected_vals)],
            label="ambiguous",
            color="orange",
        )
        ax.bar(
            x,
            no_match_vals,
            width=1.0,
            bottom=[e + c + a for e, c, a in zip(exact_vals, corrected_vals, ambiguous_vals)],
            label="no_match",
            color="grey",
        )

        ax.set_xlim(-0.5, n_bcs - 0.5)
        ax.set_ylabel("Read count")
        ax.set_xlabel("Barcode rank (sorted by exact count)")
        ax.set_title(f"QC metrics for {read_label}_{blockname}")
        ax.legend()
        ax.margins(x=0)

        # Always show rank values on x-axis
        if n_bcs <= 20:
            rank_step = 1
        elif n_bcs <= 100:
            rank_step = 5
        elif n_bcs <= 500:
            rank_step = 25
        else:
            rank_step = max(1, n_bcs // 20)

        rank_ticks = np.arange(0, n_bcs, rank_step)
        rank_ticklabels = [str(i + 1) for i in rank_ticks]
        ax.set_xticks(rank_ticks)
        ax.set_xticklabels(rank_ticklabels)

        # Optional barcode names as a second/top axis only when readable
        if n_bcs <= 20:
            ax_top = ax.twiny()
            ax_top.set_xlim(ax.get_xlim())
            ax_top.set_xticks(x)
            ax_top.set_xticklabels(labels, rotation=90, fontsize=10)
            ax_top.set_xlabel("Barcode")
        elif n_bcs <= 50:
            ax_top = ax.twiny()
            ax_top.set_xlim(ax.get_xlim())
            ax_top.set_xticks(x)
            ax_top.set_xticklabels(labels, rotation=90, fontsize=8)
            ax_top.set_xlabel("Barcode")

        plt.tight_layout()

        output_file = os.path.join(
            qc_dir,
            f"QC_metrics_barcodes_{read_label}_{blockname}.png"
        )
        plt.savefig(output_file, dpi=300)
        plt.close(fig)

        log.info(f"Wrote barcode-level stacked barplot to: {output_file}")

def make_stacked_barplot_reads(arguments):
    paths = get_qc_paths(arguments)
    summary_fpath = paths["reads_summary_tsv"]
    out_fpath = paths["reads_plot"]

    status_order = ["exact", "corrected", "ambiguous", "no_match"]
    counts = {status: 0 for status in status_order}

    with open(summary_fpath, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status = row["status"]
            count = int(row["count"])
            if status in counts:
                counts[status] = count

    total = sum(counts.values())

    if total > 0:
        percentages = {status: 100 * counts[status] / total for status in status_order}
    else:
        percentages = {status: 0 for status in status_order}

    bottom = 0
    fig, ax = plt.subplots(figsize=(4, 6))

    color_map = {
        "exact": "darkgreen",
        "corrected": "lightgreen",
        "ambiguous": "orange",
        "no_match": "grey",
    }

    for status in status_order:
        value = percentages[status]
        ax.bar(
            ["all_reads"],
            [value],
            bottom=bottom,
            color=color_map[status],
            label=status
        )
        bottom += value

    ax.set_ylabel("Reads (%)")
    ax.set_ylim(0, 100)
    ax.set_title("QC metrics across all reads")
    ax.legend()

    plt.tight_layout()
    plt.savefig(out_fpath, dpi=300)
    plt.close(fig)

    log.info(f"Wrote read-level stacked barplot to: {out_fpath}")