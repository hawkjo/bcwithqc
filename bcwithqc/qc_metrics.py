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
        "blocks_plot": os.path.join(qc_dir, "QC_stacked_barplot_blocks.png"),
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

    Each bar corresponds to one barcode from QC_metrics_bcs.tsv.
    Bars are sorted by descending exact count.

    Colors:
    - dark green: exact
    - light green: corrected
    - orange: ambiguous
    - grey: no_match
    """
    paths = get_qc_paths(arguments)
    qc_dir = paths["qc_dir"]
    input_file = paths["bcs_tsv"]

    grouped_rows = defaultdict(list)

    with open(input_file, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        for row in reader:
            barcode = row["barcode"]

            # Skip special rows from plotting.
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
        ax.set_title(f"QC metrics for {read_label}_{blockname}")
        ax.legend()

        if n_bcs <= 20:
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=10)
        elif n_bcs <= 50:
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=8)
        elif n_bcs <= 100:
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=6)
        else:
            ax.set_xticks([])
            ax.set_xlabel("Barcodes")

        ax.margins(x=0)

        plt.tight_layout()

        output_file = os.path.join(qc_dir, f"QC_metrics_barcodes_{read_label}_{blockname}.png")
        plt.savefig(output_file, dpi=300)
        plt.close(fig)

        log.info(f"Wrote barcode-level stacked barplot to: {output_file}")
