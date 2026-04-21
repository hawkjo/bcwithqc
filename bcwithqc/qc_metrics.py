import logging
import os
import pysam
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pickle
import subprocess
import shutil
import scipy
import csv
import warnings
from . import misc
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from contextlib import nullcontext
from glob import glob
from multiprocessing import Pool
from threading import BoundedSemaphore
from scipy.sparse import lil_matrix
from .bc_aligner import CustomBCAligner
from .umi import get_umi_maps_from_bam_file


log = logging.getLogger(__name__)
pysam.set_verbosity(0)

def generate_qc_metrics(arguments, read_qcs):
    """
    Generate a TSV file with barcode QC metrics in a QC_metrics subdirectory.
    """
    log.info(f"Collected {len(read_qcs):,d} read QC entries")

    qc_dir = os.path.join(arguments.output_dir, "QC_metrics")
    os.makedirs(qc_dir, exist_ok=True)

    output_file = os.path.join(qc_dir, "QC_metrics.tsv")
    status_columns = ["exact", "corrected", "ambiguous", "no_match"]

    r1_meta = get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = get_config_metadata(arguments, "barcode_struct_r2")

    counts = defaultdict(Counter)
    conflict_counts = defaultdict(Counter)
    block_summary_counts = defaultdict(Counter)
    row_keys = []

    # Pre-initialize whitelist rows plus special rows
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

    with open(output_file, "w") as out_fh:
        header = ["read", "blockname", "barcode"] + status_columns + ["total", "conflict_bcs"]
        out_fh.write("\t".join(header) + "\n")

        for key in row_keys:
            read_label, blockname, barcode_label = key
            status_counter = counts[key]
            total = sum(status_counter.values())

            conflict_counter = conflict_counts.get(key, Counter())
            if conflict_counter:
                conflict_summary = ";".join(
                    f"{entry}:{count}"
                    for entry, count in sorted(conflict_counter.items())
                )
            else:
                conflict_summary = ""

            row = [
                read_label,
                blockname,
                barcode_label,
                *[str(status_counter.get(status, 0)) for status in status_columns],
                str(total),
                conflict_summary,
            ]
            out_fh.write("\t".join(row) + "\n")

        out_fh.write("\n")
        out_fh.write("SUMMARY\n")
        out_fh.write("read\tblockname\texact\tcorrected\tambiguous\tno_match\ttotal\n")

        for read_label, blockname in sorted(block_summary_counts.keys()):
            summary = block_summary_counts[(read_label, blockname)]
            exact = summary.get("exact", 0)
            corrected = summary.get("corrected", 0)
            ambiguous = summary.get("ambiguous", 0)
            no_match = summary.get("no_match", 0)
            total = exact + corrected + ambiguous + no_match

            out_fh.write(
                "\t".join([
                    read_label,
                    blockname,
                    str(exact),
                    str(corrected),
                    str(ambiguous),
                    str(no_match),
                    str(total),
                ]) + "\n"
            )

    log.info(f"Wrote QC metrics to: {output_file}")


def make_stacked_barplot_blocks(arguments):
    """
    Read the SUMMARY section of QC_metrics.tsv and create a stacked barplot.

    Dark green  = exact
    Light green = corrected
    Orange      = ambiguous
    Grey        = no_match

    One bar per barcode block.
    Y-axis is shown as percentage.
    Output is written to the QC_metrics subdirectory.
    """
    qc_dir = os.path.join(arguments.output_dir, "QC_metrics")
    os.makedirs(qc_dir, exist_ok=True)

    input_file = os.path.join(qc_dir, "QC_metrics.tsv")
    output_file = os.path.join(qc_dir, "QC_stacked_barplot_blocks.png")

    labels = []
    exact_pct = []
    corrected_pct = []
    ambiguous_pct = []
    no_match_pct = []

    in_summary = False
    with open(input_file, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            if line == "SUMMARY":
                in_summary = True
                next(fh)  # skip summary header
                continue

            if not in_summary or not line:
                continue

            read_label, blockname, exact, corrected, ambiguous, no_match, total = line.split("\t")

            exact = int(exact)
            corrected = int(corrected)
            ambiguous = int(ambiguous)
            no_match = int(no_match)
            total = int(total)

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
    x = range(len(labels))

    ax.bar(x, exact_pct, label="exact", color="darkgreen")
    ax.bar(
        x,
        corrected_pct,
        bottom=exact_pct,
        label="corrected",
        color="lightgreen",
    )
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

    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Reads (%)")
    ax.set_ylim(0, 100)
    ax.set_title("QC metrics per barcode block")
    ax.legend()

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    log.info(f"Wrote stacked barplot to: {output_file}")

def make_stacked_barplots_bcs(arguments):
    """
    Create one stacked barplot per barcode block, with one bar per whitelist barcode.

    Each bar corresponds to one barcode from the 'barcode' column in QC_metrics.tsv.
    Bars are sorted by descending exact count.

    Colors:
    - dark green: exact
    - light green: corrected
    - orange: ambiguous
    - grey: no_match

    Plots are written into the QC_metrics subdirectory.
    """

    import csv

    qc_dir = os.path.join(arguments.output_dir, "QC_metrics")
    os.makedirs(qc_dir, exist_ok=True)

    input_file = os.path.join(qc_dir, "QC_metrics.tsv")

    # Collect rows grouped by (read, blockname)
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

        # Scale figure width mildly with barcode number, but cap it
        fig_width = min(max(8, n_bcs * 0.12), 30)
        fig, ax = plt.subplots(figsize=(fig_width, 5))

        x = np.arange(n_bcs)

        # width=1.0 removes whitespace between bars
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

        # Adaptive x tick label handling
        if n_bcs <= 20:
            fontsize = 10
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)
        elif n_bcs <= 50:
            fontsize = 8
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)
        elif n_bcs <= 100:
            fontsize = 6
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)
        else:
            # Too many labels -> hide them
            ax.set_xticks([])
            ax.set_xlabel("Barcodes")

        # Remove margins so bars touch plot edges nicely
        ax.margins(x=0)

        plt.tight_layout()

        output_file = os.path.join(
            qc_dir,
            f"QC_metrics_barcodes_{read_label}_{blockname}.png"
        )
        plt.savefig(output_file, dpi=300)
        plt.close(fig)

        log.info(f"Wrote barcode-level stacked barplot to: {output_file}")