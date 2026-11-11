"""Run count on a BAM that lacks some of the preprocessed reads.

An aligner leaves the reads it cannot align out of its BAM (STAR's default), so
count has to match each remaining BAM record to its own line in preprocess's
tags file. The placeholder BAM from ``preprocess --output-format-bam`` contains
every read, so test_05 and test_06 never exercise this. Here preprocess runs on
the se_mini example, four reads are removed from the placeholder BAM, and count
runs on the rest, with one thread (tags matched in file order) and with two
(tags matched by read name).

The expected values follow from the read names in examples/se_mini/se_mini.txt,
which state each read's barcode, UMI and errors. read_3 has two substitutions
in its barcode and fails preprocessing. Removing read_1 (the first record),
read_4 and read_5 (two in a row) and read_10 (the last record, and the only read
of TATAGGTGT) leaves these reads:

    read_2  CCTTAACAT  UMI CCCCCC
    read_6  CCTTAACAT  UMI GCGCGC
    read_7  AGCGTAGAA  UMI AAAAAA
    read_8  AGCGTAGAA  UMI AAAAAA
    read_9  AGCGTAGAA  UMI CCCCCC

No two of these UMIs are within one error of each other, so UMI correction
merges only the two AAAAAA reads.
"""

import gzip
import os
import subprocess
import sys
from glob import glob

import pysam
import pytest
from scipy.io import mmread

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPT_DIR)
SE_MINI_INPUT_DIR = os.path.join(REPO_DIR, "examples", "se_mini")
SE_MINI_CONFIG = os.path.join(REPO_DIR, "examples", "se_mini_config.json")

REMOVED_READS = {"read_1", "read_4", "read_5", "read_10"}
EXPECTED_BARCODE_BY_READ = {
    "read_2": "CCTTAACAT",
    "read_6": "CCTTAACAT",
    "read_7": "AGCGTAGAA",
    "read_8": "AGCGTAGAA",
    "read_9": "AGCGTAGAA",
}
EXPECTED_READS_PER_BARCODE = {"AGCGTAGAA": 3, "CCTTAACAT": 2}
EXPECTED_UMIS_PER_BARCODE = {"AGCGTAGAA": 2, "CCTTAACAT": 2}


def read_id(read: pysam.AlignedSegment) -> str:
    """Return the se_mini read id, the part of the read name before the first ':'."""
    return read.query_name.split(":")[0]


def run_bcwithqc(arguments: list[str]) -> None:
    """Run ``python -m bcwithqc`` from this checkout and fail the test if it fails."""
    result = subprocess.run(
        [sys.executable, "-m", "bcwithqc", *arguments],
        cwd=REPO_DIR,
        capture_output=True,
        text=True,
        check=False,
    )
    assert (
        result.returncode == 0
    ), f"bcwithqc {arguments[0]} failed:\n{result.stdout}\n{result.stderr}"


def remove_reads_from_bam(bam_fpath: str, removed_read_ids: set[str]) -> None:
    """Rewrite a BAM file without the records of the given read ids."""
    kept_fpath = bam_fpath + ".kept.bam"
    with pysam.AlignmentFile(bam_fpath, "rb") as bam_in:
        with pysam.AlignmentFile(kept_fpath, "wb", template=bam_in) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                if read_id(read) not in removed_read_ids:
                    bam_out.write(read)
    os.replace(kept_fpath, bam_fpath)


def counts_per_barcode(matrix_dir: str) -> dict[str, int]:
    """Return each barcode's count, summed over features, from a count matrix folder."""
    with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
        barcodes = fh.read().splitlines()
    # Rows are features and columns are barcodes.
    column_sums = mmread(os.path.join(matrix_dir, "matrix.mtx.gz")).sum(axis=0)
    return dict(zip(barcodes, column_sums.astype(int).tolist()[0]))


@pytest.fixture(scope="module", params=[1, 2], ids=["serial", "parallel"])
def count_output_dir(request, tmp_path_factory):
    threads = request.param
    output_dir = str(tmp_path_factory.mktemp(f"se_mini_unaligned_t{threads}"))
    run_bcwithqc(
        [
            "preprocess",
            SE_MINI_INPUT_DIR,
            f"--config={SE_MINI_CONFIG}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--output-format-bam",
        ]
    )
    (bam_fpath,) = glob(os.path.join(output_dir, "*Aligned.out.bam"))
    remove_reads_from_bam(bam_fpath, REMOVED_READS)
    run_bcwithqc(
        [
            "count",
            output_dir,
            f"--STAR-output-dir={output_dir}",
            f"--config={SE_MINI_CONFIG}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
        ]
    )
    return output_dir


def test_each_read_keeps_its_own_barcode(count_output_dir):
    bam_fpath = os.path.join(count_output_dir, "with_bc_umi.sorted.bam")
    with pysam.AlignmentFile(bam_fpath, "rb") as bam:
        barcode_by_read = [
            (read_id(read), read.get_tag("CB")) for read in bam.fetch(until_eof=True)
        ]
    assert sorted(barcode_by_read) == sorted(EXPECTED_BARCODE_BY_READ.items())


def test_reads_per_barcode(count_output_dir):
    matrix_dir = os.path.join(count_output_dir, "raw_reads_bc_matrix")
    assert counts_per_barcode(matrix_dir) == EXPECTED_READS_PER_BARCODE


def test_umis_per_barcode(count_output_dir):
    matrix_dir = os.path.join(count_output_dir, "raw_umis_bc_matrix")
    assert counts_per_barcode(matrix_dir) == EXPECTED_UMIS_PER_BARCODE
