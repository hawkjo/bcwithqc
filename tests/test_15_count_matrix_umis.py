"""Tests for UMI-aware and UMI-free count matrix construction (PR #2).

The BAM is built by hand, so the expected matrices follow directly from the
reads listed in ``_bam_reads``.
"""

import gzip
import os
from types import SimpleNamespace

import numpy as np
import pysam
import pytest
import scipy.io

from bcwithqc.count import count_matrix
from bcwithqc.misc import config_has_umis

# (cell barcode, UMI, GX, GN, flag)
# flag 0 = unpaired, 65 = paired read 1, 129 = paired read 2 with a mapped mate.
_bam_reads = [
    ("AAA", "U1", "g1", "G1", 0),
    ("AAA", "U1", "g1", "G1", 0),
    ("AAA", "U2", "g1", "G1", 0),
    ("CCC", "U1", "g1;g2", "G1;G2", 0),  # counts toward both genes
    ("CCC", "U4", None, None, 0),  # no gene assignment: barcode seen, nothing counted
    ("CCC", "U3", "g2", "G2", 65),  # paired read 1: counted
    ("CCC", "U3", "g2", "G2", 129),  # paired read 2, mate mapped: not counted
]

# Rows: features sorted [(g1, G1), (g2, G2)]; columns: barcodes sorted [AAA, CCC].
expected_read_counts = np.array([[3, 1], [0, 2]])
expected_umi_counts = np.array([[2, 1], [0, 2]])


def _write_sorted_bam(bam_fpath: str, include_umi: bool) -> None:
    """Write the reads in ``_bam_reads`` to a coordinate-sorted, indexed BAM."""
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 10000}],
    }
    unsorted_fpath = bam_fpath + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_fpath, "wb", header=header) as bam_out:
        for read_index, (cell_barcode, umi, gx, gn, flag) in enumerate(_bam_reads):
            segment = pysam.AlignedSegment()
            segment.query_name = f"read{read_index}"
            segment.query_sequence = "ACGT" * 5
            segment.flag = flag
            segment.reference_id = 0
            segment.reference_start = 100 * read_index
            segment.mapping_quality = 255
            segment.cigartuples = [(0, 20)]
            segment.query_qualities = pysam.qualitystring_to_array("I" * 20)
            if flag & 1:
                segment.next_reference_id = 0
                segment.next_reference_start = 100 * read_index
            segment.set_tag("CB", cell_barcode)
            if include_umi:
                segment.set_tag("UB", umi)
            if gx is not None:
                segment.set_tag("GX", gx)
                segment.set_tag("GN", gn)
            bam_out.write(segment)
    pysam.sort("-o", bam_fpath, unsorted_fpath)
    pysam.index(bam_fpath)
    os.remove(unsorted_fpath)


def _read_matrix_dir(matrix_dir: str):
    """Return (dense matrix, barcodes, feature ids) from a written matrix directory."""
    with gzip.open(os.path.join(matrix_dir, "matrix.mtx.gz")) as fh:
        matrix = scipy.io.mmread(fh).toarray()
    with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
        barcodes = fh.read().split("\n")
    with gzip.open(os.path.join(matrix_dir, "features.tsv.gz"), "rt") as fh:
        feature_ids = [line.split("\t")[0] for line in fh.read().split("\n")]
    return matrix, barcodes, feature_ids


@pytest.mark.parametrize("has_umis", [True, False])
def test_count_matrix_read_and_umi_counts(tmp_path, has_umis):
    bam_fpath = str(tmp_path / "with_bc_umi.sorted.bam")
    _write_sorted_bam(bam_fpath, include_umi=has_umis)
    output_dir = tmp_path / "out"
    arguments = SimpleNamespace(output_dir=str(output_dir), threads=1)

    count_matrix(arguments, bam_fpath, has_umis)

    read_matrix, barcodes, feature_ids = _read_matrix_dir(
        str(output_dir / "raw_reads_bc_matrix")
    )
    assert barcodes == ["AAA", "CCC"]
    assert feature_ids == ["g1", "g2"]
    np.testing.assert_array_equal(read_matrix, expected_read_counts)

    umi_dir = output_dir / "raw_umis_bc_matrix"
    if has_umis:
        umi_matrix, umi_barcodes, umi_feature_ids = _read_matrix_dir(str(umi_dir))
        assert umi_barcodes == barcodes
        assert umi_feature_ids == feature_ids
        np.testing.assert_array_equal(umi_matrix, expected_umi_counts)
    else:
        assert not umi_dir.exists()


def _block(blockfunction: str) -> dict:
    return {"blockfunction": blockfunction, "blocktype": "randomBarcode", "length": 6}


@pytest.mark.parametrize(
    "config, expected",
    [
        ({"barcode_struct_r1": {"blocks": [_block("barcode"), _block("UMI")]}}, True),
        (
            {
                "barcode_struct_r1": {"blocks": [_block("barcode")]},
                "barcode_struct_r2": {"blocks": [_block("UMI")]},
            },
            True,
        ),
        (
            {
                "barcode_struct_r1": {"blocks": [_block("barcode"), _block("discard")]},
                "barcode_struct_r2": {"blocks": [_block("barcode")]},
            },
            False,
        ),
        ({"barcode_struct_r1": {"blocks": [_block("barcode")]}}, False),
        (
            {"barcode_struct_r1": {"blocks": []}, "barcode_struct_r2": {"blocks": []}},
            False,
        ),
    ],
    ids=["umi_on_r1", "umi_on_r2_only", "no_umi_both_reads", "no_r2", "empty_blocks"],
)
def test_config_has_umis(config, expected):
    assert config_has_umis(config) is expected
