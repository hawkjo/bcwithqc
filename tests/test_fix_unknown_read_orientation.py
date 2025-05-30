import os
import pytest
import gzip
from bcwithqc.misc import fix_unknown_read_orientation


class DummyArgs:
    def __init__(self, output_dir, config):
        self.output_dir = output_dir
        self.config = config if config else {}

dummy_cDNA_config_r1_only = {
    "unknown_read_orientation": True,
    "barcode_struct_r1": {
        "keep_nonbarcode": False,
        "blocks": [
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "sequence": [
                    "CCTTAACAT",
                    "AGCGTAGAA",
                    "CATAGGTGT",
                    "CACGACTCA"
                ],
                "blocktype": "barcodeList",
                "blockname": "cb1"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "name": [
                    "cso1",
                    "cso2",
                    "cso3",
                    "cso4"
                ],
                "sequence": [
                    "GTCAGTACGTACGAGTC",
                    "TCAGTACGTACGAGTC",
                    "CAGTACGTACGAGTC",
                    "AGTACGTACGAGTC"
                ],
                "blocktype": "constantRegion",
                "blockname": "cso"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "sequence": [
                    "CCTTAACAT",
                    "AGCGTAGAA",
                    "CATAGGTGT",
                    "CACGACTCA",
                    "AGCATTCCA"
                ],
                "blocktype": "barcodeList",
                "blockname": "cb2"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "name": [
                    "cs2"
                ],
                "sequence": [
                    "GTACTCGCAGTAGTCGACACGTC"
                ],
                "blocktype": "constantRegion",
                "blockname": "commonseq2"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "sequence": [
                    "TCGCCTTA",
                    "CTAGTACG",
                    "TTCTGCCT",
                    "GCTCAGGA"
                ],
                "blocktype": "barcodeList",
                "blockname": "sb"
            },
            {
                "maxerrors": 0,
                "blockfunction": "UMI",
                "length": 8,
                "blocktype": "randomBarcode",
                "blockname": "umi"
            }
        ]
    }
}

dummy_gDNA_config_r1_only = {
    "unknown_read_orientation": True,
    "barcode_struct_r1": {
        "keep_nonbarcode": False,
        "blocks": [
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "sequence": [
                    "CCTTAACAT",
                    "AGCGTAGAA",
                    "CATAGGTGT"
                ],
                "blocktype": "barcodeList",
                "blockname": "bc1"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "name": [
                    "cso1",
                    "cso2",
                    "cso3",
                    "cso4"
                ],
                "sequence": [
                    "GTCAGTACGTACGAGTC",
                    "TCAGTACGTACGAGTC",
                    "CAGTACGTACGAGTC",
                    "AGTACGTACGAGTC"
                ],
                "blocktype": "constantRegion",
                "blockname": "commonseq1"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "sequence": [
                    "CCTTAACAT",
                    "AGCGTAGAA",
                    "CATAGGTGT",
                    "CACGACTCA",
                    "ACTGGACCA"
                ],
                "blocktype": "barcodeList",
                "blockname": "bc2"
            },
            {
                "maxerrors": 1,
                "blockfunction": "barcode",
                "name": [
                    "cso2"
                ],
                "sequence": [
                    "GTACTCGCAGTAGTC"
                ],
                "blocktype": "constantRegion",
                "blockname": "commonseq2"
            }
        ]
    }
}


# === Setup paths ===

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "shuffled_sequences/output")

paired_cDNA_fpaths = [
    (os.path.join(SCRIPT_DIR, "shuffled_sequences/cDNA_fastqs/cDNA_1_sequence_shuffled.txt.gz"),
     os.path.join(SCRIPT_DIR, "shuffled_sequences/cDNA_fastqs/cDNA_2_sequence_shuffled.txt.gz"))
]

paired_gDNA_fpaths = [
    (os.path.join(SCRIPT_DIR, "shuffled_sequences/gDNA_fastqs/gDNA_1_sequence_shuffled.txt.gz"),
     os.path.join(SCRIPT_DIR, "shuffled_sequences/gDNA_fastqs/gDNA_2_sequence_shuffled.txt.gz"))
]


# Expected original files
EXPECTED_FASTQS = {
    "gDNA": {
        "r1": os.path.join(SCRIPT_DIR, "shuffled_sequences/gDNA_fastqs/gDNA_1_sequence.txt.gz"),
        "r2": os.path.join(SCRIPT_DIR, "shuffled_sequences/gDNA_fastqs/gDNA_2_sequence.txt.gz"),
    },
    "cDNA": {
        "r1": os.path.join(SCRIPT_DIR, "shuffled_sequences/cDNA_fastqs/cDNA_1_sequence.txt.gz"),
        "r2": os.path.join(SCRIPT_DIR, "shuffled_sequences/cDNA_fastqs/cDNA_2_sequence.txt.gz"),
    }
}

# === Helper ===

def compare_gzipped_text_files(file1, file2, mismatch_threshold=0.002):
    mismatched_reads = 0
    total_reads = 0
    first_mismatch_index = None

    # Track whether the previous line was a read header
    last_is_read1 = False
    last_is_read2 = False

    with gzip.open(file1, 'rt', encoding='utf-8') as f1, gzip.open(file2, 'rt', encoding='utf-8') as f2:
        for i, (line1, line2) in enumerate(zip(f1, f2), 1):
            if line1.startswith('@') and line2.startswith('@'):
                total_reads += 1
                last_is_read1 = True
                last_is_read2 = True
                continue

            if last_is_read1 and last_is_read2:
                if line1 != line2:
                    mismatched_reads += 1
                    if first_mismatch_index is None:
                        first_mismatch_index = i

            last_is_read1 = line1.startswith('@')
            last_is_read2 = line2.startswith('@')

    # Compute mismatch rate
    mismatch_rate = mismatched_reads / total_reads if total_reads else 0

    # Always print a summary for test logs
    print(f"Total reads compared: {total_reads}")
    print(f"Mismatched reads: {mismatched_reads}")
    print(f"Mismatch rate: {mismatch_rate:.4%}")

    if mismatch_rate > mismatch_threshold:
        raise AssertionError(f"Mismatch rate {mismatch_rate:.4%} exceeds threshold of {mismatch_threshold*100:.2f}%\n"
                             f"Files: {file1} vs {file2}\n"
                             f"First mismatch at line {first_mismatch_index}.")



# === Test with parameterization ===
@pytest.mark.parametrize("sample_type", ["gDNA", "cDNA"])
def test_fix_unknown_read_orientation(sample_type):
    if sample_type == "gDNA":
        dummy_args = DummyArgs(OUTPUT_DIR, dummy_gDNA_config_r1_only)
        paired_fpaths = paired_gDNA_fpaths
    else:
        dummy_args = DummyArgs(OUTPUT_DIR, dummy_cDNA_config_r1_only)
        paired_fpaths = paired_cDNA_fpaths

    # Call the function with dummy arguments and input file paths
    new_fpaths = fix_unknown_read_orientation(dummy_args, paired_fpaths)

    if new_fpaths:
        first_path_tuple = new_fpaths[0]
        # Now you can use first_path_tuple in your assertions
    else:
        raise ValueError("new_fpaths is empty")

    # Expected output paths based on sample type
    expected_r1 = EXPECTED_FASTQS[sample_type]["r1"]
    expected_r2 = EXPECTED_FASTQS[sample_type]["r2"]

    # Check if one of the files in the output matches the expected files
    assert any(
        (compare_gzipped_text_files(first_path_tuple[0], expected_r1) == None and 
        compare_gzipped_text_files(first_path_tuple[1], expected_r2) == None) or
        (compare_gzipped_text_files(first_path_tuple[0], expected_r2) == None and 
        compare_gzipped_text_files(first_path_tuple[1], expected_r1) == None)
        for first_path_tuple in new_fpaths
    )

    # Check if unknown_read_orientation was set to False #
    assert not dummy_args.config["unknown_read_orientation"]
    assert os.path.exists(first_path_tuple[0])  # Ensure the first file exists
    assert os.path.exists(first_path_tuple[1])  # Ensure the second file exists
