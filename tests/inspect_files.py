from pathlib import Path
import gzip
import csv
from scipy.io import mmread


def read_tsv_gz(filepath):
    rows = []
    with gzip.open(filepath, "rt", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def load_diagonal_matrix_dir(directory):
    """
    Load a directory containing exactly:
      - 1 .mtx.gz file
      - 2 .tsv.gz files

    Only keep TSV files whose length matches the matrix diagonal length.

    Returns
    -------
    dict with:
      - diagonal: list/array of diagonal values
      - matched_tsvs: dict of {filename: rows}
      - result: list of dicts, one per diagonal entry
    """
    directory = Path(directory)

    if not directory.is_dir():
        raise ValueError(f"Not a directory: {directory}")

    files = [f for f in directory.iterdir() if f.is_file()]
    if len(files) != 3:
        raise ValueError(f"Expected exactly 3 files in {directory}, found {len(files)}")

    mtx_files = [f for f in files if f.name.endswith(".mtx.gz")]
    tsv_files = sorted([f for f in files if f.name.endswith(".tsv.gz")])

    if len(mtx_files) != 1:
        raise ValueError(f"Expected exactly 1 .mtx.gz file, found {len(mtx_files)}")
    if len(tsv_files) != 2:
        raise ValueError(f"Expected exactly 2 .tsv.gz files, found {len(tsv_files)}")

    # load matrix
    with gzip.open(mtx_files[0], "rb") as f:
        matrix = mmread(f)

    diagonal = matrix.diagonal()
    n = len(diagonal)

    # load and keep only matching TSVs
    matched_tsvs = {}
    for tsv_file in tsv_files:
        rows = read_tsv_gz(tsv_file)
        if len(rows) == n:
            matched_tsvs[tsv_file.name] = rows
        else:
            print(f"Skipping {tsv_file.name}: length {len(rows)} does not match diagonal length {n}")

    if not matched_tsvs:
        raise ValueError(f"No TSV file matched diagonal length {n}")

    # build row-wise result
    result = []
    for i, value in enumerate(diagonal):
        row = {
            "index": i,
            "diagonal_value": value,
        }
        for fname, rows in matched_tsvs.items():
            row[fname] = rows[i]
        result.append(row)

    return {
        "diagonal": diagonal,
        "matched_tsvs": matched_tsvs,
        "result": result,
    }

import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
print(SCRIPT_DIR)
gDNA_input_ser = os.path.join(SCRIPT_DIR, "../tests/gDNA_single_end/raw_reads_bc_matrix")
gDNA_input = os.path.join(SCRIPT_DIR, "../tests/gDNA/raw_reads_bc_matrix")
print("##############File Paths##############")
print(gDNA_input_ser)
print(gDNA_input)
gDNA_ser = load_diagonal_matrix_dir(gDNA_input_ser)["result"]
gDNA = load_diagonal_matrix_dir(gDNA_input)["result"]

# print first 10 rows from both
for i, (row1, row2) in enumerate(zip(gDNA_ser[:10], gDNA[:10])):
    print(f"Row {i}")
    print("ser: ", row1)
    print("inp: ", row2)
    print()

print("----")

# then only print mismatches
for i, (row1, row2) in enumerate(zip(gDNA_ser, gDNA)):
    if row1 != row2:
        print(f"Mismatch at row {i}")
        print("ser: ", row1)
        print("inp: ", row2)
        print()


from collections import Counter
import pysam

def count_bam_querynames(bam_path, n=100000):
    c = Counter()
    with pysam.AlignmentFile(bam_path, "r") as bam:
        for i, read in enumerate(bam.fetch(until_eof=True)):
            c[read.query_name] += 1
            if i + 1 >= n:
                break
    return Counter(c.values())

print(count_bam_querynames("../examples/.bam"))