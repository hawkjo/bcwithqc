import os
import sys
import shutil
import subprocess
import tempfile
import gzip
from contextlib import nullcontext

import pytest
from scipy.io import mmread


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

se_test_input = os.path.join(SCRIPT_DIR, "../examples/se_test")
se_test_config = os.path.join(SCRIPT_DIR, "../examples/se_test.json")
star_index = os.path.join(SCRIPT_DIR, "../examples/se_test_genome_index")

USE_TEMP_OUTPUT = True
USE_STAR_REF = True


def read_tsv_gz_first_col(path):
    values = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            values.append(line.split("\t")[0])
    return values


def read_barcodes(dir_path):
    return read_tsv_gz_first_col(os.path.join(dir_path, "barcodes.tsv.gz"))


def read_features(dir_path):
    return read_tsv_gz_first_col(os.path.join(dir_path, "features.tsv.gz"))


def read_matrix(dir_path):
    return mmread(os.path.join(dir_path, "matrix.mtx.gz")).tocoo()


def load_directory_files(dir_path):
    return {
        "barcodes": read_barcodes(dir_path),
        "features": read_features(dir_path),
        "matrix": read_matrix(dir_path),
    }


def matrix_to_long_tuples(x):
    mat = x["matrix"]
    rows = []
    for i, j, v in zip(mat.row, mat.col, mat.data):
        rows.append(
            (
                x["features"][i],
                x["barcodes"][j],
                int(v),
            )
        )
    return sorted(rows)


@pytest.fixture(scope="module")
def single_end_output_dir():
    env = os.environ.copy()

    if "GITHUB_ACTIONS" in env:
        github_star_dir = "/home/runner/work/bcwithqc/bcwithqc/STAR-2.7.10b/source"
        star_executable = os.path.join(github_star_dir, "STAR")
        if not (os.path.isfile(star_executable) and os.access(star_executable, os.X_OK)):
            raise RuntimeError(f"STAR executable not found at GitHub Actions path '{star_executable}'")
        env["PATH"] = f"{github_star_dir}:{env.get('PATH', '')}"
    else:
        star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"
        star_executable = os.path.join(star_dir_local, "STAR")
        if not (os.path.isfile(star_executable) and os.access(star_executable, os.X_OK)):
            raise RuntimeError(
                f"STAR executable not found or not executable at '{star_executable}'. "
                "Please set 'star_dir_local' to the correct local STAR installation path."
            )
        env["PATH"] = f"{star_dir_local}:{env.get('PATH', '')}"

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix="se_test_single_end_advanced")
    else:
        output_dir = os.path.join(SCRIPT_DIR, "se_test")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        command = [
            "python", "-m", "bcwithqc", "count",
            se_test_input,
            f"--config={se_test_config}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--single-end-reads",
            "--keep-intermediary",
            "--block-type-for-STAR-alignment=constant",
            "-vvv",
        ]
        if USE_STAR_REF:
            command.insert(4, f"--STAR-ref-dir={star_index}")

        try:
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                text=True,
            )
            print(result.stdout)
            print(result.stderr)
            print(f"Output dir: {output_dir}")
            yield output_dir
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Subprocess failed:\n")
            sys.stderr.write(f"Return code: {e.returncode}\n")
            sys.stderr.write(f"Output dir: {output_dir}\n")
            sys.stderr.write(f"Command: {' '.join(command)}\n")
            sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
            sys.stderr.write(f"STDERR:\n{e.stderr}\n")
            sys.stderr.flush()
            raise


def test_file_existence(single_end_output_dir):
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_reads_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_reads_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_reads_bc_matrix", "features.tsv.gz"))
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_umis_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_umis_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(single_end_output_dir, "raw_umis_bc_matrix", "features.tsv.gz"))


def format_long_tuples(rows, title):
    lines = [title, f"{'feature':<12} {'barcode':<10} {'count':>5}"]
    for feature, barcode, count in rows:
        lines.append(f"{feature:<12} {barcode:<10} {count:>5}")
    return "\n".join(lines)


def test_matrix_contents(single_end_output_dir):
    se_reads_dir = os.path.join(single_end_output_dir, "raw_reads_bc_matrix")
    se_umis_dir = os.path.join(single_end_output_dir, "raw_umis_bc_matrix")

    se_reads = load_directory_files(se_reads_dir)
    se_umis = load_directory_files(se_umis_dir)

    se_reads_long = matrix_to_long_tuples(se_reads)
    se_umis_long = matrix_to_long_tuples(se_umis)

    expected_reads = sorted([
        ("toy_gene", "AGCGTAGAA", 3),
        ("toy_gene", "CCTTAACAT", 4), #This needs to be changed to 5, if read_6 is allowed to pass
        ("toy_gene", "TATAGGTGT", 1),
    ])

    expected_umis = sorted([
        ("toy_gene", "AGCGTAGAA", 2),
        ("toy_gene", "CCTTAACAT", 4), #This needs to be changed to 5, if read_6 is allowed to pass
        ("toy_gene", "TATAGGTGT", 1),
    ])

    if se_reads_long != expected_reads:
        pytest.fail(
            "\n"
            + format_long_tuples(expected_reads, "Expected reads")
            + "\n\n"
            + format_long_tuples(se_reads_long, "Produced reads")
        )

    if se_umis_long != expected_umis:
        pytest.fail(
            "\n"
            + format_long_tuples(expected_umis, "Expected UMIs")
            + "\n\n"
            + format_long_tuples(se_umis_long, "Produced UMIs")
        )