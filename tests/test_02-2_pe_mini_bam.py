# This Test mirrors pe_mini but uses --output-format-bam to skip STAR
import os
import sys
import shutil
import subprocess
import tempfile
import gzip
from glob import glob
from contextlib import nullcontext

import pytest
from scipy.io import mmread


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

pe_mini_input_dir = os.path.join(SCRIPT_DIR, "../examples/pe_mini")
pe_mini_config = os.path.join(SCRIPT_DIR, "../examples/pe_mini_config.json")
star_index = os.path.join(SCRIPT_DIR, "../examples/pe_mini_genome_index")
star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"

USE_TEMP_OUTPUT = False
verbosity = "-vvv"


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


def get_env():
    env = os.environ.copy()
    return env


def run_command(command, env, workdir=None):
    try:
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            cwd=workdir,
            text=True,
        )
        print("COMMAND:", " ".join(command))
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
        return result
    except subprocess.CalledProcessError as e:
        sys.stderr.write("Subprocess failed:\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Command: {' '.join(command)}\n")
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise


@pytest.fixture(
    scope="module",
    params=[1, 2],
    ids=["serial", "parallel"]
)
def preprocess_bam_output(request):
    threads = request.param
    env = get_env()

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix=f"pe_mini_bam_star_outsourcing_threads{threads}_")
    else:
        output_dir = os.path.join(SCRIPT_DIR, "pe_mini_bam", f"pe_mini_threadsN_{threads}")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        star_dir = os.path.join(output_dir, "STAR_files")
        os.makedirs(star_dir, exist_ok=True)

        # 1. bcwithqc preprocess with output-format-bam
        preprocess_cmd = [
            "python", "-m", "bcwithqc", "preprocess",
            pe_mini_input_dir,
            f"--config={pe_mini_config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--output-format-bam",
            verbosity,
        ]
        preprocess_result = run_command(preprocess_cmd, env)

        # Expect STAR-like BAMs created by preprocess in STAR_files
        bam_glob = os.path.join(output_dir, "*Aligned.out.bam")
        matches = [p for p in glob(bam_glob)]
        if not matches:
            raise FileNotFoundError(f"Preprocess did not produce BAM in output_dir: {bam_glob}")

        # 2. run count using output_dir as STAR output dir
        count_cmd = [
            "python", "-m", "bcwithqc", "count",
            output_dir,
            f"--STAR-output-dir={output_dir}",
            f"--config={pe_mini_config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--keep-intermediary",
            verbosity,
        ]
        count_result = run_command(count_cmd, env)

        yield {
            "output_dir": output_dir,
            "preprocess_result": preprocess_result,
            "count_result": count_result,
        }


def test_preprocess_and_count_run(preprocess_bam_output):
    assert preprocess_bam_output["preprocess_result"].returncode == 0
    assert preprocess_bam_output["count_result"].returncode == 0


# 4. test if the correct files exist
def test_file_existence(preprocess_bam_output):
    output_dir = preprocess_bam_output["output_dir"]

    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "features.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "features.tsv.gz"))


# 5. test contents, same logic as test_matrix_contents
def test_matrix_contents(preprocess_bam_output):
    output_dir = preprocess_bam_output["output_dir"]

    pe_reads_dir = os.path.join(output_dir, "raw_reads_bc_matrix")
    pe_umis_dir = os.path.join(output_dir, "raw_umis_bc_matrix")

    pe_reads = load_directory_files(pe_reads_dir)
    pe_umis = load_directory_files(pe_umis_dir)

    pe_reads_long = matrix_to_long_tuples(pe_reads)
    pe_umis_long = matrix_to_long_tuples(pe_umis)

    expected_reads = sorted([
        ("toy_gene", "AGCGTAGAA.AGCGTAGAA", 6),
        ("toy_gene", "CCTTAACAT.CCTTAACAT", 10),  # change to 10 if read_6 should pass
        ("toy_gene", "TATAGGTGT.TATAGGTGT", 2),
    ])

    expected_umis = sorted([
        ("toy_gene", "AGCGTAGAA.AGCGTAGAA", 2),
        ("toy_gene", "CCTTAACAT.CCTTAACAT", 5),  # change to 5 if read_6 should pass
        ("toy_gene", "TATAGGTGT.TATAGGTGT", 1),
    ])

    if pe_reads_long != expected_reads:
        pytest.fail(
            "\n"
            + format_long_tuples(expected_reads, "Expected reads")
            + "\n\n"
            + format_long_tuples(pe_reads_long, "Produced reads")
        )

    if pe_umis_long != expected_umis:
        pytest.fail(
            "\n"
            + format_long_tuples(expected_umis, "Expected UMIs")
            + "\n\n"
            + format_long_tuples(pe_umis_long, "Produced UMIs")
        )
