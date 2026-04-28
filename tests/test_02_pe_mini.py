# This Test is for the mini paired end data set, to see if it works with the STAR outsourced method
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

pe_mini_input_dir = os.path.join(SCRIPT_DIR, "../examples/pe_mini")
pe_mini_config = os.path.join(SCRIPT_DIR, "../examples/pe_mini_config.json")
star_index = os.path.join(SCRIPT_DIR, "../examples/pe_mini_genome_index")
star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"

USE_TEMP_OUTPUT = False


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


def format_long_tuples(rows, title):
    lines = [title, f"{'feature':<12} {'barcode':<10} {'count':>5}"]
    for feature, barcode, count in rows:
        lines.append(f"{feature:<12} {barcode:<10} {count:>5}")
    return "\n".join(lines)


def get_star_env():
    env = os.environ.copy()

    if "GITHUB_ACTIONS" in env:
        github_star_dir = "/home/runner/work/bcwithqc/bcwithqc/STAR-2.7.10b/source"
        star_executable = os.path.join(github_star_dir, "STAR")
        if not (os.path.isfile(star_executable) and os.access(star_executable, os.X_OK)):
            raise RuntimeError(f"STAR executable not found at GitHub Actions path '{star_executable}'")
        env["PATH"] = f"{github_star_dir}:{env.get('PATH', '')}"
    else:
        star_executable = os.path.join(star_dir_local, "STAR")
        if not (os.path.isfile(star_executable) and os.access(star_executable, os.X_OK)):
            raise RuntimeError(
                f"STAR executable not found or not executable at '{star_executable}'. "
                "Please set 'star_dir_local' to the correct local STAR installation path."
            )
        env["PATH"] = f"{star_dir_local}:{env.get('PATH', '')}"

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
def outsourced_star_output(request):
    threads = request.param
    env = get_star_env()

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix=f"pe_mini_star_outsourcing_threads{threads}_")
    else:
        output_dir = os.path.join(SCRIPT_DIR,"pe_mini" ,f"pe_mini_threadsN_{threads}")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        star_dir = os.path.join(output_dir, "STAR_files")
        os.makedirs(star_dir, exist_ok=True)

        # 1. bcwithqc preprocess
        preprocess_cmd = [
            "python", "-m", "bcwithqc", "preprocess",
            pe_mini_input_dir,
            f"--config={pe_mini_config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "-vvv",
        ]
        preprocess_result = run_command(preprocess_cmd, env)

        # Get all four input files (two pairs)
        preprocessed_files = [
            (os.path.join(output_dir, "sans_bc_pe_mini_1_r1.fq"),
             os.path.join(output_dir, "sans_bc_pe_mini_1_r2.fq")),
            (os.path.join(output_dir, "sans_bc_pe_mini_2_r1.fq"),
             os.path.join(output_dir, "sans_bc_pe_mini_2_r2.fq")),
        ]
        
        for pair_idx, (fastq_r1, fastq_r2) in enumerate(preprocessed_files, 1):
            if not os.path.isfile(fastq_r1):
                raise FileNotFoundError(f"Preprocessed fastq not found: {fastq_r1}")
            if not os.path.isfile(fastq_r2):
                raise FileNotFoundError(f"Preprocessed fastq not found: {fastq_r2}")

        # 2. external STAR alignment - run twice, once for each pair
        aligned_bams = []
        star_results = []
        
        for pair_idx, (fastq_r1, fastq_r2) in enumerate(preprocessed_files, 1):
            star_prefix = os.path.join(star_dir, f"pe_mini_pair{pair_idx}_")
            star_cmd = [
                'STAR',
                f'--runThreadN {threads}',
                f'--genomeDir {star_index}',
                f'--readFilesIn {fastq_r1} {fastq_r2}',
                f'--outFileNamePrefix {star_prefix}',
                '--outFilterMultimapNmax 1',
                '--outSAMtype BAM Unsorted',
                '--outSAMattributes NH HI AS nM GX GN',
            ]
            
            star_result = run_command(star_cmd, env)
            star_results.append(star_result)
            
            aligned_bam = f"{star_prefix}Aligned.out.bam"
            if not os.path.isfile(aligned_bam):
                raise FileNotFoundError(f"Expected STAR BAM not found: {aligned_bam}")
            aligned_bams.append(aligned_bam)

        # 3. bcwithqc count from existing STAR result
        count_cmd = [
            "python", "-m", "bcwithqc", "count",
            output_dir,
            f"--STAR-output-dir={star_dir}",
            f"--config={pe_mini_config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--keep-intermediary",
            "-vvv",
        ]
        count_result = run_command(count_cmd, env)

        yield {
            "output_dir": output_dir,
            "star_dir": star_dir,
            "preprocess_result": preprocess_result,
            "star_results": star_results,
            "count_result": count_result,
            "preprocessed_files": preprocessed_files,
            "aligned_bams": aligned_bams,
        }


# 1. test if "bcwithqc preprocess" runs
def test_preprocess_runs(outsourced_star_output):
    assert outsourced_star_output["preprocess_result"].returncode == 0


# 2. test if aligning the results of preprocess with STAR runs
def test_star_runs_on_preprocessed_output(outsourced_star_output):
    for star_result in outsourced_star_output["star_results"]:
        assert star_result.returncode == 0


# 3. test if "bcwithqc count" runs when using the STAR results
def test_count_runs_with_external_star_results(outsourced_star_output):
    assert outsourced_star_output["count_result"].returncode == 0


# 4. test if the correct files exist
def test_file_existence(outsourced_star_output):
    output_dir = outsourced_star_output["output_dir"]

    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_reads_bc_matrix", "features.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(output_dir, "raw_umis_bc_matrix", "features.tsv.gz"))


# 5. test contents, same logic as test_matrix_contents
def test_matrix_contents(outsourced_star_output):
    output_dir = outsourced_star_output["output_dir"]

    pe_reads_dir = os.path.join(output_dir, "raw_reads_bc_matrix")
    pe_umis_dir = os.path.join(output_dir, "raw_umis_bc_matrix")

    pe_reads = load_directory_files(pe_reads_dir)
    pe_umis = load_directory_files(pe_umis_dir)

    pe_reads_long = matrix_to_long_tuples(pe_reads)
    pe_umis_long = matrix_to_long_tuples(pe_umis)

    expected_reads = sorted([
        ("toy_gene", "AGCGTAGAA.AGCGTAGAA", 6),
        ("toy_gene", "CCTTAACAT.CCTTAACAT", 6),  # change to 8 if read_6 should pass
        ("toy_gene", "TATAGGTGT.TATAGGTGT", 2),
    ])

    expected_umis = sorted([
        ("toy_gene", "AGCGTAGAA.AGCGTAGAA", 2),
        ("toy_gene", "CCTTAACAT.CCTTAACAT", 3),  # change to 8 if read_6 should pass
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
