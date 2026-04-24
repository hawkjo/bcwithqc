import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import nullcontext

import pytest

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_paired_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_paired_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA_paired_end_config.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA_paired_end_config.json")

# Absolut Path to where STAR is installed -> THIS NEEDS TO BE ADJUSTED FOR LOCAL MACHINES!
star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"

# Set to True to use a temp directory that is deleted automatically.
# Set to False to write into SCRIPT_DIR/<sample_type>_single_end_STAR_out and keep the output.
USE_TEMP_OUTPUT = True


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


@pytest.mark.parametrize(
    "sample_type,input_dir,config",
    [
        ("gDNA", gDNA_input, gDNA_config),
        ("cDNA", cDNA_input, cDNA_config),
    ],
)
def test_paired_end_pipeline_runs(sample_type, input_dir, config):
    env = get_star_env()

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix=f"{sample_type}_paired_end_basic_STAR_out_")
    else:
        output_dir = os.path.join(SCRIPT_DIR, f"{sample_type}_paired_end")
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
            input_dir,
            f"--config={config}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "-vvv",
        ]
        preprocess_result = run_command(preprocess_cmd, env)

        if sample_type == "cDNA":
            preprocessed_fastq = os.path.join(output_dir, "sans_bc_cDNA_2_sequence.fq")
            if not os.path.isfile(preprocessed_fastq):
                raise FileNotFoundError(f"Expected preprocess FASTQ not found: {preprocessed_fastq}")
            readFilesIn_line = f"--readFilesIn {preprocessed_fastq}"
        if sample_type == "gDNA":
            preprocessed_fastq_1 = os.path.join(output_dir, "sans_bc_gDNA_1_sequence.fq")
            preprocessed_fastq_2 = os.path.join(output_dir, "sans_bc_gDNA_2_sequence.fq")
            if not os.path.isfile(preprocessed_fastq_1):
                raise FileNotFoundError(f"Expected preprocess FASTQ not found: {preprocessed_fastq_1}")
            if not os.path.isfile(preprocessed_fastq_2):
                raise FileNotFoundError(f"Expected preprocess FASTQ not found: {preprocessed_fastq_2}")
            readFilesIn_line = f"--readFilesIn {preprocessed_fastq_1} {preprocessed_fastq_2}"


        # 2. external STAR alignment
        star_prefix = os.path.join(star_dir, f"{sample_type}_")
        star_cmd = [
            "STAR",
            f"--runThreadN 1",
            f"--genomeDir {star_index}",
            readFilesIn_line,
            f"--outFileNamePrefix {star_prefix}",
            "--outFilterMultimapNmax 1",
            "--outSAMtype BAM Unsorted",
            "--outSAMattributes NH HI AS nM GX GN",
        ]
        star_result = run_command(star_cmd, env)

        aligned_bam = f"{star_prefix}Aligned.out.bam"
        if not os.path.isfile(aligned_bam):
            raise FileNotFoundError(f"Expected STAR BAM not found: {aligned_bam}")

        # 3. bcwithqc count from existing STAR result
        count_cmd = [
            "python", "-m", "bcwithqc", "count",
            output_dir,
            f"--STAR-output-dir={star_dir}",
            f"--config={config}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--keep-intermediary",
            "-vvv",
        ]
        count_result = run_command(count_cmd, env)

        assert preprocess_result.returncode == 0
        assert star_result.returncode == 0
        assert count_result.returncode == 0
