import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import nullcontext

import json

import pytest


# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_single_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_single_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA.json")

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
        star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"
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
def test_single_end_pipeline_runs(sample_type, input_dir, config):
    env = get_star_env()

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix=f"{sample_type}_single_end_STAR_out_")
    else:
        output_dir = os.path.join(SCRIPT_DIR, f"{sample_type}_single_end_STAR_out")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        star_dir = os.path.join(output_dir, "STAR_files")
        os.makedirs(star_dir, exist_ok=True)

        config_to_use = config
        temp_config = None

        try:
            # We need to edit the cDNA config and change keep_nonbarcode:false to true
            # for single-end reads, otherwise everything is discarded.
            if sample_type == "cDNA":
                with open(config, "r") as fh:
                    cfg = json.load(fh)

                cfg["barcode_struct_r1"]["keep_nonbarcode"] = True

                with tempfile.NamedTemporaryFile(
                    mode="w",
                    suffix=".json",
                    prefix="temp_cDNA_single_end_config_",
                    delete=False,
                ) as tf:
                    json.dump(cfg, tf, indent=4)
                    temp_config = tf.name

                config_to_use = temp_config

            # 1. bcwithqc preprocess
            preprocess_cmd = [
                "python", "-m", "bcwithqc", "preprocess",
                input_dir,
                f"--config={config_to_use}",
                f"--output-dir={output_dir}",
                "--threads=1",
                "--single-end-reads",
                "-vvv",
            ]
            preprocess_result = run_command(preprocess_cmd, env)

            if sample_type == "cDNA":
                preprocessed_fastq = os.path.join(output_dir, "sans_bc_cDNA_1_sequence.fq")
            if sample_type == "gDNA":
                preprocessed_fastq = os.path.join(output_dir, "sans_bc_gDNA_1_sequence.fq")
            
            if not os.path.isfile(preprocessed_fastq):
                raise FileNotFoundError(f"Expected preprocess FASTQ not found: {preprocessed_fastq}")

            # 2. external STAR alignment
            star_prefix = os.path.join(star_dir, f"{sample_type}_")
            star_cmd = [
                "STAR",
                f"--runThreadN 1",
                f"--genomeDir {star_index}",
                f"--readFilesIn {preprocessed_fastq}",
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
                f"--config={config_to_use}",
                f"--output-dir={output_dir}",
                "--threads=1",
                "--keep-intermediary",
                "--single-end-reads",
                "-vvv",
            ]
            count_result = run_command(count_cmd, env)

            assert preprocess_result.returncode == 0
            assert star_result.returncode == 0
            assert count_result.returncode == 0

            print(f"Output dir: {output_dir}")

        finally:
            if temp_config is not None:
                try:
                    os.remove(temp_config)
                except FileNotFoundError:
                    pass
