import os
import sys
import re
import gzip
import glob
import shutil
import subprocess
import tempfile
from collections import Counter
from contextlib import nullcontext

import pytest

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

se_test_config = os.path.join(SCRIPT_DIR, "../examples/simulate_se_mini_config.json")
star_index = os.path.join(SCRIPT_DIR, "../examples/se_mini_genome_index")
star_dir_local = "/home/link/local/lib/STAR-2.7.11b/source"

USE_TEMP_OUTPUT = False
ERROR_SEGMENT_PATTERN = re.compile(r"S(\d+)D(\d+)I(\d+)")


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


def run_command(command, env=None, workdir=None, label="Subprocess"):
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
        sys.stderr.write(f"{label} failed:\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Command: {' '.join(command)}\n")
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise


def find_single_preprocessed_fastq(output_dir):
    candidates = glob.glob(os.path.join(output_dir, "sans_bc_*.fq"))
    candidates += glob.glob(os.path.join(output_dir, "sans_bc_*.fastq"))
    candidates += glob.glob(os.path.join(output_dir, "sans_bc_*.fq.gz"))
    candidates += glob.glob(os.path.join(output_dir, "sans_bc_*.fastq.gz"))

    assert len(candidates) == 1, (
        f"Expected exactly one preprocessed sans_bc FASTQ, "
        f"found {len(candidates)}: {candidates}"
    )
    return candidates[0]


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def write_error_combo_table(fastq_path, out_tsv):
    combo_counter = Counter()

    should_reject_2_plus_errors_in_any_block = 0
    first_block_1del_second_block_1ins = 0
    first_block_1ins_second_block_1del = 0
    should_accept_1_minus_errors_in_each_block = 0  # New counter for the 1-or-less errors condition
    
    # New counters for the extended indel checks
    first_block_1del_second_block_1error = 0
    first_block_1ins_second_block_1error = 0
    second_block_1del_first_block_1error = 0
    second_block_1ins_first_block_1error = 0

    with open_maybe_gzip(fastq_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue

            header = line.strip()
            if header.startswith("@"):
                header = header[1:]

            matches = ERROR_SEGMENT_PATTERN.findall(header)

            segments = tuple(f"S{s}D{d}I{ins}" for s, d, ins in matches)
            combo_counter[segments] += 1

            parsed = [(int(s), int(d), int(ins)) for s, d, ins in matches]

            # Update for reject condition (block has >= 2 errors)
            if any((s + d + ins) >= 2 for s, d, ins in parsed):
                should_reject_2_plus_errors_in_any_block += 1

            # Update for accept condition (each block has <= 1 error of any type)
            if all(s + d + ins <= 1 for s, d, ins in parsed):
                should_accept_1_minus_errors_in_each_block += 1

            # Check the conditions for 1 delete in the first block, 1 insert in the second block
            if len(parsed) >= 2:
                s1, d1, i1 = parsed[0]
                s2, d2, i2 = parsed[1]

                # First block: 1 delete, Second block: 1 insert
                if (d1 == 1 and s1 == 0 and i1 == 0) and (i2 == 1 and s2 == 0 and d2 == 0):
                    first_block_1del_second_block_1ins += 1

                # First block: 1 insert, Second block: 1 delete
                if (i1 == 1 and s1 == 0 and d1 == 0) and (d2 == 1 and s2 == 0 and i2 == 0):
                    first_block_1ins_second_block_1del += 1

                # First block: 1 delete, Second block: 1 or fewer errors of any type
                if (d1 == 1 and s1 == 0 and i1 == 0) and (s2 + d2 + i2 <= 1):
                    first_block_1del_second_block_1error += 1

                # First block: 1 insert, Second block: 1 or fewer errors of any type
                if (i1 == 1 and s1 == 0 and d1 == 0) and (s2 + d2 + i2 <= 1):
                    first_block_1ins_second_block_1error += 1

                # Second block: 1 delete, First block: 1 or fewer errors of any type
                if (d2 == 1 and s2 == 0 and i2 == 0) and (s1 + d1 + i1 <= 1):
                    second_block_1del_first_block_1error += 1

                # Second block: 1 insert, First block: 1 or fewer errors of any type
                if (i2 == 1 and s2 == 0 and d2 == 0) and (s1 + d1 + i1 <= 1):
                    second_block_1ins_first_block_1error += 1

    # Write the results to the output file
    with open(out_tsv, "w") as out:
        out.write("summary\tcount\n")
        out.write(f"Should Reject: 2+ errors in any block\t{should_reject_2_plus_errors_in_any_block}\n")
        out.write(f"Should Accept: 1- errors in each block\t{should_accept_1_minus_errors_in_each_block}\n")
        out.write("\n")
        out.write(f"Indel Issues\n")
        out.write(f"First block: 1D; Second block: 1I\t{first_block_1del_second_block_1ins}\n")
        out.write(f"first block: 1I; Second block: 1D\t{first_block_1ins_second_block_1del}\n")
        out.write(f"First block: 1D; Second block: 1 or fewer errors\t{first_block_1del_second_block_1error}\n")
        out.write(f"First block: 1I; Second block: 1 or fewer errors\t{first_block_1ins_second_block_1error}\n")
        out.write(f"Second block: 1D; First block: 1 or fewer errors\t{second_block_1del_first_block_1error}\n")
        out.write(f"Second block: 1I; First block: 1 or fewer errors\t{second_block_1ins_first_block_1error}\n")
        out.write("\n")
        out.write("count\tcombination\n")
        for combo, count in combo_counter.most_common():
            out.write(f"{count}\t{';'.join(combo)}\n")


@pytest.fixture(scope="module")
def simulate_se_mini_output():
    config = se_test_config
    config_dir = os.path.dirname(os.path.abspath(config))
    config_stem = os.path.splitext(os.path.basename(config))[0]

    env = get_star_env()

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix="simulate_se_mini_pipeline_")
    else:
        root_dir = os.path.join(SCRIPT_DIR, "simulate_se_mini")
        if os.path.exists(root_dir):
            shutil.rmtree(root_dir)
        os.makedirs(root_dir, exist_ok=True)
        context = nullcontext(root_dir)

    with context as root_dir:
        sim_dir = os.path.join(root_dir, "simulated_fastq")
        count_dir = os.path.join(root_dir, "count_output")
        star_dir = os.path.join(root_dir, "STAR_files")
        os.makedirs(sim_dir, exist_ok=True)
        os.makedirs(count_dir, exist_ok=True)
        os.makedirs(star_dir, exist_ok=True)

        simulated_fastq = os.path.join(sim_dir, f"{config_stem}.txt.gz")

        simulate_command = [
            "python", "-m", "bcwithqc", "simulate_reads",
            f"--config={config}",
            f"--output-dir={sim_dir}",
            "--nreads=10000",
            "--error-probability=0.1",
            "--substitution-probability=0.7",
            "-vvv",
        ]

        run_command(simulate_command, label="Simulation subprocess")

        assert os.path.isfile(simulated_fastq), f"Expected simulated FASTQ not found: {simulated_fastq}"

        # 1. bcwithqc preprocess
        # This avoids relying on the internal STAR execution path in `bcwithqc count`.
        preprocess_command = [
            "python", "-m", "bcwithqc", "preprocess",
            sim_dir,
            f"--config={config}",
            f"--output-dir={count_dir}",
            "--threads=1",
            "-vvv",
        ]
        run_command(preprocess_command, env=env, label="Preprocess subprocess")

        sans_bc_fastq = find_single_preprocessed_fastq(count_dir)

        # 2. external STAR alignment
        star_prefix = os.path.join(star_dir, "simulate_se_mini_")
        star_command = [
            "STAR",
            "--runThreadN", "1",
            "--genomeDir", star_index,
            "--readFilesIn", sans_bc_fastq,
            "--outFileNamePrefix", star_prefix,
            "--outFilterMultimapNmax", "1",
            "--outSAMtype", "BAM", "Unsorted",
            "--outSAMattributes", "NH", "HI", "AS", "nM", "GX", "GN",
        ]
        run_command(star_command, env=env, label="STAR subprocess")

        aligned_bam = f"{star_prefix}Aligned.out.bam"
        assert os.path.isfile(aligned_bam), f"Expected STAR BAM not found: {aligned_bam}"

        # 3. bcwithqc count from existing STAR result
        count_command = [
            "python", "-m", "bcwithqc", "count",
            count_dir,
            f"--STAR-output-dir={star_dir}",
            f"--config={config}",
            f"--output-dir={count_dir}",
            "--threads=1",
            "--keep-intermediary",
            "-vvv",
        ]
        run_command(count_command, env=env, label="Count subprocess")

        before_tsv = os.path.join(count_dir, "simulate_se_mini_before.tsv")
        after_tsv = os.path.join(count_dir, "simulate_se_mini_after.tsv")
        sans_bc_fastq = find_single_preprocessed_fastq(os.path.join(count_dir, "intermediary_files"))

        write_error_combo_table(simulated_fastq, before_tsv)
        write_error_combo_table(sans_bc_fastq, after_tsv)

        yield {
            "root_dir": root_dir,
            "sim_dir": sim_dir,
            "count_dir": count_dir,
            "star_dir": star_dir,
            "simulated_fastq": simulated_fastq,
            "sans_bc_fastq": sans_bc_fastq,
            "before_tsv": before_tsv,
            "after_tsv": after_tsv,
        }


def test_simulation_and_count_outputs_exist(simulate_se_mini_output):
    count_dir = simulate_se_mini_output["count_dir"]

    assert os.path.isfile(os.path.join(count_dir, "raw_reads_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(count_dir, "raw_reads_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(count_dir, "raw_reads_bc_matrix", "features.tsv.gz"))

    assert os.path.isfile(os.path.join(count_dir, "raw_umis_bc_matrix", "matrix.mtx.gz"))
    assert os.path.isfile(os.path.join(count_dir, "raw_umis_bc_matrix", "barcodes.tsv.gz"))
    assert os.path.isfile(os.path.join(count_dir, "raw_umis_bc_matrix", "features.tsv.gz"))

    assert os.path.isfile(simulate_se_mini_output["simulated_fastq"])
    assert os.path.isfile(simulate_se_mini_output["sans_bc_fastq"])
    assert os.path.isfile(simulate_se_mini_output["before_tsv"])
    assert os.path.isfile(simulate_se_mini_output["after_tsv"])


def test_error_combo_tables_are_nonempty(simulate_se_mini_output):
    for path in [simulate_se_mini_output["before_tsv"], simulate_se_mini_output["after_tsv"]]:
        with open(path, "r") as fh:
            lines = [line.rstrip("\n") for line in fh]

        assert len(lines) >= 9, f"{path} should contain summary lines, a blank line, and at least one data row"

        # Check summary section
        assert lines[0] == "summary\tcount"
        assert lines[1].startswith("Should Reject: 2+ errors in any block\t")
        assert lines[2].startswith("Should Accept: 1- errors in each block\t")
        assert lines[3] == ""
        
        # Check indel issue section
        assert lines[4] == "Indel Issues"
        assert lines[5].startswith("First block: 1D; Second block: 1I\t")
        assert lines[6].startswith("first block: 1I; Second block: 1D\t")
        assert lines[7].startswith("First block: 1D; Second block: 1 or fewer errors\t")
        assert lines[8].startswith("First block: 1I; Second block: 1 or fewer errors\t")
        assert lines[9].startswith("Second block: 1D; First block: 1 or fewer errors\t")
        assert lines[10].startswith("Second block: 1I; First block: 1 or fewer errors\t")

        # Ensure there's a blank line after the indel issues section
        assert lines[11] == ""

        # Check combination section
        assert lines[12] == "count\tcombination"
        
        assert len(lines) >= 13, f"{path} should contain at least one combination row after the header"