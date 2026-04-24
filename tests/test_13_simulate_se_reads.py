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

USE_TEMP_OUTPUT = False
USE_STAR_REF = True

ERROR_SEGMENT_PATTERN = re.compile(r"S(\d+)D(\d+)I(\d+)")


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def write_error_combo_table(fastq_path, out_tsv):
    combo_counter = Counter()

    any_block_ge_2_errors = 0
    first_block_1del_second_block_1ins = 0
    first_block_1ins_second_block_1del = 0

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

            if any((s + d + ins) >= 2 for s, d, ins in parsed):
                any_block_ge_2_errors += 1

            if len(parsed) >= 2:
                s1, d1, i1 = parsed[0]
                s2, d2, i2 = parsed[1]

                if (d1 == 1 and s1 == 0 and i1 == 0) and (i2 == 1 and s2 == 0 and d2 == 0):
                    first_block_1del_second_block_1ins += 1

                if (i1 == 1 and s1 == 0 and d1 == 0) and (d2 == 1 and s2 == 0 and i2 == 0):
                    first_block_1ins_second_block_1del += 1

    with open(out_tsv, "w") as out:
        out.write("summary\tcount\n")
        out.write(f"any_block_ge_2_total_errors\t{any_block_ge_2_errors}\n")
        out.write(f"first_block_exactly_1D_second_block_exactly_1I\t{first_block_1del_second_block_1ins}\n")
        out.write(f"first_block_exactly_1I_second_block_exactly_1D\t{first_block_1ins_second_block_1del}\n")
        out.write("\n")
        out.write("count\tcombination\n")
        for combo, count in combo_counter.most_common():
            out.write(f"{count}\t{';'.join(combo)}\n")


@pytest.fixture(scope="module")
def simulate_se_mini_output():
    config = se_test_config
    config_dir = os.path.dirname(os.path.abspath(config))
    config_stem = os.path.splitext(os.path.basename(config))[0]

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
        os.makedirs(sim_dir, exist_ok=True)
        os.makedirs(count_dir, exist_ok=True)

        simulated_fastq = os.path.join(sim_dir, f"{config_stem}.txt.gz")

        simulate_command = [
            "python", "-m", "bcwithqc", "simulate_reads",
            f"--config={config}",
            f"--output-dir={sim_dir}",
            "--nreads=10000",
            "--error-probability=0.1",
            "--substitution-probability=0.7",
            "--single-end-reads",
            "-vvv",
        ]

        try:
            result = subprocess.run(
                simulate_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            print(result.stdout)
            print(result.stderr)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Simulation subprocess failed:\n")
            sys.stderr.write(f"Return code: {e.returncode}\n")
            sys.stderr.write(f"Simulation dir: {sim_dir}\n")
            sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
            sys.stderr.write(f"STDERR:\n{e.stderr}\n")
            sys.stderr.flush()
            raise

        assert os.path.isfile(simulated_fastq), f"Expected simulated FASTQ not found: {simulated_fastq}"

        count_command = [
            "python", "-m", "bcwithqc", "count",
            sim_dir,
            f"--config={config}",
            f"--output-dir={count_dir}",
            "--threads=1",
            "--single-end-reads",
            "--keep-intermediary",
            "--block-type-for-STAR-alignment=constant",
            "-vvv",
        ]
        if USE_STAR_REF:
            count_command.insert(4, f"--STAR-ref-dir={star_index}")

        try:
            result = subprocess.run(
                count_command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=env,
                text=True,
            )
            print(result.stdout)
            print(result.stderr)
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Count subprocess failed:\n")
            sys.stderr.write(f"Return code: {e.returncode}\n")
            sys.stderr.write(f"Count output dir: {count_dir}\n")
            sys.stderr.write(f"Command: {' '.join(count_command)}\n")
            sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
            sys.stderr.write(f"STDERR:\n{e.stderr}\n")
            sys.stderr.flush()
            raise

        before_tsv = os.path.join(count_dir, "simulate_se_mini_before.tsv")
        after_tsv = os.path.join(count_dir, "simulate_se_mini_after.tsv")

        write_error_combo_table(simulated_fastq, before_tsv)

        sans_bc_candidates = glob.glob(
            os.path.join(count_dir, "intermediary_files", "sans_bc_*.fq")
        )
        assert len(sans_bc_candidates) == 1, (
            f"Expected exactly one sans_bc FASTQ, found {len(sans_bc_candidates)}: {sans_bc_candidates}"
        )
        sans_bc_fastq = sans_bc_candidates[0]

        write_error_combo_table(sans_bc_fastq, after_tsv)

        yield {
            "root_dir": root_dir,
            "sim_dir": sim_dir,
            "count_dir": count_dir,
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

        assert len(lines) >= 6, f"{path} should contain summary lines, a blank line, and at least one data row"

        assert lines[0] == "summary\tcount"
        assert lines[1].startswith("any_block_ge_2_total_errors\t")
        assert lines[2].startswith("first_block_exactly_1D_second_block_exactly_1I\t")
        assert lines[3].startswith("first_block_exactly_1I_second_block_exactly_1D\t")
        assert lines[4] == ""
        assert lines[5] == "count\tcombination"

        assert len(lines) >= 7, f"{path} should contain at least one combination row after the header"