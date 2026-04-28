import os
import sys
import time
import shutil
import subprocess
import tempfile
from contextlib import nullcontext
import gzip
import itertools

import json

import pytest

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_single_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_single_end")

gDNA_expected = os.path.join(SCRIPT_DIR, "gDNA_single_end")
cDNA_expected = os.path.join(SCRIPT_DIR, "cDNA_single_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA_single_end_config.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA_single_end_config.json")


def snapshot_tree(path):
    snapshot = []
    for root, dirs, files in os.walk(path):
        for name in sorted(dirs + files):
            p = os.path.join(root, name)
            try:
                st = os.stat(p)
                snapshot.append((
                    os.path.relpath(p, path),
                    st.st_size,
                    int(st.st_mtime_ns),
                ))
            except FileNotFoundError:
                # Something changed during scan; ignore and try again next loop
                pass
    return tuple(snapshot)


def wait_for_tree_to_stabilize(path, stable_checks=3, interval=1.0, timeout=60):
    """
    Wait until the directory tree stops changing across several consecutive checks.
    """
    deadline = time.time() + timeout
    last = None
    stable_count = 0

    while time.time() < deadline:
        current = snapshot_tree(path)
        if current == last:
            stable_count += 1
            if stable_count >= stable_checks:
                return
        else:
            stable_count = 0
            last = current
        time.sleep(interval)

    raise TimeoutError(f"Output directory did not stabilize within {timeout}s: {path}")


def rmtree_with_retries(path, retries=10, delay=1.0):
    last_err = None
    for _ in range(retries):
        try:
            shutil.rmtree(path)
            return
        except OSError as e:
            last_err = e
            time.sleep(delay)
    raise last_err


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


@pytest.fixture(
    scope="module",
    params=[
        ("gDNA", 1),
        ("gDNA", 2),
        ("cDNA", 1),
        ("cDNA", 2),
    ],
    ids=[
        "gDNA_serial",
        "gDNA_parallel",
        "cDNA_serial",
        "cDNA_parallel",
    ],
)
def sample_dirs(request):
    sample_type, threads = request.param
    input_dir = gDNA_input if sample_type == "gDNA" else cDNA_input
    config = gDNA_config if sample_type == "gDNA" else cDNA_config
    expected_dir = gDNA_expected if sample_type == "gDNA" else cDNA_expected

    tmp_dir = os.path.join(SCRIPT_DIR, f"tmp_{sample_type}_single_end_STAR_out_t{threads}")
    output_dir = tmp_dir
    if os.path.exists(tmp_dir):
        rmtree_with_retries(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

    star_dir = os.path.join(tmp_dir, "STAR_files")
    os.makedirs(star_dir, exist_ok=True)

    env = get_star_env()
    try:
        # 1. preprocess
        preprocess_cmd = [
            "python", "-m", "bcwithqc", "preprocess",
            input_dir,
            f"--config={config}",
            f"--output-dir={tmp_dir}",
            f"--threads={threads}",
            "-vvv",
        ]
        preprocess_result = run_command(preprocess_cmd, env)

        if sample_type == "cDNA":
            preprocessed_fastq = os.path.join(output_dir, "sans_bc_cDNA_1_sequence.fq")
        if sample_type == "gDNA":
            preprocessed_fastq = os.path.join(output_dir, "sans_bc_gDNA_1_sequence.fq")
        
        if not os.path.isfile(preprocessed_fastq):
            raise FileNotFoundError(f"Expected preprocess FASTQ not found: {preprocessed_fastq}")

        # 2. run STAR
        star_prefix = os.path.join(star_dir, f"{sample_type}_")
        star_cmd = [
            "STAR",
            f"--runThreadN {threads}",
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

        # 3. count from STAR output
        count_cmd = [
            "python", "-m", "bcwithqc", "count",
            tmp_dir,
            f"--STAR-output-dir={star_dir}",
            f"--config={config}",
            f"--output-dir={tmp_dir}",
            f"--threads={threads}",
            "--keep-intermediary",
            "-vvv",
        ]
        count_result = run_command(count_cmd, env)

        assert preprocess_result.returncode == 0
        assert star_result.returncode == 0
        assert count_result.returncode == 0

        wait_for_tree_to_stabilize(tmp_dir, stable_checks=3, interval=1.0, timeout=60)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Subprocess failed:\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Output dir: {tmp_dir}\n")
        sys.stderr.write(f"Command: {' '.join(count_cmd if 'count_cmd' in locals() else preprocess_cmd)}\n")
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise

    yield sample_type, threads, output_dir, expected_dir

    if os.path.exists(tmp_dir):
        rmtree_with_retries(tmp_dir)



def should_ignore_relpath(rel_path):
    parts = rel_path.split(os.sep)
    base = os.path.basename(rel_path)

    if "STAR_files" in parts:
        return True

    if base.startswith("sans_bc_") and base.endswith(".fq"):
        return True

    return False


def test_file_existence(sample_dirs):
    sample_type, threads, output_dir, expected_dir = sample_dirs

    def get_all_paths(base_dir):
        paths = set()
        for root, dirs, files in os.walk(base_dir):
            rel_root = os.path.relpath(root, base_dir)
            if rel_root == ".":
                rel_root = ""
            for d in dirs:
                rel_path = os.path.normpath(os.path.join(rel_root, d))
                if not should_ignore_relpath(rel_path):
                    paths.add(rel_path)
            for f in files:
                rel_path = os.path.normpath(os.path.join(rel_root, f))
                if not should_ignore_relpath(rel_path):
                    paths.add(rel_path)
        return paths

    expected_paths = get_all_paths(expected_dir)
    tmp_paths = get_all_paths(output_dir)

    missing_in_tmp = expected_paths - tmp_paths
    extra_in_tmp = tmp_paths - expected_paths

    assert expected_paths == tmp_paths, (
        f"Directory structures differ for sample {sample_type}:\n"
        f"Missing in tmp: {missing_in_tmp}\n"
        f"Extra in tmp: {extra_in_tmp}"
    )


def open_maybe_gzip(path, mode="rt", encoding="utf-8"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding=encoding)
    return open(path, mode, encoding=encoding)

def compare_gzipped_text_files(file1, file2, preview_lines=10):
    """
    Compare two gzipped text files line by line.

    Returns a dict with:
      - equal: bool
      - first_diff_line: 1-based line number of first difference, or None
      - expected_line: line content from file1 at first difference
      - output_line: line content from file2 at first difference
      - expected_preview: up to preview_lines lines from file1 starting at first difference
      - output_preview: up to preview_lines lines from file2 starting at first difference
      - expected_total_lines: total lines in file1
      - output_total_lines: total lines in file2
    """
    with open_maybe_gzip(file1, "rt", encoding="utf-8") as f1, open_maybe_gzip(file2, "rt", encoding="utf-8") as f2:
        line_num = 0

        while True:
            l1 = f1.readline()
            l2 = f2.readline()

            if not l1 and not l2:
                return {
                    "equal": True,
                    "first_diff_line": None,
                    "expected_line": None,
                    "output_line": None,
                    "expected_preview": [],
                    "output_preview": [],
                    "expected_total_lines": line_num,
                    "output_total_lines": line_num,
                }

            line_num += 1

            if l1 != l2:
                expected_preview = [l1.rstrip("\n")]
                output_preview = [l2.rstrip("\n")]

                for _ in range(preview_lines - 1):
                    next1 = f1.readline()
                    next2 = f2.readline()

                    if next1:
                        expected_preview.append(next1.rstrip("\n"))
                    if next2:
                        output_preview.append(next2.rstrip("\n"))

                    if not next1 and not next2:
                        break

                expected_total = line_num + sum(1 for _ in f1)
                output_total = line_num + sum(1 for _ in f2)

                return {
                    "equal": False,
                    "first_diff_line": line_num,
                    "expected_line": l1.rstrip("\n"),
                    "output_line": l2.rstrip("\n"),
                    "expected_preview": expected_preview,
                    "output_preview": output_preview,
                    "expected_total_lines": expected_total,
                    "output_total_lines": output_total,
                }

def test_file_contents(sample_dirs):
    sample_type, threads, output_dir, expected_dir = sample_dirs

    def is_in_intermediary(path):
        return "intermediary_files" in path.split(os.sep)

    def format_preview(lines, title):
        out = [title]
        for i, line in enumerate(lines, start=1):
            out.append(f"{i:>2}: {line}")
        return "\n".join(out)

    files_that_passed = []

    for root, dirs, files in os.walk(expected_dir):
        dirs.sort()
        files.sort()

        rel_root = os.path.relpath(root, expected_dir)
        if is_in_intermediary(rel_root):
            continue

        for f in files:
            if not (f.endswith(".tsv.gz") or f.endswith(".mtx.gz") or f.endswith(".tsv")):
                continue

            expected_path = os.path.join(root, f)
            output_path = os.path.join(output_dir, rel_root, f)
            rel_file = os.path.normpath(os.path.join(rel_root, f))

            assert os.path.exists(output_path), f"File missing in output_dir: {output_path}"

            cmp_result = compare_gzipped_text_files(expected_path, output_path, preview_lines=10)

            if not cmp_result["equal"]:
                raise AssertionError(
                    f"File contents differ for sample {sample_type}, threads={threads}:\n"
                    f"Expected: {expected_path}\n"
                    f"Output:   {output_path}\n"
                    f"First difference at line: {cmp_result['first_diff_line']}\n"
                    f"Expected total lines: {cmp_result['expected_total_lines']}\n"
                    f"Output total lines:   {cmp_result['output_total_lines']}\n\n"
                    f"Expected line:\n{cmp_result['expected_line']}\n\n"
                    f"Output line:\n{cmp_result['output_line']}\n\n"
                    f"{format_preview(cmp_result['expected_preview'], 'Expected preview')}\n\n"
                    f"{format_preview(cmp_result['output_preview'], 'Output preview')}\n\n"
                    f"{format_preview(files_that_passed, 'Files that passed content check before first failure:')}"
                )

            files_that_passed.append(rel_file)
