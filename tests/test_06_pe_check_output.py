"""Compare the paired-end examples' preprocess and count outputs with stored copies.

count runs on the placeholder BAM written by ``preprocess --output-format-bam``,
in which every read pair is aligned to one feature, ``toy_gene``, rather than on
STAR output: from identical alignments, STAR builds on different platforms
assign genes (the GX/GN tags) to different reads, so stored outputs made with
STAR do not reproduce on every machine. The placeholder BAM contains every
read; test_07 runs count on a BAM that lacks some reads, as an aligner's does.
The reads preprocess writes for alignment (sans_bc_*.fq) are too large to
store, so they are compared by SHA-256 checksum.

The stored cDNA UMI matrix keeps a known bug: umi.get_umi_maps_from_bam_file
counts both reads of a pair, while the count matrices count one read per pair.
With every count doubled, UMIs that should merge, such as two single-read UMIs
one error apart, stay separate. The matrix has 16,649 UMIs where the single-end
example, with the same reads, has 16,589. test_umi_counts_match_single_end
asserts the intended result and is marked xfail until the bug is fixed.

The stored outputs in tests/gDNA_paired_end/ and tests/cDNA_paired_end/ are
written by tests/regenerate_expected.py.
"""
import os
import sys
import time
import shutil
import subprocess
import tempfile
from contextlib import nullcontext
from glob import glob
import gzip
import hashlib
import itertools

import pytest

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_paired_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_paired_end")

gDNA_expected = os.path.join(SCRIPT_DIR, "gDNA_paired_end")
cDNA_expected = os.path.join(SCRIPT_DIR, "cDNA_paired_end")

# Stored outputs of the single-end examples (test_05), which hold the same reads.
single_end_expected = {
    "gDNA": os.path.join(SCRIPT_DIR, "gDNA_single_end"),
    "cDNA": os.path.join(SCRIPT_DIR, "cDNA_single_end"),
}

gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA_paired_end_config.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA_paired_end_config.json")

verbosity = "-v"

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
                pass
    return tuple(snapshot)


def wait_for_tree_to_stabilize(path, stable_checks=3, interval=1.0, timeout=60):
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

    env = os.environ.copy()

    tmp_dir = os.path.join(SCRIPT_DIR, f"tmp_{sample_type}_paired_end_out_t{threads}")
    if os.path.exists(tmp_dir):
        rmtree_with_retries(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    output_dir = tmp_dir

    try:
        # 1. preprocess, also writing the placeholder BAM (see the module docstring)
        preprocess_cmd = [
            "python", "-m", "bcwithqc", "preprocess",
            input_dir,
            f"--config={config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--output-format-bam",
            verbosity,
        ]
        preprocess_result = run_command(preprocess_cmd, env)

        # 2. bcwithqc count from the placeholder BAM
        count_cmd = [
            "python", "-m", "bcwithqc", "count",
            output_dir,
            f"--STAR-output-dir={output_dir}",
            f"--config={config}",
            f"--output-dir={output_dir}",
            f"--threads={threads}",
            "--keep-intermediary",
            verbosity,
        ]
        count_result = run_command(count_cmd, env)

        wait_for_tree_to_stabilize(output_dir, stable_checks=3, interval=1.0, timeout=60)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Subprocess failed:\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Output dir: {output_dir}\n")
        sys.stderr.write(f"Command: {' '.join(e.cmd)}\n")
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise

    yield sample_type, threads, output_dir, expected_dir

    if os.path.exists(output_dir):
        rmtree_with_retries(output_dir)


def get_all_paths_filtered(base_dir):
    paths = set()
    for root, dirs, files in os.walk(base_dir):
        rel_root = os.path.relpath(root, base_dir)
        if rel_root == ".":
            rel_root = ""

        for d in dirs:
            paths.add(os.path.normpath(os.path.join(rel_root, d)))

        for f in files:
            rel_path = os.path.normpath(os.path.join(rel_root, f))
            if f.startswith("sans_bc_") and f.endswith(".fq"):
                continue
            paths.add(rel_path)

    return paths


def read_expected_paths(expected_dir):
    """Return the relative paths listed in ``expected_dir/expected_paths.txt``."""
    with open(os.path.join(expected_dir, "expected_paths.txt")) as fh:
        return {os.path.normpath(line.strip()) for line in fh if line.strip()}


def test_file_existence(sample_dirs):
    sample_type, threads, output_dir, expected_dir = sample_dirs

    # Only the content-compared files of the expected output are committed;
    # the full directory structure is recorded in expected_paths.txt.
    expected_paths = read_expected_paths(expected_dir)
    output_paths = get_all_paths_filtered(output_dir)

    missing_in_output = expected_paths - output_paths
    extra_in_output = output_paths - expected_paths

    assert expected_paths == output_paths, (
        f"Directory structures differ for sample {sample_type}:\n"
        f"Missing in output: {missing_in_output}\n"
        f"Extra in output: {extra_in_output}"
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


def sha256_of_file(fpath):
    """Return the hex SHA-256 digest of a file's bytes."""
    with open(fpath, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def test_preprocessed_fastqs(sample_dirs):
    """The reads preprocess writes for alignment match the stored checksums."""
    sample_type, threads, output_dir, expected_dir = sample_dirs
    checksums_fpath = os.path.join(expected_dir, "sans_bc_sha256.txt")
    with open(checksums_fpath) as fh:
        expected = {fname: digest for digest, fname in (line.split() for line in fh)}
    fastq_fpaths = glob(os.path.join(output_dir, "intermediary_files", "sans_bc_*.fq"))
    observed = {os.path.basename(fpath): sha256_of_file(fpath) for fpath in fastq_fpaths}
    assert observed == expected, f"Preprocessed reads differ from {checksums_fpath}"


def assert_same_matrix_as_single_end(sample_type, output_dir, matrix_dir):
    """Assert that a count matrix folder equals the stored single-end one.

    The single-end example FASTQs hold the same 20,000 reads as the paired-end
    ones, each made of the start of read 1 followed by read 2. preprocess finds
    the same barcodes and UMIs in both layouts, and the placeholder BAM has one
    feature, so both layouts should give the same counts.
    """
    for fname in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
        cmp_result = compare_gzipped_text_files(
            os.path.join(single_end_expected[sample_type], matrix_dir, fname),
            os.path.join(output_dir, matrix_dir, fname),
        )
        assert cmp_result["equal"], (
            f"{matrix_dir}/{fname} differs from the single-end example's "
            f"from line {cmp_result['first_diff_line']}"
        )


def test_read_counts_match_single_end(sample_dirs):
    sample_type, threads, output_dir, expected_dir = sample_dirs
    assert_same_matrix_as_single_end(sample_type, output_dir, "raw_reads_bc_matrix")


@pytest.mark.xfail(
    strict=True,
    reason="UMI correction counts both reads of a pair (see the module docstring)",
)
def test_umi_counts_match_single_end(sample_dirs):
    sample_type, threads, output_dir, expected_dir = sample_dirs
    if sample_type == "gDNA":
        pytest.skip("the gDNA config has no UMI block")
    assert_same_matrix_as_single_end(sample_type, output_dir, "raw_umis_bc_matrix")