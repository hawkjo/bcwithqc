import os
import sys
import shutil
import subprocess
import tempfile
import warnings
import gzip
from contextlib import nullcontext

import pytest

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA.json")

# Set to True to use a temp directory that is deleted automatically.
# Set to False to write into SCRIPT_DIR/gDNA and SCRIPT_DIR/cDNA and keep the output.
#
# WARNING:
# When set to False, expected_dir == output_dir, so the file existence and file
# content comparison tests become tautological and will always pass as long as
# the pipeline run itself succeeds.
USE_TEMP_OUTPUT = True


@pytest.fixture(scope="module", params=["gDNA", "cDNA"])
def sample_dirs(request):
    sample_type = request.param

    input_dir = gDNA_input if sample_type == "gDNA" else cDNA_input
    config = gDNA_config if sample_type == "gDNA" else cDNA_config
    expected_dir = os.path.join(SCRIPT_DIR, sample_type)

    env = os.environ.copy()
    if "GITHUB_ACTIONS" in env:
        # On GitHub Actions, add the installed STAR path explicitly
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
        context = tempfile.TemporaryDirectory(prefix=f"{sample_type}_")
        cleanup_output = True
    else:
        output_dir = os.path.join(SCRIPT_DIR, sample_type)

        warnings.warn(
            f"USE_TEMP_OUTPUT=False: writing output directly into '{output_dir}'. "
            "This overwrites the expected output directory, so the file existence "
            "and content comparison tests will always pass if the pipeline run succeeds.",
            UserWarning,
        )

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        context = nullcontext(output_dir)
        cleanup_output = False

    with context as output_dir:
        command = [
            "python", "-m", "bcwithqc", "count",
            input_dir,
            f"--STAR-ref-dir={star_index}",
            f"--config={config}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--keep-intermediary",
            "-vvv",
        ]

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
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Subprocess failed:\n")
            sys.stderr.write(f"Return code: {e.returncode}\n")
            sys.stderr.write(f"Output dir: {output_dir}\n")
            sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
            sys.stderr.write(f"STDERR:\n{e.stderr}\n")
            sys.stderr.flush()
            raise

        yield sample_type, output_dir, expected_dir

    if USE_TEMP_OUTPUT and cleanup_output and os.path.exists(output_dir):
        shutil.rmtree(output_dir)


def test_file_existence(sample_dirs):
    sample_type, output_dir, expected_dir = sample_dirs

    def get_all_paths(base_dir):
        paths = set()
        for root, dirs, files in os.walk(base_dir):
            rel_root = os.path.relpath(root, base_dir)
            if rel_root == ".":
                rel_root = ""
            for d in dirs:
                paths.add(os.path.normpath(os.path.join(rel_root, d)))
            for f in files:
                paths.add(os.path.normpath(os.path.join(rel_root, f)))
        return paths

    expected_paths = get_all_paths(expected_dir)
    output_paths = get_all_paths(output_dir)

    missing_in_output = expected_paths - output_paths
    extra_in_output = output_paths - expected_paths

    assert expected_paths == output_paths, (
        f"Directory structures differ for sample {sample_type}:\n"
        f"Missing in output: {missing_in_output}\n"
        f"Extra in output: {extra_in_output}"
    )


def test_file_contents(sample_dirs):
    sample_type, output_dir, expected_dir = sample_dirs

    def is_in_intermediary(path):
        return "intermediary_files" in path.split(os.sep)

    def compare_gzipped_text_files(file1, file2):
        with gzip.open(file1, "rt", encoding="utf-8") as f1, gzip.open(file2, "rt", encoding="utf-8") as f2:
            return f1.read() == f2.read()

    for root, dirs, files in os.walk(expected_dir):
        rel_root = os.path.relpath(root, expected_dir)
        if is_in_intermediary(rel_root):
            continue

        for f in files:
            if not (f.endswith(".tsv.gz") or f.endswith(".mtx.gz")):
                continue

            expected_path = os.path.join(root, f)
            output_path = os.path.join(output_dir, rel_root, f)

            assert os.path.exists(output_path), f"File missing in output_dir: {output_path}"

            if not compare_gzipped_text_files(expected_path, output_path):
                raise AssertionError(
                    f"File contents differ: {expected_path} vs {output_path}"
                )