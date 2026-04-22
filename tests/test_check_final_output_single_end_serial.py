import os
import sys
import time
import shutil
import subprocess
import pytest
import gzip
import json
import tempfile
from filecmp import cmp

# Get the directory this test script lives in
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Paths to input directories (relative to script)
gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_single_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_single_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA.json")


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


@pytest.fixture(scope="module", params=["gDNA", "cDNA"])
def sample_dirs(request):
    sample_type = request.param
    tmp_dir = os.path.join(SCRIPT_DIR, f"tmp_{sample_type}_single_end")
    expected_dir = os.path.join(SCRIPT_DIR, f"{sample_type}_single_end")
    input_dir = gDNA_input if sample_type == "gDNA" else cDNA_input
    config = gDNA_config if sample_type == "gDNA" else cDNA_config

    if os.path.exists(tmp_dir):
        rmtree_with_retries(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)

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

    config_to_use = config
    temp_config = None

    try:
        # For cDNA single-end, keep_nonbarcode must be True or everything is discarded
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

        command = [
            "python", "-m", "bcwithqc", "count",
            input_dir,
            f"--STAR-ref-dir={star_index}",
            f"--config={config_to_use}",
            f"--output-dir={tmp_dir}",
            "--threads=1",
            "--keep-intermediary",
            "--single-end-reads",
            "-vvv"
        ]

        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            text=True
        )
        print(result.stdout)
        print(result.stderr)

        wait_for_tree_to_stabilize(tmp_dir, stable_checks=3, interval=1.0, timeout=60)

    except subprocess.CalledProcessError as e:
        sys.stderr.write("Subprocess failed:\n")
        sys.stderr.write(f"Return code: {e.returncode}\n")
        sys.stderr.write(f"Output dir: {tmp_dir}\n")
        sys.stderr.write(f"Command: {' '.join(command)}\n")
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise

    yield sample_type, tmp_dir, expected_dir

    if os.path.exists(tmp_dir):
        rmtree_with_retries(tmp_dir)

    if temp_config is not None:
        try:
            os.remove(temp_config)
        except FileNotFoundError:
            pass


def test_file_existence(sample_dirs):
    sample_type, tmp_dir, expected_dir = sample_dirs

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
    tmp_paths = get_all_paths(tmp_dir)

    missing_in_tmp = expected_paths - tmp_paths
    extra_in_tmp = tmp_paths - expected_paths

    assert expected_paths == tmp_paths, (
        f"Directory structures differ for sample {sample_type}:\n"
        f"Missing in tmp: {missing_in_tmp}\n"
        f"Extra in tmp: {extra_in_tmp}"
    )


def test_file_contents(sample_dirs):
    sample_type, tmp_dir, expected_dir = sample_dirs

    def is_in_intermediary(path):
        return 'intermediary_files' in path.split(os.sep)

    def compare_gzipped_text_files(file1, file2):
        with gzip.open(file1, 'rt', encoding='utf-8') as f1, gzip.open(file2, 'rt', encoding='utf-8') as f2:
            return f1.read() == f2.read()

    def assert_tsv_mtx_files_equal(expected_dir, tmp_dir):
        for root, dirs, files in os.walk(expected_dir):
            rel_root = os.path.relpath(root, expected_dir)
            if is_in_intermediary(rel_root):
                continue

            for f in files:
                if not (f.endswith('.tsv.gz') or f.endswith('.mtx.gz')):
                    continue

                expected_path = os.path.join(root, f)
                tmp_path = os.path.join(tmp_dir, rel_root, f)

                assert os.path.exists(tmp_path), f"File missing in tmp_dir: {tmp_path}"

                if not compare_gzipped_text_files(expected_path, tmp_path):
                    raise AssertionError(f"File contents differ: {expected_path} vs {tmp_path}")

    assert_tsv_mtx_files_equal(expected_dir, tmp_dir)