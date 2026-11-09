"""Regenerate the stored outputs that test_05 and test_06 compare against.

test_05_se_check_output.py and test_06_pe_check_output.py run
``bcwithqc preprocess --output-format-bam`` and ``bcwithqc count`` on the four
example datasets and compare the results with the files stored in
``tests/{gDNA,cDNA}_{single,paired}_end/``. This script runs the same two
commands with one thread and writes, for each dataset:

- the files those tests compare (``.tsv``, ``.tsv.gz`` and ``.mtx.gz`` outside
  ``intermediary_files/``);
- ``expected_paths.txt``, every path the run creates except the preprocessed
  reads (``sans_bc_*.fq``);
- ``sans_bc_sha256.txt``, the SHA-256 checksums of those preprocessed reads.

A stored file is rewritten only when its decompressed content changes, so
unchanged outputs keep their bytes; rewritten gzip files carry no timestamp or
file name.

Run it from the repository root with the test environment activated, then
check ``git status tests/`` before committing the changes:

    python tests/regenerate_expected.py
"""

import gzip
import hashlib
import os
import subprocess
import sys
import tempfile
from glob import glob

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(TESTS_DIR)
EXAMPLES_DIR = os.path.join(REPO_DIR, "examples")

SAMPLE_TYPES = ["gDNA", "cDNA"]
READ_LAYOUTS = ["single_end", "paired_end"]
COMPARED_SUFFIXES = (".tsv", ".tsv.gz", ".mtx.gz")


def run_pipeline(input_dir: str, config_fpath: str, output_dir: str) -> None:
    """Run preprocess with the placeholder BAM, then count, as the tests do.

    Args:
        input_dir: Folder with the example FASTQ files.
        config_fpath: The example's bcwithqc config.
        output_dir: Folder for all outputs; preprocess and count share it.

    Raises:
        subprocess.CalledProcessError: If either command fails. Its output is
            printed first.
    """
    commands = [
        [
            sys.executable,
            "-m",
            "bcwithqc",
            "preprocess",
            input_dir,
            f"--config={config_fpath}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--output-format-bam",
            "-v",
        ],
        [
            sys.executable,
            "-m",
            "bcwithqc",
            "count",
            output_dir,
            f"--STAR-output-dir={output_dir}",
            f"--config={config_fpath}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--keep-intermediary",
            "-v",
        ],
    ]
    for command in commands:
        # cwd=REPO_DIR makes "python -m bcwithqc" use this checkout.
        result = subprocess.run(
            command, cwd=REPO_DIR, capture_output=True, text=True, check=False
        )
        if result.returncode != 0:
            print(result.stdout, result.stderr, sep="\n", file=sys.stderr)
            result.check_returncode()


def list_output_paths(output_dir: str) -> list[str]:
    """Return every folder and file under ``output_dir`` as a sorted relative path.

    Preprocessed FASTQs (``sans_bc_*.fq``) are left out, as in the tests'
    ``test_file_existence``.
    """
    relative_paths = []
    for root, dirs, files in os.walk(output_dir):
        for name in dirs + files:
            if name.startswith("sans_bc_") and name.endswith(".fq"):
                continue
            relative_paths.append(os.path.relpath(os.path.join(root, name), output_dir))
    return sorted(relative_paths)


def list_compared_files(base_dir: str) -> list[str]:
    """Return the relative paths of the files that ``test_file_contents`` compares."""
    compared = []
    for relative_path in list_output_paths(base_dir):
        if "intermediary_files" in relative_path.split(os.sep):
            continue
        if relative_path.endswith(COMPARED_SUFFIXES) and os.path.isfile(
            os.path.join(base_dir, relative_path)
        ):
            compared.append(relative_path)
    return compared


def read_content(fpath: str) -> bytes:
    """Return a file's bytes, decompressed if it is gzipped."""
    with open(fpath, "rb") as fh:
        data = fh.read()
    return gzip.decompress(data) if fpath.endswith(".gz") else data


def write_content(fpath: str, content: bytes) -> None:
    """Write ``content`` to ``fpath``; gzip it, without a timestamp, for .gz names."""
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "wb") as fh:
        fh.write(gzip.compress(content, mtime=0) if fpath.endswith(".gz") else content)


def sha256_checksums(output_dir: str) -> str:
    """Return ``sha256sum``-style lines for the preprocessed reads of a run.

    count with ``--keep-intermediary`` leaves the ``sans_bc_*.fq`` files in
    ``intermediary_files/``.
    """
    lines = []
    for fpath in sorted(
        glob(os.path.join(output_dir, "intermediary_files", "sans_bc_*.fq"))
    ):
        with open(fpath, "rb") as fh:
            lines.append(
                f"{hashlib.sha256(fh.read()).hexdigest()}  {os.path.basename(fpath)}\n"
            )
    return "".join(lines)


def update_stored_file(
    stored_fpath: str, new_content: bytes, relative_path: str, changes: list[str]
) -> None:
    """Write ``new_content`` to ``stored_fpath`` if it differs; record the change."""
    if not os.path.exists(stored_fpath):
        changes.append(f"added {relative_path}")
    elif read_content(stored_fpath) != new_content:
        changes.append(f"changed {relative_path}")
    else:
        return
    write_content(stored_fpath, new_content)


def update_expected_dir(output_dir: str, expected_dir: str) -> list[str]:
    """Bring ``expected_dir`` in line with a fresh pipeline run in ``output_dir``.

    Args:
        output_dir: Outputs of ``run_pipeline``.
        expected_dir: The stored outputs for the same example, in ``tests/``.

    Returns:
        One line per change, such as ``"changed raw_reads_bc_matrix/matrix.mtx.gz"``.
    """
    changes = []
    new_files = list_compared_files(output_dir)
    for relative_path in new_files:
        update_stored_file(
            os.path.join(expected_dir, relative_path),
            read_content(os.path.join(output_dir, relative_path)),
            relative_path,
            changes,
        )

    for relative_path in list_compared_files(expected_dir):
        if relative_path not in new_files:
            os.remove(os.path.join(expected_dir, relative_path))
            changes.append(f"removed {relative_path}")

    listings = {
        "expected_paths.txt": "".join(
            f"{path}\n" for path in list_output_paths(output_dir)
        ),
        "sans_bc_sha256.txt": sha256_checksums(output_dir),
    }
    for name, text in listings.items():
        update_stored_file(
            os.path.join(expected_dir, name), text.encode(), name, changes
        )
    return changes


def main() -> None:
    """Regenerate the stored outputs of all four example datasets."""
    for sample_type in SAMPLE_TYPES:
        for read_layout in READ_LAYOUTS:
            name = f"{sample_type}_{read_layout}"
            input_dir = os.path.join(
                EXAMPLES_DIR, f"{sample_type}_fastqs_{read_layout}"
            )
            config_fpath = os.path.join(EXAMPLES_DIR, f"{name}_config.json")
            with tempfile.TemporaryDirectory(prefix=f"bcwithqc_{name}_") as tmp_dir:
                output_dir = os.path.join(tmp_dir, "out")
                os.makedirs(output_dir)
                run_pipeline(input_dir, config_fpath, output_dir)
                changes = update_expected_dir(output_dir, os.path.join(TESTS_DIR, name))
            print(f"{name}: {'; '.join(changes) if changes else 'unchanged'}")


if __name__ == "__main__":
    main()
