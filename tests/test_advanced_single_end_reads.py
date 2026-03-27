import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import nullcontext

import pytest


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

se_test_input = os.path.join(SCRIPT_DIR, "../examples/se_test")
se_test_config = os.path.join(SCRIPT_DIR, "../examples/se_test.json")

# Set this to an existing STAR index if you want the test to run with a prebuilt index
star_index = os.path.join(SCRIPT_DIR, "../examples/se_test_genome_index")

# True = temporary output dir, automatically deleted
# False = persistent output dir for manual inspection
USE_TEMP_OUTPUT = True
USE_STAR_REF = True


def test_single_end_pipeline_runs():
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
        context = tempfile.TemporaryDirectory(prefix="se_test_single_end_")
    else:
        output_dir = os.path.join(SCRIPT_DIR, "se_test")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        command = [
            "python", "-m", "bcwithqc", "count",
            se_test_input,
            f"--config={se_test_config}",
            f"--output-dir={output_dir}",
            "--threads=1",
            "--single-end-reads",
            "--keep-intermediary",
            "--block-type-for-STAR-alignment=constant",
            "-vvv",
        ]
        if USE_STAR_REF:
            # Add STAR index only if you want to test with a prebuilt one
            command.insert(4, f"--STAR-ref-dir={star_index}")

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
            print(f"Output dir: {output_dir}")
        except subprocess.CalledProcessError as e:
            sys.stderr.write("Subprocess failed:\n")
            sys.stderr.write(f"Return code: {e.returncode}\n")
            sys.stderr.write(f"Output dir: {output_dir}\n")
            sys.stderr.write(f"Command: {' '.join(command)}\n")
            sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
            sys.stderr.write(f"STDERR:\n{e.stderr}\n")
            sys.stderr.flush()
            raise


############ ADD ADDITIONAL TEST FUNCTION FOR CHECKING THE CORRECT VALUES IN CREATED FILES HERE. 