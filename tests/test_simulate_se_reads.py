import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import nullcontext

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

se_test_config = os.path.join(SCRIPT_DIR, "../examples/se_test.json")

# Set to True to use a temp directory that is deleted automatically.
# Set to False to write into the config directory under "se_test_lc" and keep the output.
USE_TEMP_OUTPUT = True


def test_single_end_simulation_runs():
    config = se_test_config
    config_dir = os.path.dirname(os.path.abspath(config))

    if USE_TEMP_OUTPUT:
        context = tempfile.TemporaryDirectory(prefix="se_test_lc_")
    else:
        output_dir = os.path.join(config_dir, "se_test_lc")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        command = [
            "python", "-m", "bcwithqc", "simulate_reads",
            f"--config={config}",
            f"--output-dir={output_dir}",
            "--nreads=10000",
            "--error-probability=0.3",
            "--substitution-probability=0.3",
            "--single-end-reads",
            "-vvv",
        ]

        try:
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
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

        config_stem = os.path.splitext(os.path.basename(config))[0]
        expected_fastq = os.path.join(output_dir, f"{config_stem}.txt.gz")

        assert os.path.isdir(output_dir)
        assert os.path.isfile(expected_fastq), f"Expected FASTQ not found: {expected_fastq}"