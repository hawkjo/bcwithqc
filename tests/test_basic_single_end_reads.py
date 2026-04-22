import os
import sys
import shutil
import subprocess
import tempfile
from contextlib import nullcontext

import pytest
import json


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_single_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_single_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA.json")

# Set to True to use a temp directory that is deleted automatically.
# Set to False to write into SCRIPT_DIR/<sample_type>_single_end and keep the output.
USE_TEMP_OUTPUT = True


@pytest.mark.parametrize(
    "sample_type,input_dir,config",
    [
        ("gDNA", gDNA_input, gDNA_config),
        ("cDNA", cDNA_input, cDNA_config),
    ],
)
def test_single_end_pipeline_runs(sample_type, input_dir, config):
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
        context = tempfile.TemporaryDirectory(prefix=f"{sample_type}_single_end_basic_")
    else:
        output_dir = os.path.join(SCRIPT_DIR, f"{sample_type}_single_end")
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        context = nullcontext(output_dir)

    with context as output_dir:
        config_to_use = config
        temp_config = None
        
        # We need to edit the cDNA config, and change keep_nonbarcode:false to true, for single end reads, otherwise everything is discarded.
        try:
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
                f"--output-dir={output_dir}",
                "--threads=1",
                "--keep-intermediary",
                "--single-end-reads",
                "-vvv",
            ]

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

        finally:
            if temp_config is not None:
                try:
                    os.remove(temp_config)
                except FileNotFoundError:
                    pass