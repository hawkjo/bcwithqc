import os
import sys
import shutil
import subprocess
import pytest

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

gDNA_input = os.path.join(SCRIPT_DIR, "../examples/gDNA_fastqs_single_end")
cDNA_input = os.path.join(SCRIPT_DIR, "../examples/cDNA_fastqs_single_end")

star_index = os.path.join(SCRIPT_DIR, "../examples/SDR001_REF_index")
gDNA_config = os.path.join(SCRIPT_DIR, "../examples/gDNA.json")
cDNA_config = os.path.join(SCRIPT_DIR, "../examples/cDNA.json")


@pytest.mark.parametrize(
    "sample_type,input_dir,config",
    [
        ("gDNA", gDNA_input, gDNA_config),
        ("cDNA", cDNA_input, cDNA_config),
    ],
)
def test_single_end_pipeline_runs(sample_type, input_dir, config):
    output_dir = os.path.join(SCRIPT_DIR, f"{sample_type}_single_end")

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

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

    command = [
        "python", "-m", "bcwithqc", "count",
        input_dir,
        f"--STAR-ref-dir={star_index}",
        f"--config={config}",
        f"--output-dir={output_dir}",
        "--threads=1",
        "--keep-intermediary",
        "--single-end-reads",
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
        sys.stderr.write(f"STDOUT:\n{e.stdout}\n")
        sys.stderr.write(f"STDERR:\n{e.stderr}\n")
        sys.stderr.flush()
        raise