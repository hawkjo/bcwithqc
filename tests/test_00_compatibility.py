import os
import sys
import subprocess
import pytest


def run_python_snippet(snippet):
    """
    Run a Python snippet in a subprocess so that SIGILL/SIGSEGV does not kill pytest itself.
    """
    return subprocess.run(
        [sys.executable, "-X", "faulthandler", "-c", snippet],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def get_cpu_flags():
    """
    Return CPU flags from /proc/cpuinfo if available.
    """
    cpuinfo_path = "/proc/cpuinfo"

    if not os.path.exists(cpuinfo_path):
        return set()

    with open(cpuinfo_path) as f:
        for line in f:
            if line.startswith("flags"):
                return set(line.strip().split(":")[1].split())

    return set()


def running_on_github_actions():
    return os.environ.get("GITHUB_ACTIONS") == "true"


def env_var_is_true(name):
    return os.environ.get(name, "").lower() in {"1", "true", "yes", "y"}


def env_var_is_false(name):
    return os.environ.get(name, "").lower() in {"0", "false", "no", "n"}


def require_avx512():
    """
    Decide whether this environment should strictly require AVX-512.

    Default behavior:
    - GitHub Actions: do not require AVX-512
    - local/HPC: require AVX-512

    Overrides:
    - BCWITHQC_REQUIRE_AVX512=1 forces requirement
    - BCWITHQC_REQUIRE_AVX512=0 disables requirement
    """
    if env_var_is_true("BCWITHQC_REQUIRE_AVX512"):
        return True

    if env_var_is_false("BCWITHQC_REQUIRE_AVX512"):
        return False

    if running_on_github_actions():
        return False

    return True


def test_00_cpu_avx512_compatibility():
    """
    Check CPU compatibility.

    In local/HPC environments we require AVX-512 by default because pywfa may
    crash with SIGILL if it was compiled with AVX-512 instructions but runs on
    a node without AVX-512.

    In GitHub Actions, AVX-512 is not guaranteed, so we only require AVX2 and
    let the subprocess import/runtime tests catch real crashes.
    """
    flags = get_cpu_flags()

    if not flags:
        pytest.skip("Could not read CPU flags from /proc/cpuinfo")

    assert "avx2" in flags, (
        "This node does not report AVX2 support. "
        "This may be too old for some compiled dependencies. "
        f"Node: {subprocess.getoutput('hostname')}"
    )

    if not require_avx512():
        if "avx512f" not in flags:
            pytest.skip(
                "AVX-512 not available, but this environment does not require it. "
                "This is expected on GitHub Actions. "
                f"Node: {subprocess.getoutput('hostname')}"
            )
        return

    assert "avx512f" in flags, (
        "This node does not support AVX-512, but this environment requires it. "
        "If pywfa was compiled with AVX-512 instructions, it may crash with SIGILL. "
        f"Node: {subprocess.getoutput('hostname')}"
    )


@pytest.mark.parametrize(
    "module_name",
    [
        "numpy",
        "scipy",
        "scipy.sparse",
        "scipy.sparse.csgraph",
        "pysam",
        "Bio",
        "Bio.SeqIO",
        "Bio.Align",
        "regex",
        "matplotlib.pyplot",
        "pywfa",
        "freebarcodes.decode",
        "freebarcodes.editmeasures",
        "bcwithqc",
        "bcwithqc.misc",
        "bcwithqc.bc_decoders",
        "bcwithqc.bc_lookup",
        "bcwithqc.bc_aligner",
        "bcwithqc.bc_parser",
        "bcwithqc.count",
        "bcwithqc.count_matrix",
        "bcwithqc.simulate",
        "bcwithqc.qc_metrics",
    ],
)
def test_01_imports_do_not_crash(module_name):
    """
    Import each relevant module in a subprocess.

    This catches native import crashes such as SIGILL without killing pytest.
    """
    result = run_python_snippet(
        f"import {module_name}; print('{module_name} ok')"
    )

    assert result.returncode == 0, (
        f"Import failed or crashed for module: {module_name}\n"
        f"Return code: {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )


def test_02_pywfa_aligner_construction_does_not_crash():
    """
    Test pywfa construction separately from import.

    Some CPU-instruction problems only appear when the compiled code is used,
    not when it is imported.
    """
    snippet = """
from pywfa import WavefrontAligner

print("before WavefrontAligner")
aligner = WavefrontAligner()
print("after WavefrontAligner")
"""

    result = run_python_snippet(snippet)

    assert result.returncode == 0, (
        "pywfa WavefrontAligner construction failed or crashed.\n"
        f"Return code: {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )


def test_03_distance_thresh_levenshtein_does_not_crash():
    """
    This is the minimal check for the crash seen during preprocessing.

    The previous failure happened at:
        DistanceThresh('levenshtein', 2)
    """
    snippet = """
from bcwithqc.misc import DistanceThresh

print("before DistanceThresh")
distfun = DistanceThresh("levenshtein", 2)
print("after DistanceThresh")
print(distfun("AAAA", "AAAT"))
"""

    result = run_python_snippet(snippet)

    assert result.returncode == 0, (
        "DistanceThresh('levenshtein', 2) failed or crashed.\n"
        "This likely indicates an incompatible compiled backend, e.g. pywfa "
        "compiled with AVX-512 on a node without AVX-512.\n"
        f"Return code: {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )