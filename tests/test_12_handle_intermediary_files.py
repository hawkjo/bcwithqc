import sys
import os
import shutil
import tempfile
from pathlib import Path
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from bcwithqc.count import handle_intermediary_files  # Adjust this import to match your structure


class DummyArgs:
    def __init__(self, output_dir, keep_intermediary_files):
        self.output_dir = output_dir
        self.keep_intermediary_files = keep_intermediary_files


@pytest.mark.parametrize("keep_files", [False, True])
def test_handle_intermediary_files(keep_files):
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        args = DummyArgs(output_dir=tmpdir, keep_intermediary_files=keep_files)

        # Simulate final output files (should be preserved)
        final_bam = output_dir / "with_bc_umi.sorted.bam"
        final_bam.touch()
        (output_dir / "with_bc_umi.sorted.bam.bai").touch()
        (output_dir / "raw_reads_bc_matrix").mkdir()
        (output_dir / "raw_umis_bc_matrix").mkdir()

        # Simulate intermediary files (should be deleted or moved)
        intermediary_file = output_dir / "tempfile.tmp"
        intermediary_folder = output_dir / "tempdir"
        intermediary_file.touch()
        intermediary_folder.mkdir()

        # Run the function
        handle_intermediary_files(args, str(final_bam))

        # Check preservation
        assert final_bam.exists()
        assert (output_dir / "with_bc_umi.sorted.bam.bai").exists()
        assert (output_dir / "raw_reads_bc_matrix").exists()
        assert (output_dir / "raw_umis_bc_matrix").exists()

        if keep_files:
            moved_file = output_dir / "intermediary_files" / "tempfile.tmp"
            moved_dir = output_dir / "intermediary_files" / "tempdir"
            assert moved_file.exists()
            assert moved_dir.exists()
        else:
            assert not intermediary_file.exists()
            assert not intermediary_folder.exists()
