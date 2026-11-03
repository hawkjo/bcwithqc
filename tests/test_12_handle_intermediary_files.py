import sys
import os
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from bcwithqc.count import handle_intermediary_files


class DummyArgs:
    def __init__(self, output_dir, keep_intermediary_files):
        self.output_dir = output_dir
        self.keep_intermediary_files = keep_intermediary_files


@pytest.mark.parametrize("keep_files", [False, True])
def test_handle_intermediary_files(keep_files):
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir)
        args = DummyArgs(output_dir=tmpdir, keep_intermediary_files=keep_files)

        # Simulate final output files/directories, which should be preserved
        final_bam = output_dir / "with_bc_umi.sorted.bam"
        final_bam.touch()
        (output_dir / "with_bc_umi.sorted.bam.bai").touch()
        (output_dir / "raw_reads_bc_matrix").mkdir()
        (output_dir / "raw_umis_bc_matrix").mkdir()
        (output_dir / "QC_metrics").mkdir()

        # Simulate intermediary files/directories
        intermediary_file = output_dir / "tempfile.tmp"
        intermediary_folder = output_dir / "tempdir"
        intermediary_file.touch()
        intermediary_folder.mkdir()

        # Simulate log/output files
        log_file = output_dir / "bcwithqc.log"
        out_file = output_dir / "Log.final.out"
        already_prefixed_out_file = output_dir / "STAR_Log.progress.out"

        log_file.touch()
        out_file.touch()
        already_prefixed_out_file.touch()

        # Run the function
        handle_intermediary_files(args, str(final_bam))

        # Final outputs should still exist
        assert final_bam.exists()
        assert (output_dir / "with_bc_umi.sorted.bam.bai").exists()
        assert (output_dir / "raw_reads_bc_matrix").exists()
        assert (output_dir / "raw_umis_bc_matrix").exists()
        assert (output_dir / "QC_metrics").exists()

        # logs/ should always exist and be preserved
        logs_dir = output_dir / "logs"
        assert logs_dir.exists()
        assert logs_dir.is_dir()

        # .log files should be moved to logs/
        assert not log_file.exists()
        assert (logs_dir / "bcwithqc.log").exists()

        # .out files should be moved to logs/ and renamed with STAR_ prefix
        assert not out_file.exists()
        assert (logs_dir / "STAR_Log.final.out").exists()

        # Already prefixed .out files should not get STAR_STAR_ prefix
        assert not already_prefixed_out_file.exists()
        assert (logs_dir / "STAR_Log.progress.out").exists()
        assert not (logs_dir / "STAR_STAR_Log.progress.out").exists()

        if keep_files:
            intermediary_dir = output_dir / "intermediary_files"

            moved_file = intermediary_dir / "tempfile.tmp"
            moved_dir = intermediary_dir / "tempdir"

            assert moved_file.exists()
            assert moved_dir.exists()

            # Log/output files should not be moved into intermediary_files/
            assert not (intermediary_dir / "bcwithqc.log").exists()
            assert not (intermediary_dir / "Log.final.out").exists()
            assert not (intermediary_dir / "STAR_Log.final.out").exists()

        else:
            # Regular intermediary files should be deleted
            assert not intermediary_file.exists()
            assert not intermediary_folder.exists()

            # No intermediary_files/ directory should be created in delete mode
            assert not (output_dir / "intermediary_files").exists()