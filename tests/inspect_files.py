#!/usr/bin/env python3

import re
from pathlib import Path

INPUT_FASTQ = Path(
    "/home/link/John_UMIs/bcwithqc/tests/simulate_se/simulate_se_1k_reads_default_error/"
    "count_output/QC_metrics/simulate_se_mini_config_no_match_reads.fq"
)

OUTPUT_FASTQ = Path(
    "/home/link/John_UMIs/bcwithqc/examples/se_mini/se_mini.txt"
)

ERROR_SEGMENT_PATTERN = re.compile(r"S(\d+)D(\d+)I(\d+)")


def parse_error_segments(read_name):
    """
    Extract all SxDyIz error segments from a FASTQ read name.

    Example:
    @150_GGCGACTCA_S1D0I1_GGGTCAGTACGTACGAGTCCC_S1D0I0_CCCTTT_randTail26 ...
    returns:
    [(1, 0, 1), (1, 0, 0)]
    """
    return [
        tuple(map(int, match.groups()))
        for match in ERROR_SEGMENT_PATTERN.finditer(read_name)
    ]


def keep_read(read_name):
    error_segments = parse_error_segments(read_name)

    if len(error_segments) < 2:
        return False

    first_block = error_segments[0]
    second_block = error_segments[1]

    first_block_matches = first_block == (0, 0, 1)
    second_block_has_at_most_one_error = sum(second_block) <= 1

    return first_block_matches and second_block_has_at_most_one_error


def main():
    OUTPUT_FASTQ.parent.mkdir(parents=True, exist_ok=True)

    total_reads = 0
    kept_reads = 0

    with INPUT_FASTQ.open("r") as infile, OUTPUT_FASTQ.open("w") as outfile:
        while True:
            header = infile.readline()
            if not header:
                break

            seq = infile.readline()
            plus = infile.readline()
            qual = infile.readline()

            if not qual:
                raise ValueError("Input FASTQ appears truncated.")

            total_reads += 1

            if keep_read(header.strip()):
                outfile.write(header)
                outfile.write(seq)
                outfile.write(plus)
                outfile.write(qual)
                kept_reads += 1

    print(f"Total reads: {total_reads}")
    print(f"Kept reads:  {kept_reads}")
    print(f"Wrote:       {OUTPUT_FASTQ}")


if __name__ == "__main__":
    main()