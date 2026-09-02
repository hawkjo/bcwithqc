import re
from collections import Counter

ERROR_SEGMENT_PATTERN = re.compile(r"S(\d+)D(\d+)I(\d+)")

def update_tsv_to_new_format(input_path, output_path):
    combo_counter = Counter()

    # Initialize counters for the new summary and indel issues
    should_reject_2_plus_errors_in_any_block = 0
    first_block_1del_second_block_1ins = 0
    first_block_1ins_second_block_1del = 0
    should_accept_1_minus_errors_in_each_block = 0  # New counter for the 1-or-less errors condition
    first_block_1del_second_block_1error = 0
    first_block_1ins_second_block_1error = 0
    second_block_1del_first_block_1error = 0
    second_block_1ins_first_block_1error = 0

    # Read the input file
    with open(input_path, "r") as fh:
        lines = [line.rstrip("\n") for line in fh]

    # Read the sections and find the count/combination part
    summary_section_end = lines.index("")  # Find the blank line between summary and combinations
    combinations_start = summary_section_end + 2  # The combinations section starts after "count\tcombination"

    # Read the combination lines to preserve them
    combination_lines = lines[combinations_start:]

    # Process combinations and count occurrences
    for line in combination_lines:
        count, combo = line.split("\t")
        count = int(count)
        combo_counter[tuple(combo.split(";"))] += count  # Store the combination and count

    # Compute the summary statistics
    for line in combination_lines:
        count, combo = line.split("\t")
        count = int(count)
        combos = combo.split(";")

        # Parse the error patterns in both blocks
        s1, d1, i1 = map(int, re.match(r"S(\d+)D(\d+)I(\d+)", combos[0]).groups())
        s2, d2, i2 = map(int, re.match(r"S(\d+)D(\d+)I(\d+)", combos[1]).groups())

        # Update counters based on the new logic
        if (s1 + d1 + i1) >= 2 or (s2 + d2 + i2) >= 2:
            should_reject_2_plus_errors_in_any_block += count

        if (s1 + d1 + i1 <= 1) and (s2 + d2 + i2 <= 1):
            should_accept_1_minus_errors_in_each_block += count

        if d1 == 1 and s1 == 0 and i1 == 0 and i2 == 1 and s2 == 0 and d2 == 0:
            first_block_1del_second_block_1ins += count
        if i1 == 1 and s1 == 0 and d1 == 0 and d2 == 1 and s2 == 0 and i2 == 0:
            first_block_1ins_second_block_1del += count

        if d1 == 1 and s1 == 0 and i1 == 0 and (s2 + d2 + i2 <= 1):
            first_block_1del_second_block_1error += count
        if i1 == 1 and s1 == 0 and d1 == 0 and (s2 + d2 + i2 <= 1):
            first_block_1ins_second_block_1error += count

        if d2 == 1 and s2 == 0 and i2 == 0 and (s1 + d1 + i1 <= 1):
            second_block_1del_first_block_1error += count
        if i2 == 1 and s2 == 0 and d2 == 0 and (s1 + d1 + i1 <= 1):
            second_block_1ins_first_block_1error += count

    # Write the output file with the updated summary and indel issues
    with open(output_path, "w") as out:
        # Write the new summary section
        out.write("summary\tcount\n")
        out.write(f"Should Reject: 2+ errors in any block\t{should_reject_2_plus_errors_in_any_block}\n")
        out.write(f"Should Accept: 1- errors in each block\t{should_accept_1_minus_errors_in_each_block}\n")
        out.write("\n")

        # Write the updated indel issues section
        out.write("Indel Issues\n")
        out.write(f"First block: 1D; Second block: 1I\t{first_block_1del_second_block_1ins}\n")
        out.write(f"first block: 1I; Second block: 1D\t{first_block_1ins_second_block_1del}\n")
        out.write(f"First block: 1D; Second block: 1 or fewer errors\t{first_block_1del_second_block_1error}\n")
        out.write(f"First block: 1I; Second block: 1 or fewer errors\t{first_block_1ins_second_block_1error}\n")
        out.write(f"Second block: 1D; First block: 1 or fewer errors\t{second_block_1del_first_block_1error}\n")
        out.write(f"Second block: 1I; First block: 1 or fewer errors\t{second_block_1ins_first_block_1error}\n")
        out.write("\n")

        # Write the combinations section (preserving the existing combinations)
        out.write("count\tcombination\n")
        for combo, count in combo_counter.most_common():
            out.write(f"{count}\t{';'.join(combo)}\n")


# Usage: Call this function for each old TSV file
update_tsv_to_new_format("/home/link/John_UMIs/bcwithqc/tests/simulate_se_mini/count_output/simulate_se_mini_before.tsv", "/home/link/John_UMIs/bcwithqc/tests/simulate_se_mini/count_output/simulate_se_mini_before.tsv")
update_tsv_to_new_format("/home/link/John_UMIs/bcwithqc/tests/simulate_se_mini/count_output/simulate_se_mini_after.tsv", "/home/link/John_UMIs/bcwithqc/tests/simulate_se_mini/count_output/simulate_se_mini_after.tsv")