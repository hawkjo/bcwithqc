import gzip
import os
import re
from pathlib import Path
import regex
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

fastq_gz = os.path.join(SCRIPT_DIR, "../examples/se_test_lc/se_test_lc.txt.gz")
out_tsv = os.path.join(SCRIPT_DIR, "../test/se_test_lc/se_test_lc_before.tsv")
pattern = re.compile(r"S\d+D\d+I\d+")

combo_counter = Counter()

with gzip.open(fastq_gz, "rt") as fh:
    for i, line in enumerate(fh):
        if i % 4 == 0:
            header = line.strip()
            segments = tuple(pattern.findall(header))
            combo_counter[segments] += 1

with open(out_tsv, "w") as out:
    out.write("count\tcombination\n")
    for combo, count in combo_counter.most_common():
        out.write(f"{count}\t{';'.join(combo)}\n")