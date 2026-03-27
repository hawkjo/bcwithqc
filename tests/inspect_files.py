import gzip
import os
import re
from pathlib import Path
import regex


BASE_DIR = Path("/home/link/John_UMIs/bcwithqc/examples")
GDNA_MOTIF = "GTACTCGCAGTAGTC"


def normalize_header(header: str) -> str:
    """
    Normalize FASTQ header so R1 and R2 from the same cluster compare equal.
    Example:
      @... 1:N:0:INDEX
      @... 2:N:0:INDEX
    becomes:
      @... N:0:INDEX
    """
    header = header.rstrip("\n")
    m = re.match(r"^(.*)\s+([12]):(.*)$", header)
    if not m:
        return header
    return f"{m.group(1)} {m.group(3)}"


def trim_r1_after_motif(seq: str, qual: str, motif: str = GDNA_MOTIF, max_errors: int = 1):
    """
    Keep R1 only up to and including the first fuzzy match of motif.
    Allows up to max_errors total edits.
    If no match is found, keep the full read unchanged.
    """
    pattern = f"({motif}){{e<={max_errors}}}"
    match = regex.search(pattern, seq, regex.BESTMATCH)

    if match is None:
        return seq, qual, False

    end = match.end()
    return seq[:end], qual[:end], True


def merge_fastq_pair(r1_path: Path, r2_path: Path, out_path: Path, trim_gdna_r1: bool = False) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2, open(out_path, "w") as out:
        rec_num = 0

        while True:
            h1 = f1.readline()
            if not h1:
                break

            s1 = f1.readline()
            p1 = f1.readline()
            q1 = f1.readline()

            h2 = f2.readline()
            s2 = f2.readline()
            p2 = f2.readline()
            q2 = f2.readline()

            rec_num += 1

            if not all([s1, p1, q1, h2, s2, p2, q2]):
                raise ValueError(
                    f"Incomplete FASTQ record near record {rec_num}:\n"
                    f"R1={r1_path}\nR2={r2_path}"
                )

            h1s = h1.rstrip("\n")
            s1s = s1.rstrip("\n")
            p1s = p1.rstrip("\n")
            q1s = q1.rstrip("\n")

            h2s = h2.rstrip("\n")
            s2s = s2.rstrip("\n")
            p2s = p2.rstrip("\n")
            q2s = q2.rstrip("\n")

            if not h1s.startswith("@") or not h2s.startswith("@"):
                raise ValueError(f"Bad FASTQ header at record {rec_num}")

            if not p1s.startswith("+") or not p2s.startswith("+"):
                raise ValueError(f"Bad FASTQ plus line at record {rec_num}")

            if len(s1s) != len(q1s):
                raise ValueError(f"R1 seq/qual length mismatch at record {rec_num}")
            if len(s2s) != len(q2s):
                raise ValueError(f"R2 seq/qual length mismatch at record {rec_num}")

            if normalize_header(h1s) != normalize_header(h2s):
                raise ValueError(
                    f"Header mismatch at record {rec_num}\n"
                    f"R1: {h1s}\n"
                    f"R2: {h2s}"
                )

            if trim_gdna_r1:
                s1s, q1s, found = trim_r1_after_motif(s1s, q1s)
                # Uncomment for debugging:
                # if not found:
                #     print(f"Warning: motif not found in {r1_path.name}, record {rec_num}")

            merged_header = h1s
            merged_seq = s1s + s2s
            merged_qual = q1s + q2s

            if len(merged_seq) != len(merged_qual):
                raise ValueError(
                    f"Merged seq/qual length mismatch at record {rec_num}: "
                    f"{len(merged_seq)} vs {len(merged_qual)}"
                )

            out.write(merged_header + "\n")
            out.write(merged_seq + "\n")
            out.write("+\n")
            out.write(merged_qual + "\n")


def find_pairs(input_dir: Path):
    gz_files = sorted(input_dir.glob("*.txt.gz"))
    r1_files = [p for p in gz_files if "_1" in p.name or "R1" in p.name or "1.txt.gz" in p.name]

    pairs = []
    for r1 in r1_files:
        candidates = [
            Path(str(r1).replace("_1", "_2")),
            Path(str(r1).replace("R1", "R2")),
            Path(str(r1).replace("1.txt.gz", "2.txt.gz")),
        ]
        r2 = next((c for c in candidates if c.exists()), None)
        if r2 is None:
            raise FileNotFoundError(f"Could not find matching R2 for {r1}")
        pairs.append((r1, r2))
    return pairs


def process_dataset(dataset_name: str):
    input_dir = BASE_DIR / dataset_name
    output_dir = BASE_DIR / f"{dataset_name}_single_end"

    print(f"\nProcessing {dataset_name}")
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")

    pairs = find_pairs(input_dir)
    if not pairs:
        raise FileNotFoundError(f"No FASTQ pairs found in {input_dir}")

    trim_gdna_r1 = dataset_name == "gDNA_fastqs"

    for r1, r2 in pairs:
        out_name = r1.name.replace(".gz", "")
        out_path = output_dir / out_name
        print(f"Merging:\n  {r1.name}\n  {r2.name}\n  -> {out_path.name}")
        merge_fastq_pair(r1, r2, out_path, trim_gdna_r1=trim_gdna_r1)


if __name__ == "__main__":
    process_dataset("cDNA_fastqs")
    process_dataset("gDNA_fastqs")
    print("\nDone.")