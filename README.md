# bcwithqc
An end-to-end tool for processing fastq files with barcoded sequences into annotated fastq, BAM, and count matrix files. Inspired by single-cell DNA and RNA sequencing data, but much more general.

<p align="center">
  <img src="doc/SDRLoneRanger.png" width=50% alt="The Single-Cell DNA RNA Lone Ranger">
</p>
<p align="center">The Single-Cell DNA RNA Lone Ranger</p>

## Installation

bcwithqc works in Linux, and has been tested on el8 and CentOS 7. It requires Python 3.10 or 3.11 (`biopython==1.79` is pinned).

You can install bcwithqc from github directly into your local (conda or virtual) environment using pip:
```
pip install git+https://github.com/hawkjo/bcwithqc.git
```
This typically takes a few minutes. 

bcwithqc does not align reads itself. Install an aligner separately; the examples and tests use the [STAR aligner](https://github.com/alexdobin/STAR) (CI uses STAR 2.7.10b). To run the test suite, see [Running the tests](#running-the-tests).

## Usage

The basic usage for bcwithqc can be displayed at any time via `bcwithqc --help`:
```
Usage:
  bcwithqc preprocess       <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  bcwithqc count            <fastq_dir> --STAR-output-dir=<> --config=<> --output-dir=<> [--threads=<>] [--keep-intermediary] [-v | -vv | -vvv]
  bcwithqc count_matrix     <bcwithqc_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  bcwithqc simulate_reads   --config=<> --output-dir=<> --nreads=<> [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [--random-tail-length=<>] [-v | -vv | -vvv]

Options:
  --STAR-output-dir=<>:               Path to STAR output directory. All BAM/SAM files with the suffix "*Aligned.out.bam" will be processed in lexicographic order.
  --config=<>:                        Path to JSON configuration.
  --output-dir=<>:                    Path to output directory [default: .]. Required for count, which removes
                                        (or, with --keep-intermediary, moves) everything in it except its final outputs.
  --threads=<>:                       Number of threads [default: 1].
  -v:                                 Verbosity of output. Can be -v, -vv, or -vvv, with more v's leading to more verbose output.
  --nreads=<>:                        Number of reads to simulate.
  --unique-umis=<>:                   Fraction of all reads that have unique UMIs [default: 0.5].
  --seed=<>:                          Random seed [default: 42].
  --error-probability=<>:             Probability of an error occurring per base [default: 0.1]. Set to a negative number to
                                        always introduce as many errors as allowed by the configuration.
  --substitution-probability=<>:      Probability of generating a substitution as opposed to an indel [default: 0.7].
  --insertion-probability=<>:         Probability of generating an insertion as opposed to a deletion when generating an indel [default: 0.5].
  --random-tail-length=<>:            Mean (poisson) length of the random nucleotide tail [default: 20]. Set to a negative number to
                                        generate reads without tails.
  --keep-intermediary                 Keep intermediary files instead of deleting them [default: False].
  -h --help                           Show this screen.
  --version                           Show version.

Commands:
  preprocess       Preprocess files such that STAR can be run on the output.
  count            Process and count input files using an existing STAR output directory.
  count_matrix     Build a count matrix (or matrices) from an existing bam file.
  simulate_reads   Generate synthetic sequencing reads given a barcode configuration.
```
The barcode details are input via a json configuration file, of which standard gDNA and RNA versions can be found in the `examples` folder. [`doc/config_format.md`](doc/config_format.md) describes every config key, its default and its meaning.

Alignment is done by the user between `preprocess` and `count` (see Basic Workflow).

## Example Workflow

The `examples` folder contains several small example gDNA and RNA datasets for both paired-end and single-end reads and corresponding example scripts (in the `example scripts` subdirectory) and configuration `.json` files. The example script `.sh` files demonstrate proper syntax for their respective datasets and are runnable directly from within the examples folder. 

### Basic Workflow:
1. Run `bcwithqc preprocess` on your fastq files while providing a `config.json` file and specifying an output directory.
    (Important: The config file specifies which parts of the reads will be kept for alignment, via `keep_nonbarcode` and `keepblock`.)
    `single_end_reads` in the config selects single-end or paired-end processing.
2. Align the `sans_bc_*.fq` files from the preprocess output directory yourself, producing one unsorted BAM per input FASTQ (pair), named `*Aligned.out.bam` and placed in one directory. Counting needs the `GX`/`GN` gene tags. The tests use:
    ```
    STAR --genomeDir <index> --readFilesIn <sans_bc fastq(s)> --outFileNamePrefix <star_dir>/<sample>_ \
         --outFilterMultimapNmax 1 --outSAMtype BAM Unsorted --outSAMattributes NH HI AS nM GX GN
    ```
3. Run `bcwithqc count <preprocess output dir> --STAR-output-dir=<star_dir> --config=<same config>` with the same `--output-dir` as `preprocess`. BAMs and tag files are matched in lexicographic order.

## Outputs
The primary outputs from bcwithqc are:
* An annotated BAM file
* A read count matrix
* A UMI count matrix (only when the config contains a `UMI` block)

```text
<sample_output_folder>/
├── intermediary_files/          # optional
├── logs/
├── QC_metrics/
├── raw_reads_bc_matrix/
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
├── raw_umis_bc_matrix/         # only with UMI blocks
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
├── with_bc_umi.sorted.bam
└── with_bc_umi.sorted.bam.bai
```

### BAM file tags
The BAM file `with_bc_umi.sorted.bam` is annotated with custom tags that have been created in the style of current community standards. These are:

| Tag  | Meaning                                 |
| ---- | --------------------------------------- |
| `CB` | Cell barcode                            |
| `CR` | Raw, uncorrected cell barcode           |
| `UB` | Corrected UMI (only with UMI blocks)    |
| `UR` | Raw, uncorrected UMI                    |
| `FL` | Combined length of the linker sequences |


The cell barcode tag contains all pieces of the cell barcode, including the sample barcode, concatenated with periods.

### Read/UMI count matrices
The count matrices `matrix.mtx.gz` are written in sparse Matrix Market format.

| Matrix component | Meaning | Defined by |
|---|---|---|
| Rows | Features | `features.tsv.gz` |
| Columns | Barcodes | `barcodes.tsv.gz` |
| Values | Read or UMI counts | `matrix.mtx.gz` |

`matrix[i, j]` is the count for feature `i` and barcode `j`.

### QC_metrics 
QC metrics contains:
1. Filtered FASTQ files containing reads that were excluded from further processing because:
- A sequence was decoded to multiple possible barcodes                   -> `*_ambiguous_reads.fq`
- A sequence could not be decoded to any barcode                         -> `*_no_match_reads.fq`
- A sequence was decoded, but the overall score of the read was too low  -> `*_threshold_fail_reads.fq`

2. QC tables and graphics, both as summary and detailed versions

- `bcs_summary.tsv`  
  Summary per barcode block: how many barcode observations in each block were exact matches, corrected, below threshold despite being corrected, ambiguous, or not matched.

- `reads_summary.tsv`  
  Summary per read/read pair: how many reads in total were exact matches, corrected, below threshold despite being corrected, ambiguous, or not matched.

- `reads_and_blocks_summary.png`  
  Combined summary figure: read-level summary on top and barcode-block-level summary below.

- `bcs.tsv` and `barcodes_<read>_<block>.png`  
  Detailed barcode-level QC: the same status information as above, but for each individual expected/provided barcode.  
  The corresponding barcode plots show the non-normalized counts on top and the normalized percentages below.

- `reads.tsv`  
  Detailed read-level QC: one row per read/read pair, including the total read status and the status of each barcode block.

## Running the tests

The tests use pytest. From the repository root, create a Python 3.10 or 3.11 environment with the `test` extra:
```
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```
With [uv](https://docs.astral.sh/uv/): `uv venv -p 3.11 .venv`, `source .venv/bin/activate`, then `uv pip install -e ".[test]"`.

The pipeline tests also run the STAR aligner (bcwithqc itself never does), so STAR must be on your `PATH`. CI uses STAR 2.7.10b. To get the same version under `~/opt`:
```
mkdir -p ~/opt
curl -L https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz | tar -xzf - -C ~/opt
```
- Linux: build it as CI does, `make -C ~/opt/STAR-2.7.10b/source STAR`; the program is then in `~/opt/STAR-2.7.10b/source`.
- macOS: use the included binary in `~/opt/STAR-2.7.10b/bin/MacOSX_x86_64` (it runs under Rosetta on Apple Silicon). Homebrew's `rna-star` (2.7.11b) fails with the test genome indexes.

Add that directory to `PATH` in your shell profile and open a new terminal, for example on macOS:
```
echo 'export PATH="$HOME/opt/STAR-2.7.10b/bin/MacOSX_x86_64:$PATH"' >> ~/.zshrc
```
Adding it by hand after activating the environment does not last: re-activating a virtual environment restores the `PATH` it saved when it was first activated.

Then, with the environment activated:
```
pytest
```
The full suite takes about 7 minutes. Other ways to run it:
- one file: `pytest tests/test_15_count_matrix_umis.py`
- stop at the first failure: `pytest -x`
- as CI does: `pytest --capture=tee-sys -vv`

Run `pytest` with the environment activated rather than as `.venv/bin/pytest`: the pipeline tests start `python -m bcwithqc` in a subprocess, which has to find the environment's Python.

On macOS, expect 8 failures, all `test_file_contents` in `test_05` and `test_06`. The stored outputs they compare against were made on Linux, and STAR's macOS build assigns gene tags (`GX`) to different reads from identical alignments. CI runs these comparisons on Linux. To leave them out locally, run `pytest -k "not test_file_contents"`. `test_00_cpu_avx512_compatibility` is skipped where `/proc/cpuinfo` is missing (macOS).
