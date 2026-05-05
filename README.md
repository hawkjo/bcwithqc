# bcwithqc
An end-to-end tool for processing fastq files with barcoded sequences into annotated fastq, BAM, and count matrix files. Inspired by single-cell DNA and RNA sequencing data, but much more general.

<p align="center">
  <img src="doc/SDRLoneRanger.png" width=50% alt="The Single-Cell DNA RNA Lone Ranger">
</p>
<p align="center">The Single-Cell DNA RNA Lone Ranger</p>

## Installation

bcwithqc works in Linux, and has been tested on el8 and CentOS 7. 

You can install bcwithqc from github directly into your local (conda or virtual) environment using pip:
```
pip install git+https://github.com/hawkjo/bcwithqc.git
```
This typically takes a few minutes. 

[STAR aligner](https://github.com/alexdobin/STAR)  needs to be installed independently as well.

## Usage

The basic usage for bcwithqc can be displayed at any time via `bcwithqc --help`:
```
Usage:
  bcwithqc preprocess       <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [-v | -vv | -vvv]
  bcwithqc count            <fastq_dir> --STAR-output-dir=<> --config=<> [--output-dir=<>] [--threads=<>] [--keep-intermediary] [-v | -vv | -vvv]
  bcwithqc count_matrix     <bcwithqc_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  bcwithqc simulate_reads   --config=<> --output-dir=<> --nreads=<> [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [--random-tail-length=<>] [-v | -vv | -vvv]

Options:
  --STAR-output-dir=<>:               Path to STAR output directory. All BAM/SAM files with the suffix "*Aligned.out.bam" will be processed in lexicographic order.
  --config=<>:                        Path to JSON configuration.
  --output-dir=<>:                    Path to output directory [default: .].
  --threads=<>:                       Number of threads [default: 1].
  -v:                                 Verbose output.
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
The barcode details are input via a json configuration file, of which standard gDNA and RNA versions can be found in the `examples` folder. 

STAR references need to be prebuilt and their top directory input as a parameter.

### Outputs
The primary outputs from bcwithqc are:
* An annotated BAM file
* A read count matrix
* A UMI count matrix

#### BAM file tags
The BAM file is annotated with custom tags that have been created in the style of current community standards. These are:

| Tag | Meaning |
|--|--|
| CB | Cell barcode |
| CR | Raw, uncorrected cell barcode |
| UB | UMI |
| UR | Raw, uncorrected UMI |
| FL | Combined length of the linker sequences |

The cell barcode tag contains all pieces of the cell barcode, including the sample barcode, concatenated with periods.

## Example Workflow

The `examples` folder contains several small example gDNA and RNA datasets for both paired-end and single-end reads and corresponding example scripts (in the `example scripts` subdirectory) and configuration files. The example script .sh files demonstrate proper syntax for their respective datasets and are runnable directly from within the examples folder. 

Basic principle:
1. Run `bcwithqc preprocess` on your fastq files while providing a `config.json` file and specifying an output directory.
    (Important: The config file specifies which parts of the reads will be kept for aligment with STAR)
2. Run `STAR` (BAM unsorted) on the `sans_bc_*.fq` files in the preprocess output directory while providing a `STAR genome index`.
3. Run `bcwithqc count` on the `sans_bc_*.fq` files while providing the STAR output directory containing the `*Aligned.out.bam` files. 
