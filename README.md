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
  bcwithqc count            <fastq_dir> (--STAR-ref-dir=<> | --STAR-output=<>...) --config=<> [--output-dir=<>] [--threads=<>] [--keep-intermediary] [--single-end-reads] [--block-type-for-STAR-alignment=<>] [-v | -vv | -vvv]
  bcwithqc preprocess       <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [--single-end-reads] [--block-type-for-STAR-alignment=<>] [-v | -vv | -vvv]
  bcwithqc count_matrix     <bcwithqc_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  bcwithqc simulate_reads   --config=<> --output-dir=<> --nreads=<> [--single-end-reads] [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [--random-tail-length=<>] [-v | -vv | -vvv]

Options:
  --STAR-ref-dir=<>:                  Path to directory with STAR index.
  --STAR-output=<>:                   Path to STAR output file (BAM/SAM). Can be repeated multiple times,
                                        in which case the order must correspond to the lexicographic ordering
                                        of paired FASTQ files in <fastq_dir>.
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
  --single-end-reads:                 Option for using single end instead of paired end reads [default: False]. 
  --block-type-for-STAR-alignment=<>: What part of the read is retained for STAR alignment [default: undefined_only] 
                                        'undefined_only' concatenates all parts that are not specified in the config.json
                                        'constant'  concatenates all 'constantRegion' blocktypes and unspecified parts
                                        'constant_mask' retains 'constantRegion' blocktypes and unspecified parts and replaces the rest with Ns
  -h --help                           Show this screen.
  --version                           Show version.

Commands:
  preprocess       Preprocess files such that STAR can be run on the output.
  count            Process and count input files.
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

## Examples

The `examples` folder contains small example gDNA and RNA datasets and corresponding example scripts and configuration files. The `cmd_gDNA_json.sh` and `cmd_cDNA_json.sh` files demonstrate proper syntax for their respective datasets and are runnable directly from within the examples folder. They each take about 20 seconds to run.
