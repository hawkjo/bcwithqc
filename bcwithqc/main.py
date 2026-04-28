"""
bcwithqc: Process sequencing barcodes and UMIs

Usage:
  bcwithqc count            <fastq_dir> (--STAR-ref-dir=<> | --STAR-output-dir=<>) --config=<> [--output-dir=<>] [--threads=<>] [--keep-intermediary] [--block-type-for-STAR-alignment=<>] [-v | -vv | -vvv]
  bcwithqc preprocess       <fastq_dir> --config=<> [--output-dir=<>] [--threads=<>] [--block-type-for-STAR-alignment=<>] [-v | -vv | -vvv]
  bcwithqc count_matrix     <bcwithqc_bam_file> --output-dir=<> [--threads=<>] [-v | -vv | -vvv]
  bcwithqc simulate_reads   --config=<> --output-dir=<> --nreads=<> [--unique-umis=<>] [--seed=<>] [--error-probability=<>] [--substitution-probability=<>] [--insertion-probability=<>] [--random-tail-length=<>] [-v | -vv | -vvv]

Options:
  --STAR-ref-dir=<>:                  Path to directory with STAR index.
  --STAR-output-dir=<>:                   Path to STAR output directory. All (BAM/SAM) files with the suffix "*Aligned.out.bam" will be processed in lexographic order. 
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
"""
import logging
import os
from docopt import docopt
from .__init__ import __version__
from .config import CommandLineArguments
from .count import process_fastqs, preprocess_fastqs
from .count_matrix import build_count_matrices_from_bam
from .simulate import simulate_reads

def main(**kwargs):
    docopt_args = docopt(__doc__, version=__version__)
    arguments = CommandLineArguments(docopt_args)

    os.makedirs(arguments.output_dir, exist_ok=True)

    log = logging.getLogger()
    log.handlers.clear()
    log.setLevel(arguments.log_level)

    formatter = logging.Formatter("%(asctime)s   %(message)s", "%Y-%m-%d %H:%M:%S")

    # console output
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)

    # per-command log file
    log_fpath = os.path.join(arguments.output_dir, f"bcwithqc_{arguments.command}.log")
    file_handler = logging.FileHandler(log_fpath, mode="w")
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)

    log.debug(docopt_args)
    log.info(f"Writing log to: {log_fpath}")

    commands = {
        'count': process_fastqs,
        'preprocess': preprocess_fastqs,
        'count_matrix': build_count_matrices_from_bam,
        'simulate_reads': simulate_reads
    }

    try:
        commands[arguments.command](arguments)
    except Exception:
        log.exception("bcwithqc crashed with an unhandled exception")
        raise


if __name__ == '__main__':
    main()
