import logging
import pysam
from .count import count_matrix

log = logging.getLogger(__name__)
pysam.set_verbosity(0)

def build_count_matrices_from_bam(arguments):
    with pysam.AlignmentFile(arguments.bcwithqc_bam_file) as bam:
        for read in bam.fetch(until_eof=True):
            break
        else:
            raise ValueError(
                f"No reads found in BAM file: {arguments.bcwithqc_bam_file}"
            )

    if read.has_tag("UB"):
        log.info("Detected bam file with UMIs...")
        count_matrix(arguments, arguments.bcwithqc_bam_file)
    else:
        log.info("Detected bam file without UMIs...")
        log.warning("THE NO-UMI build_count_matrix_from_bam BRANCH IS CURRENTLY UNTESTED!")
        gDNA_count_matrix(arguments, arguments.bcwithqc_bam_file)
