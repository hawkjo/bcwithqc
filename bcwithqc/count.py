import logging
import os
import pysam
import itertools
import numpy as np
import matplotlib.pyplot as plt
import pickle
import subprocess
import shutil
import scipy
import warnings
from . import misc
from . import qc_metrics
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from contextlib import nullcontext
from glob import glob
from multiprocessing import Pool
from threading import BoundedSemaphore
from scipy.sparse import lil_matrix
from .bc_aligner import CustomBCAligner
from .umi import get_umi_maps_from_bam_file


log = logging.getLogger(__name__)
pysam.set_verbosity(0)


n_first_seqs = 10000  # n seqs for finding score threshold

def preprocess_fastqs(arguments):
    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)

    # For single end reads
    if arguments.single_end_reads:
        single_fpaths = misc.find_single_fastqs_in_dir(arguments.fastq_dir)
        if arguments.config["unknown_read_orientation"]:
            warnings.warn(
                "`unknown_read_orientation` has no effect for single-end reads; ignoring it.",
                category=UserWarning)
        if misc.config_has_barcodes_on_both_reads(arguments.config):
            warnings.warn(
                "`Barcodes on both reads` option has no effect for single-end reads; ignoring it.",
                category=UserWarning)  

        log.info('Files to process:')
        for fpath in single_fpaths:
            log.info(f'  {fpath}')

        single_align_fq_and_tags_fpaths = []

        process_fastqs_func = (
            serial_process_fastqs_single_end
            if arguments.threads == 1
            else parallel_process_fastqs_single_end
        )

        all_reads_qcs = []
        for fpath in single_fpaths:
            bc_fq1_fpath = fpath
            bc_fq2_fpath = None

            sans_bc_fq1_fpath = os.path.join(
                arguments.output_dir,
                f'sans_bc_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.fq'
            )

            tags1_fpath = os.path.join(
                arguments.output_dir,
                f'rec_names_and_tags1_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.txt'
            )

            log.info('Processing file:')
            log.info(f'  barcode fastq: {bc_fq1_fpath}')
            log.info(f'  sans-bc fastq: {sans_bc_fq1_fpath}')
            log.info(f'  tags file:     {tags1_fpath}')

            if any(os.path.exists(fpath) for fpath in [sans_bc_fq1_fpath, tags1_fpath]):
                raise FileExistsError('Partial results found: Remove partial results and restart.')

            file_read_qcs = process_fastqs_func(
                arguments,
                bc_fq1_fpath,
                sans_bc_fq1_fpath,
                tags1_fpath,
            )
            all_reads_qcs.extend(file_read_qcs)
            single_align_fq_and_tags_fpaths.append((sans_bc_fq1_fpath, tags1_fpath))

        qc_metrics.generate_qc_metrics(arguments, all_reads_qcs)
        qc_metrics.generate_qc_metrics_reads(arguments, all_reads_qcs)
        return single_align_fq_and_tags_fpaths

    # For paired reads
    else:
        paired_fpaths = misc.find_paired_fastqs_in_dir(arguments.fastq_dir)

        if arguments.config["unknown_read_orientation"]:
            paired_fpaths = misc.fix_unknown_read_orientation(arguments, paired_fpaths)

        log.info('Files to process:')
        for i, (fpath1, fpath2) in enumerate(paired_fpaths):
            log.info(f'  {fpath1}')
            log.info(f'  {fpath2}')
            if i < len(paired_fpaths)-1:
                log.info('  -')

        paired_align_fqs_and_tags_fpaths = []
        namepairidxs = []

        process_fastqs_func = serial_process_fastqs if arguments.threads == 1 else parallel_process_fastqs
        bcs_on_both_reads = misc.config_has_barcodes_on_both_reads(arguments.config)

        if not bcs_on_both_reads:
            bc_fq_idx, paired_fq_idx = misc.determine_bc_and_paired_fastq_idxs(paired_fpaths, arguments.config["barcode_struct_r1"]["blocks"], arguments.config["unknown_read_orientation"])
            log.info(f'Detected barcodes in read{bc_fq_idx+1} files')
        else:
            bc_fq_idx, paired_fq_idx = 0, 1

        all_reads_qcs = []
        for fpath_tup in paired_fpaths:
            bc_fq1_fpath = fpath_tup[bc_fq_idx]
            bc_fq2_fpath = fpath_tup[paired_fq_idx]

            sans_bc_fq1_fpath = os.path.join(
                arguments.output_dir,
                f'sans_bc_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.fq'
            )
            sans_bc_fq2_fpath = os.path.join(
                arguments.output_dir,
                f'sans_bc_{misc.file_prefix_from_fpath(bc_fq2_fpath)}.fq'
            )

            tags1_fpath = os.path.join(
                arguments.output_dir,
                f'rec_names_and_tags1_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.txt'
            )
            tags2_fpath = os.path.join(
                arguments.output_dir,
                f'rec_names_and_tags2_{misc.file_prefix_from_fpath(bc_fq2_fpath)}.txt'
            )

            namepairidx_fpath = os.path.join(
                arguments.output_dir,
                f'namepairidx_{misc.file_prefix_from_fpath(bc_fq1_fpath)}.pkl'
            )

            namepairidx = misc.get_namepair_index(bc_fq1_fpath, bc_fq2_fpath)
            namepairidxs.append(namepairidx)
            with open(namepairidx_fpath, "wb") as pkl:
                pickle.dump(namepairidx, pkl)

            log.info('Processing files:')
            if bcs_on_both_reads:
                log.info(f'  barcode fastq 1: {bc_fq1_fpath}')
                log.info(f'  barcode fastq 2: {bc_fq2_fpath}')
                log.info(f'  sans-bc fastq 1: {sans_bc_fq1_fpath}')
                log.info(f'  sans-bc fastq 2: {sans_bc_fq2_fpath}')
                log.info(f'  tags file 1:     {tags1_fpath}')
                log.info(f'  tags file 2:     {tags2_fpath}')
            else:
                log.info(f'  barcode fastq: {bc_fq1_fpath}')
                log.info(f'  paired fastq:  {bc_fq2_fpath}')
                log.info(f'  sans-bc fastq: {sans_bc_fq1_fpath}')
                log.info(f'  sans-bc fastq: {sans_bc_fq2_fpath}')
                log.info(f'  tags file:     {tags1_fpath}')

            if any(os.path.exists(fpath) for fpath in [sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath]):
                raise FileExistsError('Partial results found: Remove partial results and restart.')

            file_read_qcs = process_fastqs_func(
                arguments,
                bc_fq1_fpath,
                bc_fq2_fpath,
                sans_bc_fq1_fpath,
                sans_bc_fq2_fpath,
                tags1_fpath,
                tags2_fpath,
                bcs_on_both_reads
            )

            all_reads_qcs.extend(file_read_qcs)

            if bc_fq_idx == 0:
                paired_align_fqs_and_tags_fpaths.append((sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath))
            else:
                paired_align_fqs_and_tags_fpaths.append((sans_bc_fq2_fpath, sans_bc_fq1_fpath, tags2_fpath, tags1_fpath))

        qc_metrics.generate_qc_metrics(arguments, all_reads_qcs)
        qc_metrics.generate_qc_metrics_reads(arguments, all_reads_qcs)
        return paired_align_fqs_and_tags_fpaths, namepairidxs

def process_fastqs(arguments):
    """
    Output single file with parsed bcs from bc_fastq in read names and seqs from paired_fastq.
    """
    star_w_bc_fpath = os.path.join(arguments.output_dir, 'with_bc.bam')
    star_w_bc_sorted_fpath = os.path.join(arguments.output_dir, 'with_bc.sorted.bam')
    star_w_bc_umi_fpath = os.path.join(arguments.output_dir, 'with_bc_umi.bam')
    star_w_bc_umi_sorted_fpath = os.path.join(arguments.output_dir, 'with_bc_umi.sorted.bam')
    if os.path.exists(star_w_bc_umi_sorted_fpath):
        raise FileExistsError('Sorted bam with UMIs found. Run possibly previously completed.')
    elif os.path.exists(star_w_bc_umi_fpath):
        raise FileExistsError('Partial results found: Corrected UMI file. Remove partial results and restart.')
    elif os.path.exists(star_w_bc_sorted_fpath):
        raise FileExistsError('Partial results found: Sorted STAR results. Remove partial results and restart.')
    elif os.path.exists(star_w_bc_fpath):
        raise FileExistsError('Partial results found: STAR results. Remove partial results and restart.')

    if not os.path.exists(arguments.output_dir):
        os.makedirs(arguments.output_dir)
    # find barcodes, write to file, remove from reads before mapping
    # map reads
    # Add barcodes etc to bam file
    # Sort, index
    log.info('Writing output to:')
    log.info(f'  {star_w_bc_fpath}')
    do_single_file_pairing = False

    if arguments.single_end_reads:
        if not arguments.star_output_path: # This branch was deprecated, but the code still remains for future use. 
            raise ValueError(
                "Internal Processing with STAR not available."
                "Please provide a path to the unsorted aligned BAM/SAM file with --STAR-output-dir"
            ) 
            single_align_fqs_and_tags_fpaths = preprocess_fastqs(arguments)

            log.info('Running STAR alignment...')
            star_out_dirs = set()
            star_bam_and_tags_fpaths = []

            for R1_fpath, tags1_fpath in single_align_fqs_and_tags_fpaths:
                if R1_fpath:
                    log.info(f'  {R1_fpath}')
                if R1_fpath != single_align_fqs_and_tags_fpaths[-1][0]:
                    log.info('  -')

                star_out_dir, star_out_fpath = run_STAR_single_end(arguments, R1_fpath)
                star_out_dirs.add(star_out_dir)
                star_bam_and_tags_fpaths.append((star_out_fpath, tags1_fpath))
        else:
            log.info(f'Using STAR results from {arguments.star_output_path}')
            bam_fpaths = sorted(glob(os.path.join(arguments.star_output_path, "*Aligned.out.bam")))
            tag_fpaths = sorted(glob(os.path.join(arguments.fastq_dir, "rec_names_and_tags1_*.txt")))

            if len(bam_fpaths) == 0:
                raise FileNotFoundError(
                    f"No STAR BAM files found in STAR output directory: {arguments.star_output_path}"
                )
            if len(tag_fpaths) == 0:
                raise FileNotFoundError(
                    f"No tag files found in fastq/output directory: {arguments.fastq_dir}"
                )
            if len(bam_fpaths) != len(tag_fpaths):
                raise ValueError(
                    f"Unequal number of BAM files and rec_names_and_tags files found. \
                    {len(bam_fpaths)} BAM files found: {bam_fpaths}. \
                    {len(tag_fpaths)} rec_names_and_tags files found: {tag_fpaths}. "
                    )
            star_bam_and_tags_fpaths = list(zip(bam_fpaths, tag_fpaths))

    else:        
        if not arguments.star_output_path: # This branch was deprecated, but the code still remains for future use. 
            raise ValueError(
                "Internal Processing with STAR not available."
                "Please provide a path to the unsorted aligned BAM/SAM file with --STAR-output-dir"
            )
            paired_align_fqs_and_tags_fpaths, namepairidxs = preprocess_fastqs(arguments)

            log.info(f'Running STAR alignment...')
            star_out_dirs = set()
            star_bam_and_tags_fpaths = []
            for R1_fpath, R2_fpath, tags1_fpath, tags2_fpath in paired_align_fqs_and_tags_fpaths:
                r1_empty = misc.check_if_sans_bc_is_empty(R1_fpath)
                r2_empty = misc.check_if_sans_bc_is_empty(R2_fpath)

                if r1_empty and r2_empty:
                    raise ValueError(
                        "After preprocessing, both paired-end sans_bc FASTQ files contain no alignable sequence. "
                        "This usually means all sequence content was removed by barcode trimming / keep_nonbarcode settings. "
                        "No reads remain that STAR can align."
                    )

                if r1_empty:
                    warnings.warn(
                        f"Preprocessed read 1 FASTQ has no alignable sequence left and will be omitted from STAR alignment: {R1_fpath}",
                        category=UserWarning
                    )
                    R1_fpath = None

                if r2_empty:
                    warnings.warn(
                        f"Preprocessed read 2 FASTQ has no alignable sequence left and will be omitted from STAR alignment: {R2_fpath}",
                        category=UserWarning
                    )
                    R2_fpath = None

                if R1_fpath:
                    log.info(f'  {R1_fpath}')
                log.info(f'  {R2_fpath}')
                if R1_fpath != paired_align_fqs_and_tags_fpaths[-1][0]:
                    log.info('  -')
                star_out_dir, star_out_fpath = run_STAR(arguments, R1_fpath, R2_fpath)
                star_out_dirs.add(star_out_dir)
                star_bam_and_tags_fpaths.append((star_out_fpath, tags1_fpath, tags2_fpath))
        else:
            log.info(f'Using STAR results from {arguments.star_output_path}')
            bam_fpaths = sorted(glob(os.path.join(arguments.star_output_path, "*Aligned.out.bam")))
            tag1_fpaths = sorted(glob(os.path.join(arguments.fastq_dir, "rec_names_and_tags1_*.txt")))
            tag2_fpaths = sorted(glob(os.path.join(arguments.fastq_dir, "rec_names_and_tags2_*.txt")))
            namepairidx_fpaths = sorted(glob(os.path.join(arguments.fastq_dir, "namepairidx_*.pkl")))

            if len(tag1_fpaths) == 0 and len(tag2_fpaths) == 0:
                raise FileNotFoundError(
                    f"No tag files found in fastq/output directory: {arguments.fastq_dir}"
                )
            elif len(tag2_fpaths) == 0 and len(tag1_fpaths) > 0:
                log.info(
                    "Found rec_names_and_tags1_*.txt files but no rec_names_and_tags2_*.txt files. "
                    "This is normal for barcodes only on the forward read. "
                    "Proceeding with only the rec_names_and_tags1_*.txt files."
                )
                tag_fpaths = tag1_fpaths
            elif len(tag1_fpaths) == 0 and len(tag2_fpaths) > 0:
                log.info(
                    "Found rec_names_and_tags2_*.txt files but no rec_names_and_tags1_*.txt files. "
                    "This is normal for barcodes only on the reverse read. "
                    "Proceeding with only the rec_names_and_tags2_*.txt files."
                )
                tag_fpaths = tag2_fpaths
            elif len(tag1_fpaths) > 0 and len(tag2_fpaths) > 0 and len(tag1_fpaths) != len(tag2_fpaths):
                raise ValueError(
                    f"Unequal number of rec_names_and_tags1_*.txt files and rec_names_and_tags2_*.txt files found. "
                    f"{len(tag1_fpaths)} rec_names_and_tags1_*.txt files found: {tag1_fpaths}. "
                    f"{len(tag2_fpaths)} rec_names_and_tags2_*.txt files found: {tag2_fpaths}."
                )
            elif  len(tag1_fpaths) > 0 and len(tag2_fpaths) > 0 and len(tag1_fpaths) == len(tag2_fpaths):
                log.info(
                    f"Found equal number of rec_names_and_tags1_*.txt files and rec_names_and_tags2_*.txt files. "
                    f"Proceeding with paired tag files."
                )
                tag_fpaths = [fpath for pair in zip(tag1_fpaths, tag2_fpaths) for fpath in pair]
            else:
                raise ValueError(
                    f"Unexpected combination of rec_names_and_tags1_*.txt files and rec_names_and_tags2_*.txt files found. "
                    f"This should never trigger if the above cases are correctly specified, but was added as a sanity check. "
                    f"{len(tag1_fpaths)} rec_names_and_tags1_*.txt files found: {tag1_fpaths}. "
                    f"{len(tag2_fpaths)} rec_names_and_tags2_*.txt files found: {tag2_fpaths}."
                )
            
            if len(bam_fpaths) == 0:
                raise FileNotFoundError(
                    f"No STAR BAM files found in STAR output directory: {arguments.star_output_path}"
                )
            if len(namepairidx_fpaths) == 0:
                raise FileNotFoundError(
                    f"No namepairidx .pkl files found in fastq/output directory: {arguments.fastq_dir}"
                )
            if len(namepairidx_fpaths) != len(bam_fpaths):
                raise ValueError(
                    f"Unequal number of BAM files and namepairidx files found. "
                    f"{len(bam_fpaths)} BAM files found: {bam_fpaths}. "
                    f"{len(namepairidx_fpaths)} namepairidx files found: {namepairidx_fpaths}."
                )
            if len(bam_fpaths) == len(tag_fpaths):
                log.info(
                    "Found equal number of BAM and tag files; treating outsourced STAR input like single-file pairing."
                )
                do_single_file_pairing = True
                star_bam_and_tags_fpaths = list(zip(bam_fpaths, tag_fpaths))
            elif len(tag_fpaths) == 2 * len(bam_fpaths):
                log.info(
                    f"Found {len(bam_fpaths)} BAM file(s) and {len(tag_fpaths)} tag file(s); "
                    "treating tags as paired entries per BAM."
                )
                star_bam_and_tags_fpaths = list(zip(bam_fpaths, tag1_fpaths, tag2_fpaths))
            else:
                raise ValueError(
                    "Mismatch between STAR BAM files and tag files in outsourced paired-end mode. "
                    f"Found {len(bam_fpaths)} BAM file(s): {bam_fpaths} "
                    f"and {len(tag_fpaths)} tag file(s): {tag_fpaths}. "
                    "Expected either:\n"
                    "  - exactly 1 tag per BAM file (barcodes only on one read direction), or\n"
                    "  - exactly 2 tag files per BAM (barcodes on both read directions)."
                )                            
            namepairidxs = [pickle.load(open(fpath, "rb")) for fpath in namepairidx_fpaths]


    if arguments.single_end_reads or do_single_file_pairing:
        log.info('Adding barcode tags to mapped reads... ')
        template_bam_fpath = star_bam_and_tags_fpaths[0][0]
        with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as bam_out:
            for star_out_fpath, tags1_fpath in star_bam_and_tags_fpaths:
                if arguments.threads == 1:
                    readsiter = iter(serial_add_tags_to_reads_single_end(tags1_fpath, star_out_fpath))
                else:
                    readsiter = iter(parallel_add_tags_to_reads_single_end(tags1_fpath, star_out_fpath, arguments.threads))

                for read in readsiter:
                    bam_out.write(read)
    else:
        log.info('Adding barcode tags to mapped reads...')
        template_bam_fpath = star_bam_and_tags_fpaths[0][0]
        with pysam.AlignmentFile(star_w_bc_fpath, 'wb', template=pysam.AlignmentFile(template_bam_fpath)) as bam_out:
            for (star_out_fpath, tags1_fpath, tags2_fpath), namepairidx in zip(star_bam_and_tags_fpaths, namepairidxs):
                if arguments.threads == 1:
                    readsiter = iter(serial_add_tags_to_reads(tags1_fpath, tags2_fpath, star_out_fpath))
                else:
                    readsiter = iter(parallel_add_tags_to_reads(tags1_fpath, tags2_fpath, star_out_fpath, namepairidx, arguments.threads))
                for read in readsiter:
                    bam_out.write(read)



    log.info('Sorting bam...')
    pysam.sort('-@', str(arguments.threads), '-o', star_w_bc_sorted_fpath, star_w_bc_fpath)
    log.info('Indexing bam...')
    pysam.index(star_w_bc_sorted_fpath)

    log.info('Correcting UMIs...')
    correct_UMIs(star_w_bc_sorted_fpath, star_w_bc_umi_fpath, arguments.threads)

    log.info('Sorting bam...')
    pysam.sort('-@', str(arguments.threads), '-o', star_w_bc_umi_sorted_fpath, star_w_bc_umi_fpath)
    log.info('Indexing bam...')
    pysam.index(star_w_bc_umi_sorted_fpath)
    count_matrix(arguments, star_w_bc_umi_sorted_fpath)

    handle_intermediary_files(arguments, star_w_bc_umi_sorted_fpath)
    qc_metrics.make_stacked_barplot_blocks(arguments)
    qc_metrics.make_stacked_barplots_bcs(arguments)
    qc_metrics.make_stacked_barplot_reads(arguments)

    log.info('Done')


def handle_intermediary_files(arguments, star_w_bc_umi_sorted_fpath):
    """
    Delete or move all files/subdirectories in output_dir except the final outputs.

    Final outputs currently preserved:
    - raw_reads_bc_matrix
    - raw_umis_bc_matrix
    - with_bc_umi.sorted.bam
    - with_bc_umi.sorted.bam.bai
    - QC_metrics/   (contains QC_metrics.tsv and QC_metrics_stacked_barplot.png)

    If the final outputs are updated, this function must be updated as well!
    """

    output_dir = arguments.output_dir
    intermediary_dir = os.path.join(output_dir, "intermediary_files")

    qc_dir = os.path.join(output_dir, "QC_metrics")

    # Define paths to keep
    paths_to_keep = {
        star_w_bc_umi_sorted_fpath,
        star_w_bc_umi_sorted_fpath + ".bai",
        os.path.join(output_dir, "raw_reads_bc_matrix"),
        os.path.join(output_dir, "raw_umis_bc_matrix"),
        qc_dir,
    }

    # List all top-level items in output_dir
    all_items = [os.path.join(output_dir, item) for item in os.listdir(output_dir)]

    if not arguments.keep_intermediary_files:
        log.info("Deleting intermediary files...")
        for item in all_items:
            # Skip if item is in the keep list
            if item in paths_to_keep:
                continue

            if os.path.isfile(item):
                os.remove(item)
            elif os.path.isdir(item):
                shutil.rmtree(item)

    else:
        log.info("Preserving intermediary files in: %s", intermediary_dir)
        os.makedirs(intermediary_dir, exist_ok=True)

        for item in all_items:
            # Skip if item is in the keep list
            if item in paths_to_keep:
                continue
            # Move item to intermediary directory
            shutil.move(item, os.path.join(intermediary_dir, os.path.basename(item)))


def run_STAR(arguments, R1_fpath, R2_fpath):
    """
    Run STAR aligner for gDNA files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    R2_bname = misc.file_prefix_from_fpath(R2_fpath)
    out_prefix = os.path.join(star_out_dir, f'{R2_bname}_')
    cmd_star = [
        'STAR',
        f'--runThreadN {arguments.threads}',
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {R1_fpath if R1_fpath else ""} {R2_fpath if R2_fpath else ""}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1',
        '--outSAMtype BAM Unsorted',
        '--outSAMattributes NH HI AS nM GX GN',
    ]
    if R2_fpath.endswith('.gz'):
        if R1_fpath and not R1_fpath.endswith('.gz'):
            raise ValueError('Paired read files must be both zipped or both unzipped')
        cmd_star.append('--readFilesCommand zcat')
    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    if os.path.exists(star_out_fpath):
        log.info("STAR results found. Skipping alignment")
    else:
        subprocess.run(cmd_star, check=True)
    return star_out_dir, star_out_fpath

def run_STAR_single_end(arguments, R1_fpath):
    """
    Run STAR aligner for single-end files.

    Returns STAR output directory and bam path.
    """
    star_out_dir = os.path.join(arguments.output_dir, 'STAR_files')
    R1_bname = misc.file_prefix_from_fpath(R1_fpath)
    out_prefix = os.path.join(star_out_dir, f'{R1_bname}_')

    cmd_star = [
        'STAR',
        f'--runThreadN {arguments.threads}',
        f'--genomeDir {arguments.star_ref_dir}',
        f'--readFilesIn {R1_fpath}',
        f'--outFileNamePrefix {out_prefix}',
        '--outFilterMultimapNmax 1',
        '--outSAMtype BAM Unsorted',
        '--outSAMattributes NH HI AS nM GX GN',
    ]

    if R1_fpath.endswith('.gz'):
        cmd_star.append('--readFilesCommand zcat')

    star_out_fpath = f'{out_prefix}Aligned.out.bam'
    if os.path.exists(star_out_fpath):
        log.info("STAR results found. Skipping alignment")
    else:
        subprocess.run(cmd_star, check=True)

    return star_out_dir, star_out_fpath

def process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders):
    """
    Find barcodes etc in bc_rec 
    """
    
    scores_pieces_end_bounds = [al.find_norm_score_pieces_and_boundaries(bc_rec.seq, return_seq=True) for al in aligners]
    scores_pieces_end_pos = [(s, p, bounds[-1], sq) for s, p, bounds, sq in scores_pieces_end_bounds]       

    raw_score, raw_pieces, raw_end_pos, raw_seq = max(scores_pieces_end_pos)
    raw_bounds = next(
        bounds for s, p, bounds, sq in scores_pieces_end_bounds
        if s == raw_score and p == raw_pieces
    )

    raw_bcs = [raw_piece.upper().replace('N', 'A') for raw_piece, block in zip(raw_pieces, blocks) if block["blocktype"] == "barcodeList"]
    decode_results = [decoder.decode_with_status(raw_bc) for raw_bc, decoder in zip(raw_bcs, decoders)]
    bcs, statuses, overlapping_bcs = map(list, zip(*decode_results))

    best_aligner = next(al for al, (s, p, e, sq) in zip(aligners, scores_pieces_end_pos) if s == raw_score)

    commonseqs = [best_aligner.prefixes[i] for i, block in enumerate(blocks) if block["blocktype"] == "constantRegion"]

    # QC Metrics
    bc_qc = {
        "read_name": str(bc_rec.id),
        "raw_bcs": raw_bcs,
        "decoded_bcs": bcs,
        "statuses": statuses,
        "conflict_bcs": overlapping_bcs,
    }

    bc_idx = commonseq_idx = 0
    corrected_pieces = []

    for i, block in enumerate(blocks):
        if block["blocktype"] == "barcodeList":
            if bcs[bc_idx] is None:
                return raw_score, None, None, bc_qc
            corrected_pieces.append(bcs[bc_idx])
            bc_idx += 1

        elif block["blocktype"] == "constantRegion":
            if commonseqs[commonseq_idx] is None:
                return raw_score, None, None, bc_qc

            observed_const = raw_pieces[i].upper().replace("N", "A")
            expected_const = commonseqs[commonseq_idx]

            distfun = misc.DistanceThresh("levenshtein", block["maxerrors"])
            dist = distfun(observed_const, expected_const)
            if dist is False:
                return raw_score, None, None, bc_qc

            corrected_pieces.append(expected_const)
            commonseq_idx += 1

        else:
            corrected_pieces.append('N' * block["length"])

    new_aligner = CustomBCAligner(*corrected_pieces)
    new_score, new_pieces = new_aligner.find_norm_score_and_pieces(raw_seq)

    bam_bcs = []
    bam_raw_bcs = []
    bam_umis = []
    bam_commonseqs = []

    bc_idx = commonseq_idx = 0
    for i, block in enumerate(blocks):
        if block["blockfunction"] != "discard":
            if block["blocktype"] == "barcodeList":
                bam_bcs.append(bcs[bc_idx])
                bam_raw_bcs.append(raw_pieces[i])
                bc_idx += 1
            elif block["blocktype"] == "constantRegion":
                bam_commonseqs.append(commonseqs[commonseq_idx])
                commonseq_idx += 1
            else:
                if block["blockfunction"] == "UMI":
                    bam_umis.append(new_pieces[i])
                else:
                    bam_bcs.append(new_pieces[i])
                    bam_raw_bcs.append(new_pieces[i])


    tags = []
    # Cell barcode
    tags.append(('CB', '.'.join(bam_bcs)))
    tags.append(('CR', '.'.join(bam_raw_bcs)))
    # Filler sequences
    tags.append(('FL', sum(len(seq) for seq in bam_commonseqs)))
    # And raw UMI
    tags.append(("UR", ".".join(bam_umis)))

    starts = [0] + list(raw_bounds[:-1])
    quals = bc_rec.letter_annotations["phred_quality"]

    keep_seq_parts = []
    keep_qual_parts = []

    # Keep selected defined blocks
    for block, piece, start, end in zip(blocks, raw_pieces, starts, raw_bounds):
        if block.get("keepblock", False):
            keep_seq_parts.append(piece)
            keep_qual_parts.extend(quals[start:end])

    # Keep undefined / trailing sequence if requested
    if keep_nonbarcode:
        trailing_rec = bc_rec[raw_end_pos:]
        keep_seq_parts.append(str(trailing_rec.seq))
        keep_qual_parts.extend(trailing_rec.letter_annotations["phred_quality"])

    new_seq = "".join(keep_seq_parts)

    new_rec = SeqRecord(
        Seq(new_seq),
        id=bc_rec.id,
        description=bc_rec.description,
    )
    new_rec.letter_annotations["phred_quality"] = keep_qual_parts

    return raw_score, new_rec, tags, bc_qc

def umi_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    return ref, get_umi_maps_from_bam_file(input_bam_fpath, chrm=ref)

def correct_UMIs(input_bam_fpath, out_bam_fpath, threads=1):
    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references
    reference_names_with_input_bam = [(ref, input_bam_fpath) for ref in reference_names]

    with pysam.AlignmentFile(out_bam_fpath, 'wb', template=pysam.AlignmentFile(input_bam_fpath)) as bam_out, \
            Pool(threads) as pool:
        for i, (ref, umi_map_given_bc_then_feature) in enumerate(pool.imap_unordered(
                umi_parallel_wrapper,
                reference_names_with_input_bam)):
            log.info(f'  {ref}')
            for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
                for gx_gn_tup in misc.gx_gn_tups_from_read(read):
                    corrected_umi = umi_map_given_bc_then_feature[read.get_tag('CB')][gx_gn_tup][read.get_tag('UR')]
                    read.set_tag('UB', corrected_umi)
                    bam_out.write(read)
                    break # use the first one. Ideally same across all


def build_tags_iter(tags_fpath):
    for line in open(tags_fpath):
        words = line.strip().split('\t')
        read_name = words[0]
        tags = [parse_tag_str(word) for word in words[1:]]
        yield read_name, tags


def parse_tag_str(tag_str):
    tag, tag_type, val = tag_str.split(':')
    if tag_type == 'i':
        # only worry about strings and ints
        val = int(val)
    return tag, val


def merge_tags(tags1, tags2):
    """
    Merge tags from read 1 and read 2.

    String-like duplicated tags are concatenated with a dot.
    Integer duplicated tags, such as FL, are summed.
    Ask John if summing FL is the right thing to do,
    but it seems reasonable since FL is just the total length of constant regions,
    so summing should give the correct total length if both reads have constant regions.
    """
    merged = dict(tags1)

    for name, val in dict(tags2).items():
        if name not in merged:
            merged[name] = val
        elif isinstance(merged[name], int) and isinstance(val, int):
            merged[name] += val
        else:
            merged[name] = ".".join([str(merged[name]), str(val)])

    return merged.items()

def serial_add_tags_to_reads(tags1_fpath, tags2_fpath, bam_fpath):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    tags1_iter = build_tags_iter(tags1_fpath)
    tags2_iter = build_tags_iter(tags2_fpath) if os.path.isfile(tags2_fpath) else itertools.repeat(("", ()))
    read_name = 'Lorem ipsum'
    for read in pysam.AlignmentFile(bam_fpath).fetch(until_eof=True):
        while not misc.names_pair(read_name, str(read.query_name)):
            read_name, tags1 = next(tags1_iter)
            _, tags2 = next(tags2_iter)
        for name, val in merge_tags(tags1, tags2):
            read.set_tag(name, val)
        yield read

def serial_add_tags_to_reads_single_end(tags1_fpath, bam_fpath):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    tags1_iter = build_tags_iter(tags1_fpath)
    read_name = None
    tags1 = None

    for read in pysam.AlignmentFile(bam_fpath).fetch(until_eof=True):
        while read_name != str(read.query_name):
            read_name, tags1 = next(tags1_iter)
        for name, val in tags1:
            read.set_tag(name, val)
        yield read

def parallel_add_tags_to_reads(tags1_fpath, tags2_fpath, bam_fpath, namepairidx, threads):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    sortedbampath = bam_fpath + "_readname_sorted.bam"
    bamidx = misc.sort_and_index_readname_bam(bam_fpath, sortedbampath, namepairidx, threads)
    tags1_iter = build_tags_iter(tags1_fpath)
    tags2_iter = build_tags_iter(tags2_fpath) if os.path.isfile(tags2_fpath) else itertools.repeat(("", ()))
    with pysam.AlignmentFile(sortedbampath) as bam:
        for (read_name, tags1), (_, tags2) in zip(tags1_iter, tags2_iter):
            for read in misc.get_bam_read_by_name(read_name, bam, bamidx):
                for name, val in merge_tags(tags1, tags2):
                    read.set_tag(name, val)
                yield read
    os.remove(sortedbampath)

def parallel_add_tags_to_reads_single_end(tags1_fpath, bam_fpath, threads):
    """
    Iterates bam reads and matching tags records and adds tags to reads.
    """
    sortedbampath = bam_fpath + "_readname_sorted.bam"
    bamidx = misc.sort_and_index_readname_bam_single_end(bam_fpath, sortedbampath, threads)
    tags1_iter = build_tags_iter(tags1_fpath)

    with pysam.AlignmentFile(sortedbampath) as bam:
        for read_name, tags1 in tags1_iter:
            for read in misc.get_bam_read_by_name_single_end(read_name, bam, bamidx):
                for name, val in tags1:
                    read.set_tag(name, val)
                yield read

    os.remove(sortedbampath)

def tag_type_from_val(val):
    if isinstance(val, str):
        return 'Z'
    elif isinstance(val, int):
        return 'i'
    else:
        # no other types implemented
        raise ValueError('Unexpected tag type')


def output_rec_name_and_tags(rec, tags, out_fh):
    out_fh.write('\t'.join([f'{str(rec.id)}'] + [f'{tag}:{tag_type_from_val(val)}:{val}' for tag, val in tags]) + '\n')


def get_per_fastq_qc_fastq_paths(arguments, fq_fpath):
    """
    Return ambiguous/no_match FASTQ paths for one input FASTQ file.
    """
    qc_dir = qc_metrics.get_qc_paths(arguments)["qc_dir"]
    prefix = misc.file_prefix_from_fpath(fq_fpath)
    return {
        "ambiguous": os.path.join(qc_dir, f"{prefix}_ambiguous_reads.fq"),
        "no_match": os.path.join(qc_dir, f"{prefix}_no_match_reads.fq"),
    }


def clone_rec_with_new_description(rec, description):
    """
    Make a FASTQ-safe copy of a SeqRecord with the same seq/quals and a new header.
    """
    new_rec = SeqRecord(
        rec.seq,
        id=rec.id,
        description=description,
    )
    if "phred_quality" in rec.letter_annotations:
        new_rec.letter_annotations["phred_quality"] = list(
            rec.letter_annotations["phred_quality"]
        )
    return new_rec


def build_qc_failure_records(rec, read_qc, meta_list):
    """
    From one raw read + its QC info, build:
      - ambiguous FASTQ record (or None)
      - no_match FASTQ record (or None)

    A read can appear in both outputs if it has both ambiguous and no_match blocks.
    """
    ambiguous_details = []
    no_match_blocks = []

    for meta, status, conflicts in zip(
        meta_list,
        read_qc["statuses"],
        read_qc["conflict_bcs"],
    ):
        blockname = meta["blockname"]

        if status == "ambiguous":
            candidate_bcs, conflict_meta = qc_metrics.parse_conflicts(conflicts)

            if candidate_bcs:
                detail = f"{blockname}:{','.join(candidate_bcs)}"
            elif conflict_meta:
                detail = f"{blockname}:{','.join(map(str, conflict_meta))}"
            else:
                detail = f"{blockname}:unidentified"

            ambiguous_details.append(detail)

        elif status == "no_match":
            no_match_blocks.append(blockname)

    ambiguous_rec = None
    no_match_rec = None

    if ambiguous_details:
        ambiguous_desc = (
            f"{rec.id} QC_status=ambiguous "
            f"conflicting_barcodes={'|'.join(ambiguous_details)}"
        )
        ambiguous_rec = clone_rec_with_new_description(rec, ambiguous_desc)

    if no_match_blocks:
        no_match_desc = (
            f"{rec.id} QC_status=no_match "
            f"failed_blocks={'|'.join(no_match_blocks)}"
        )
        no_match_rec = clone_rec_with_new_description(rec, no_match_desc)

    return ambiguous_rec, no_match_rec


def write_qc_failure_records(rec, read_qc, meta_list, ambiguous_fh, no_match_fh):
    """
    Write raw reads with QC annotations into the appropriate FASTQ files.
    """
    ambiguous_rec, no_match_rec = build_qc_failure_records(rec, read_qc, meta_list)

    if ambiguous_rec is not None and ambiguous_fh is not None:
        SeqIO.write(ambiguous_rec, ambiguous_fh, "fastq")

    if no_match_rec is not None and no_match_fh is not None:
        SeqIO.write(no_match_rec, no_match_fh, "fastq")


def serial_process_fastqs(arguments, fq1_fpath, fq2_fpath, sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath, bcs_on_both_reads):
    blocks = [arguments.config["barcode_struct_r1"]["blocks"]]
    keep_nonbarcodes = [arguments.config["barcode_struct_r1"]["keep_nonbarcode"]]
    if bcs_on_both_reads:
        blocks.append(arguments.config["barcode_struct_r2"]["blocks"])
        keep_nonbarcodes.append(arguments.config["barcode_struct_r2"]["keep_nonbarcode"])

    log.info('Building aligners and barcode decoders')
    aligners = tuple(misc.build_bc_aligners(bblocks, arguments.config["unknown_read_orientation"]) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    thresh = [np.inf, -np.inf]
    for i, (fq_fpath, bblocks, bkeep_nonbarcodes, baligners, bdecoders) in enumerate(zip((fq1_fpath, fq2_fpath), blocks, keep_nonbarcodes, aligners, decoders, strict=False)):
        first_scores_recs_tags = []
        for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq_fpath), 'fastq')):
            first_scores_recs_tags.append(process_bc_rec(arguments, bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders)[:-1]) # Append the tags (excluding QC in the last tupel)
            if j >= n_first_seqs:
                break

        scores = [score for score, rec, tags in first_scores_recs_tags]
        thresh[i] = np.average(scores) - 2 * np.std(scores)
    
    r1_meta = qc_metrics.get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = qc_metrics.get_config_metadata(arguments, "barcode_struct_r2") if bcs_on_both_reads else []

    r1_qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq1_fpath)
    r2_qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq2_fpath) if bcs_on_both_reads else None

    log.info(f'Score threshold: {thresh[0]:.2f}{f"{thresh[1]:.2f}" if thresh[1] != -np.inf else ""}')

    total_out = 0
    all_reads_qcs = []

    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            (open(sans_bc_fq2_fpath, 'w') if sans_bc_fq2_fpath else nullcontext(None)) as bc_fq2_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            (open(tags2_fpath, 'w') if bcs_on_both_reads else nullcontext(None)) as tag2_fh, \
            open(r1_qc_fastq_paths["ambiguous"], "w") as r1_ambiguous_fh, \
            open(r1_qc_fastq_paths["no_match"], "w") as r1_no_match_fh, \
            (open(r2_qc_fastq_paths["ambiguous"], "w") if bcs_on_both_reads else nullcontext(None)) as r2_ambiguous_fh, \
            (open(r2_qc_fastq_paths["no_match"], "w") if bcs_on_both_reads else nullcontext(None)) as r2_no_match_fh:

        log.info('Continuing...')
        for i, (bc_rec1, bc_rec2) in enumerate(zip(
            SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq'),
            SeqIO.parse(misc.gzip_friendly_open(fq2_fpath), 'fastq'))):

            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d} processed,  {total_out:,d} output')

            scores = [0, 0]
            sans_bc_rec = [None, None]
            tags = [None, None]
            read_qcs = [None, None]

            for j, (bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders) in enumerate(
                zip(blocks, keep_nonbarcodes, (bc_rec1, bc_rec2), aligners, decoders, strict=False)
            ):
                read_first_scores_recs_tags_plus_qc = process_bc_rec(arguments, bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders)
                read_qc = read_first_scores_recs_tags_plus_qc[-1] # get our QC (last item of the tuple)
                read_qcs[j] = read_qc
                scores[j], sans_bc_rec[j], tags[j] = read_first_scores_recs_tags_plus_qc[:-1] # get the first three tuple entries

            # write raw failure reads
            write_qc_failure_records(
                bc_rec1,
                read_qcs[0],
                r1_meta,
                r1_ambiguous_fh,
                r1_no_match_fh,
            )

            if bcs_on_both_reads and read_qcs[1] is not None:
                write_qc_failure_records(
                    bc_rec2,
                    read_qcs[1],
                    r2_meta,
                    r2_ambiguous_fh,
                    r2_no_match_fh,
                )

            if all(score >= cthresh for score, cthresh in zip(scores, thresh)) and all(sans_bc_rec):
                total_out += 1
                for bc_rec, bc_fq_fh, csans_bc_rec, ctags, tag_fh in zip(
                    (bc_rec1, bc_rec2),
                    (bc_fq1_fh, bc_fq2_fh),
                    sans_bc_rec,
                    tags,
                    (tag1_fh, tag2_fh)
                ):
                    if bc_fq_fh:
                        if csans_bc_rec and csans_bc_rec is not True:
                            SeqIO.write(csans_bc_rec, bc_fq_fh, 'fastq')
                        else:
                            SeqIO.write(bc_rec, bc_fq_fh, 'fastq')
                    if ctags:
                        output_rec_name_and_tags(bc_rec, ctags, tag_fh)

            all_reads_qcs.append(tuple(read_qcs[:len(blocks)]))

    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')
    return all_reads_qcs

def serial_process_fastqs_single_end(arguments, fq1_fpath, sans_bc_fq1_fpath, tags1_fpath):
    blocks = arguments.config["barcode_struct_r1"]["blocks"]
    keep_nonbarcode = arguments.config["barcode_struct_r1"]["keep_nonbarcode"]

    log.info('Building aligners and barcode decoders')
    aligners = misc.build_bc_aligners(blocks, unknown_read_orientation=False)
    decoders = misc.build_bc_decoders(blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    first_scores_recs_tags = []
    for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq')):
        first_scores_recs_tags.append(process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders)[:-1]) # get the first three tuple entries 
        if j >= n_first_seqs:
            break

    scores = [score for score, rec, tags in first_scores_recs_tags]
    thresh = np.average(scores) - 2 * np.std(scores)

    log.info(f'Score threshold: {thresh:.2f}')

    total_out = 0
    n_processed = 0
    all_reads_qcs = []

    r1_meta = qc_metrics.get_config_metadata(arguments, "barcode_struct_r1")
    qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq1_fpath)

    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            open(qc_fastq_paths["ambiguous"], "w") as ambiguous_fh, \
            open(qc_fastq_paths["no_match"], "w") as no_match_fh:
        log.info('Continuing...')
        for i, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq')):
            n_processed = i + 1

            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d} processed,  {total_out:,d} output')

            read_first_scores_recs_tags_plus_qc = process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders) # get our tags for each read and trim the read
            read_qc = read_first_scores_recs_tags_plus_qc[-1] # get our QC (last item of the tuple)
            all_reads_qcs.append(read_qc)
            score, sans_bc_rec, tags = read_first_scores_recs_tags_plus_qc[:-1] # get the first three tuple entries

            # write raw failure reads
            write_qc_failure_records(
                bc_rec,
                read_qc,
                r1_meta,
                ambiguous_fh,
                no_match_fh,
            )

            if score >= thresh and sans_bc_rec:
                total_out += 1

                if bc_fq1_fh:
                    SeqIO.write(sans_bc_rec, bc_fq1_fh, 'fastq')

                if tags:
                    output_rec_name_and_tags(bc_rec, tags, tag1_fh)

    log.info(f'{n_processed:,d} barcode records processed')
    log.info(f'{total_out:,d} records output')
    return all_reads_qcs

def worker_build_aligners(aarguments, ablocks, akeep_nonbarcodes, unknown_read_orientation, athresh):
    global arguments, blocks, keep_nonbarcodes, aligners, decoders, thresh
    arguments = aarguments
    blocks = ablocks
    keep_nonbarcodes = akeep_nonbarcodes
    thresh = athresh
    aligners = tuple(misc.build_bc_aligners(bblocks, unknown_read_orientation) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

def worker_build_aligners_single_end(aarguments, ablocks, akeep_nonbarcode, athresh):
    global arguments, blocks, keep_nonbarcode, aligners, decoders, thresh
    arguments = aarguments
    blocks = ablocks
    keep_nonbarcode = akeep_nonbarcode
    thresh = athresh
    aligners = misc.build_bc_aligners(blocks, unknown_read_orientation=False)
    decoders = misc.build_bc_decoders(blocks)

def worker_process_read(bc_recs):
    bc_rec1, bc_rec2 = bc_recs

    scores = [0, 0]
    sans_bc_rec = [None, None]
    tags = [None, None]
    read_qcs = [None, None]

    for i, (bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders) in enumerate(
        zip(blocks, keep_nonbarcodes, (bc_rec1, bc_rec2), aligners, decoders, strict=False)
    ):
        read_first_scores_recs_tags_plus_qc = process_bc_rec(arguments, bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders) # get our tags for each read and trim the read
        read_qc = read_first_scores_recs_tags_plus_qc[-1] # get our QC (last item of the tuple)
        read_qcs[i] = read_qc
        scores[i], sans_bc_rec[i], tags[i] = read_first_scores_recs_tags_plus_qc[:-1] # get the first three tuple entries

    output_recs = [None, None]
    output_tags = [None, None]
    if all(score >= cthresh for score, cthresh in zip(scores, thresh)) and all(sans_bc_rec):
        # TODO: don't write second fastq if barcodes only in first. Need to do this currently in case input
        # fastqs are gzipped: STAR input must be either all gzipped or all uncompressed
        for i, (bc_rec, csans_bc_rec, ctags) in enumerate(zip((bc_rec1, bc_rec2), sans_bc_rec, tags)):
            output_recs[i] = csans_bc_rec if csans_bc_rec and csans_bc_rec is not True else bc_rec
            output_tags[i] = ctags

    r1_meta = qc_metrics.get_config_metadata(arguments, "barcode_struct_r1")
    r2_meta = qc_metrics.get_config_metadata(arguments, "barcode_struct_r2") if len(read_qcs) > 1 else []

    r1_ambiguous_rec, r1_no_match_rec = build_qc_failure_records(bc_rec1, read_qcs[0], r1_meta)

    r2_ambiguous_rec, r2_no_match_rec = (None, None)
    if len(blocks) > 1 and read_qcs[1] is not None:
        r2_ambiguous_rec, r2_no_match_rec = build_qc_failure_records(bc_rec2, read_qcs[1], r2_meta)

    qc_failure_recs = {
        "r1_ambiguous": r1_ambiguous_rec,
        "r1_no_match": r1_no_match_rec,
        "r2_ambiguous": r2_ambiguous_rec,
        "r2_no_match": r2_no_match_rec,
    }

    return output_recs, output_tags, read_qcs, qc_failure_recs

def worker_process_read_single_end(bc_rec):
    read_first_scores_recs_tags_plus_qc = process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders) # get our tags for each read and trim the read
    read_qc = read_first_scores_recs_tags_plus_qc[-1] # get our QC (last item of the tuple)
    score, sans_bc_rec, tags = read_first_scores_recs_tags_plus_qc[:-1] # get the first three tuple entries

    output_rec = None
    output_tags = None

    meta_list = qc_metrics.get_config_metadata(arguments, "barcode_struct_r1")
    ambiguous_rec, no_match_rec = build_qc_failure_records(
        bc_rec,
        read_qc,
        meta_list,
    )

    if score >= thresh and sans_bc_rec:
        output_rec = sans_bc_rec if sans_bc_rec is not True else bc_rec
        output_tags = tags

    return output_rec, output_tags, read_qc, ambiguous_rec, no_match_rec

def read_iterator(fq1_fpath, fq2_fpath, sem):
    for bc_rec1, bc_rec2 in zip(SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq'), SeqIO.parse(misc.gzip_friendly_open(fq2_fpath), 'fastq')):
        sem.acquire()
        yield bc_rec1, bc_rec2

def read_iterator_single_end(fq1_fpath, sem):
    for bc_rec1 in SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq'):
        sem.acquire()
        yield bc_rec1

def parallel_process_fastqs(arguments, fq1_fpath, fq2_fpath, sans_bc_fq1_fpath, sans_bc_fq2_fpath, tags1_fpath, tags2_fpath, bcs_on_both_reads):
    chunksize = 100

    blocks = [arguments.config["barcode_struct_r1"]["blocks"]]
    keep_nonbarcodes = [arguments.config["barcode_struct_r1"]["keep_nonbarcode"]]
    if bcs_on_both_reads:
        blocks.append(arguments.config["barcode_struct_r2"]["blocks"])
        keep_nonbarcodes.append(arguments.config["barcode_struct_r2"]["keep_nonbarcode"])

    log.info('Building aligners and barcode decoders')
    aligners = tuple(misc.build_bc_aligners(bblocks, arguments.config["unknown_read_orientation"]) for bblocks in blocks)
    decoders = tuple(misc.build_bc_decoders(bblocks) for bblocks in blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    thresh = [np.inf, -np.inf]
    for i, (fq_fpath, bblocks, bkeep_nonbarcodes, baligners, bdecoders) in enumerate(zip((fq1_fpath, fq2_fpath), blocks, keep_nonbarcodes, aligners, decoders, strict=False)):
        first_scores_recs_tags = []
        for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq_fpath), 'fastq')):
            first_scores_recs_tags.append(process_bc_rec(arguments, bblocks, bkeep_nonbarcodes, bc_rec, baligners, bdecoders)[:-1]) # get the first three tuple entries
            if j >= n_first_seqs:
                break

        scores = [score for score, rec, tags in first_scores_recs_tags]
        thresh[i] = np.average(scores) - 2 * np.std(scores)

    log.info(f'Score threshold: {thresh[0]:.2f}{f"{thresh[1]:.2f}" if thresh[1] != -np.inf else ""}')

    r1_qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq1_fpath)
    r2_qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq2_fpath) if bcs_on_both_reads else None

    all_reads_qcs = []
    # multiprocessing.pool.imap tends to prefetch too many values into memory. The semaphore prevents that.
    sem = BoundedSemaphore(3 * chunksize * arguments.threads)
    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            (open(sans_bc_fq2_fpath, 'w') if sans_bc_fq2_fpath else nullcontext(None)) as bc_fq2_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            (open(tags2_fpath, 'w') if bcs_on_both_reads else nullcontext(None)) as tag2_fh, \
            open(r1_qc_fastq_paths["ambiguous"], "w") as r1_ambiguous_fh, \
            open(r1_qc_fastq_paths["no_match"], "w") as r1_no_match_fh, \
            (open(r2_qc_fastq_paths["ambiguous"], "w") if bcs_on_both_reads else nullcontext(None)) as r2_ambiguous_fh, \
            (open(r2_qc_fastq_paths["no_match"], "w") if bcs_on_both_reads else nullcontext(None)) as r2_no_match_fh, \
            Pool(arguments.threads, initializer=worker_build_aligners, initargs=(arguments, blocks, keep_nonbarcodes, arguments.config["unknown_read_orientation"], thresh)) as pool:

        log.info('Continuing...')
        total_out = 0
        for i, (output_recs, output_tags, read_qcs, qc_failure_recs) in enumerate(
            pool.imap(worker_process_read, read_iterator(fq1_fpath, fq2_fpath, sem), chunksize=chunksize)
        ):
            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d} processed,  {total_out:,d} output')
            sem.release()

            all_reads_qcs.append(tuple(read_qcs[:len(blocks)]))

            if qc_failure_recs["r1_ambiguous"] is not None:
                SeqIO.write(qc_failure_recs["r1_ambiguous"], r1_ambiguous_fh, "fastq")
            if qc_failure_recs["r1_no_match"] is not None:
                SeqIO.write(qc_failure_recs["r1_no_match"], r1_no_match_fh, "fastq")

            if bcs_on_both_reads:
                if qc_failure_recs["r2_ambiguous"] is not None:
                    SeqIO.write(qc_failure_recs["r2_ambiguous"], r2_ambiguous_fh, "fastq")
                if qc_failure_recs["r2_no_match"] is not None:
                    SeqIO.write(qc_failure_recs["r2_no_match"], r2_no_match_fh, "fastq")

            processed = False
            for rec, bc_fh, tags, tag_fh in zip(output_recs, (bc_fq1_fh, bc_fq2_fh), output_tags, (tag1_fh, tag2_fh)):
                if rec and bc_fh:
                    processed = True
                    SeqIO.write(rec, bc_fh, 'fastq')
                if tags:
                    output_rec_name_and_tags(rec, tags, tag_fh)
            total_out += processed
    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} pairs of records output')
    return all_reads_qcs

def parallel_process_fastqs_single_end(arguments, fq1_fpath, sans_bc_fq1_fpath, tags1_fpath):
    chunksize = 100

    blocks = arguments.config["barcode_struct_r1"]["blocks"]
    keep_nonbarcode = [arguments.config["barcode_struct_r1"]["keep_nonbarcode"]]

    log.info('Building aligners and barcode decoders')
    aligners = misc.build_bc_aligners(blocks, unknown_read_orientation=False)
    decoders = misc.build_bc_decoders(blocks)

    log.info(f'Processing first {n_first_seqs:,d} for score threshold...')

    first_scores_recs_tags = []
    for j, bc_rec in enumerate(SeqIO.parse(misc.gzip_friendly_open(fq1_fpath), 'fastq')):
        first_scores_recs_tags.append(process_bc_rec(arguments, blocks, keep_nonbarcode, bc_rec, aligners, decoders)[:-1])# get the first three tuple entries
        if j >= n_first_seqs:
            break

    scores = [score for score, rec, tags in first_scores_recs_tags]
    thresh = np.average(scores) - 2 * np.std(scores)

    log.info(f'Score threshold: {thresh:.2f}')

    qc_fastq_paths = get_per_fastq_qc_fastq_paths(arguments, fq1_fpath)

    # multiprocessing.pool.imap tends to prefetch too many values into memory.
    # The semaphore prevents that.
    sem = BoundedSemaphore(3 * chunksize * arguments.threads)
    with (open(sans_bc_fq1_fpath, 'w') if sans_bc_fq1_fpath else nullcontext(None)) as bc_fq1_fh, \
            open(tags1_fpath, 'w') as tag1_fh, \
            open(qc_fastq_paths["ambiguous"], "w") as ambiguous_fh, \
            open(qc_fastq_paths["no_match"], "w") as no_match_fh, \
            Pool(
                arguments.threads,
                initializer=worker_build_aligners_single_end,
                initargs=(arguments, blocks, keep_nonbarcode, thresh)
            ) as pool:
        log.info('Continuing...')
        total_out = 0
        all_reads_qcs = []

        for i, (output_rec, output_tags, read_qc, ambiguous_rec, no_match_rec) in enumerate(
            pool.imap(worker_process_read_single_end, read_iterator_single_end(fq1_fpath, sem), chunksize=chunksize)
        ):
            if i % 100000 == 0 and i > 0:
                log.info(f'  {i:,d} processed,  {total_out:,d} output')
            sem.release()

            all_reads_qcs.append(read_qc)

            if ambiguous_rec is not None:
                SeqIO.write(ambiguous_rec, ambiguous_fh, "fastq")
            if no_match_rec is not None:
                SeqIO.write(no_match_rec, no_match_fh, "fastq")

            processed = False
            if output_rec and bc_fq1_fh:
                processed = True
                SeqIO.write(output_rec, bc_fq1_fh, 'fastq')
            if output_tags:
                output_rec_name_and_tags(output_rec, output_tags, tag1_fh)

            total_out += processed

    log.info(f'{i+1:,d} barcode records processed')
    log.info(f'{total_out:,d} records output')
    return all_reads_qcs

def count_parallel_wrapper(ref_and_input_bam_fpath):
    ref, input_bam_fpath = ref_and_input_bam_fpath
    read_count_given_bc_then_feature_then_umi = defaultdict(lambda : defaultdict(Counter))
    for read in pysam.AlignmentFile(input_bam_fpath).fetch(ref):
        if not read.is_paired or read.is_read1 or (read.is_read2 and read.mate_is_unmapped):
            for gx_gn_tup in misc.gx_gn_tups_from_read(read): # count read toward all compatible genes
                read_count_given_bc_then_feature_then_umi[read.get_tag('CB')][gx_gn_tup][read.get_tag('UB')] += 1
    for k in read_count_given_bc_then_feature_then_umi.keys():
        read_count_given_bc_then_feature_then_umi[k] = dict(read_count_given_bc_then_feature_then_umi[k])
    read_count_given_bc_then_feature_then_umi = dict(read_count_given_bc_then_feature_then_umi)
    return ref, read_count_given_bc_then_feature_then_umi


def count_matrix(arguments, input_bam_fpath):
    """
    Counts the reads from the input bam file and outputs sparse matrices of read and UMI counts.
    """
    raw_reads_output_dir = os.path.join(arguments.output_dir, 'raw_reads_bc_matrix')
    raw_umis_output_dir = os.path.join(arguments.output_dir, 'raw_umis_bc_matrix')
    for out_dir in [raw_reads_output_dir, raw_umis_output_dir]:
        if os.path.exists(out_dir):
            log.info('Matrix output folder exists. Skipping count matrix build')
            return
        else:
            os.makedirs(out_dir)

    log.info('Finding all barcodes present...')
    sorted_complete_bcs, sorted_features = misc.get_bcs_and_features_from_bam(
            input_bam_fpath,
            arguments.threads-1
            )
    i_given_feature = {feat: i for i, feat in enumerate(sorted_features)}
    j_given_complete_bc = {comp_bc: j for j, comp_bc in enumerate(sorted_complete_bcs)}


    with pysam.AlignmentFile(input_bam_fpath) as bamfile:
        reference_names = bamfile.references
    reference_names_with_input_bam = [(ref, input_bam_fpath) for ref in reference_names]


    log.info('Counting reads...')
    M_reads = lil_matrix((len(sorted_features), len(sorted_complete_bcs)), dtype=int)
    M_umis = lil_matrix((len(sorted_features), len(sorted_complete_bcs)), dtype=int)
    with Pool(arguments.threads) as pool:
        for ref, read_count_given_bc_then_feature_then_umi in pool.imap_unordered(
                count_parallel_wrapper,
                reference_names_with_input_bam):
            for comp_bc, read_count_given_feature_then_umi in read_count_given_bc_then_feature_then_umi.items():
                j = j_given_complete_bc[comp_bc]
                for gx_gn_tup, umi_cntr in read_count_given_feature_then_umi.items():
                    i = i_given_feature[gx_gn_tup]
                    for umi, count in umi_cntr.items():
                        M_reads[i, j] += count
                        M_umis[i, j] += 1

    log.info('Writing raw read count matrix...')
    for out_dir, M in [(raw_reads_output_dir, M_reads), (raw_umis_output_dir, M_umis)]:
        misc.write_matrix(M, sorted_complete_bcs, sorted_features, out_dir)
