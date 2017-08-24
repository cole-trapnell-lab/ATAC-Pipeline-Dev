#!/usr/bin/python
# Python pipeline to convert NextSeq BCL files to fastq,
# fix errors in barcodes and trim adapters from fastqs.

# Cannibalized from Darren with a bunch of stuff stolen from Andrew and
# compiled/maintained by Hannah


import argparse
import os
import subprocess
import gzip
import os.path
import sys
import logging

# Construct paths to pipeline scripts as constants
PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))

BARCODE_CORRECTER = os.path.join(PIPELINE_PATH, 'src/barcode_correct_scatac.jl')
TRIMMOMATIC = os.path.join(PIPELINE_PATH,
    'Trimmomatic-0.36/trimmomatic-0.36.jar')


def initialize_directories(pipeline_root_directory):
    """Initialize necessary sub-directories in output folder"""
    fastqs = os.path.join(args.outdir, 'fastqs')

    for directory in [fastqs]:
        if not os.path.exists(directory):
            os.mkdir(directory)


# main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to convert'
        ' BCL files to cleaned, corrected and mapped BAM files for sci-ATAC-seq '
        'analysis.')
    parser.add_argument('-R','--rundir', help='Run directory containing BCL '
        'files', dest='rundir', required=True)
    parser.add_argument('-O','--outdir', help='Output directory',
        dest='outdir', required=True)
    parser.add_argument('-P','--prefix',help='Output file prefix, otherwise '
        'default(out)', default = "out", dest='prefix', required=False)
    parser.add_argument('--miseq', help='Run on MiSeq, otherwise assumes '
        'NextSeq', action='store_true', required=False)
    parser.add_argument('-E','--maxedit', help='Maximum allowed edit distance '
        '(default = 3)', default=3, dest='maxedit', required=False)
    parser.add_argument('-G','--genome', help='Path to genome annotation '
        '(default=/net/trapnell/vol1/genomes/human/hg19/bowtie2/)',
        default='/net/trapnell/vol1/genomes/human/hg19/bowtie2/hg19',
        dest='genome', required=False)
    parser.add_argument('-p','--nthreads', help='Number of cores available '
        '(default=1)', default=1, dest='nthreads', required=False)
    parser.add_argument('--force_overwrite_all', action='store_true',
        help='Force overwrite of all steps of pipeline regardless of files '
        'already present.')
    parser.add_argument('--force_overwrite_bcl2fastq', action='store_true',
        help='Force overwrite bcl2fastq regardless of files already present.')
    parser.add_argument('--force_overwrite_cleanfastq', action='store_true',
        help='Force overwrite fastq cleaning regardless of files already '
        'present.')
    parser.add_argument('--force_overwrite_barcodecorrect',
        action='store_true', help='Force overwrite barcode correction '
        'regardless of files already present.')
    parser.add_argument('--force_overwrite_trimming', action='store_true',
        help='Force overwrite trimming regardless of files already present.')
    parser.add_argument('--force_overwrite_mapping', action='store_true',
        help='Force overwrite mapping regardless of files already present.')
    parser.add_argument('--keep_intermediates',
        action='store_true', help='Skip clean up steps to keep intermediate '
        'files.')
    args = parser.parse_args()

    # Initialize directories and file prefixes required for the pipeline
    initialize_directories(args.outdir)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
    FASTQ_DIRECTORY = os.path.join(args.outdir, 'fastqs')

    bcl_out1 = FASTQ_DIRECTORY + '/Undetermined_S0_R1_001.fastq.gz'
    bcl_out2 = FASTQ_DIRECTORY + '/Undetermined_S0_R2_001.fastq.gz'
    bar_out1 = FASTQ_DIRECTORY + '/Undetermined_S0_R1_001.fastq.gz.out.fq.gz'
    bar_out2 = FASTQ_DIRECTORY + '/Undetermined_S0_R2_001.fastq.gz.out.fq.gz'
    trimmer_out1 = OUTPUT_PREFIX + '.1.trimmed.paired.fastq.gz'
    trimmer_out2 = OUTPUT_PREFIX + '.2.trimmed.paired.fastq.gz'
    trimmer_un_out1 = OUTPUT_PREFIX + '.1.trimmed.unpaired.fastq.gz'
    trimmer_un_out2 = OUTPUT_PREFIX + '.2.trimmed.unpaired.fastq.gz'

    # Configure logger
    logging.basicConfig(filename= OUTPUT_PREFIX + '.log',format='%(asctime)s '
        '%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('Pipeline started.')

    # Only continue with initial cleaning process if cleaned and trimmed files
    # do not exist, or user specifies. This line exists because after the clean
    # step, the intermediate files checked for below are removed.
    if not os.path.exists(trimmer_out1) or \
        args.force_overwrite_all or \
        args.force_overwrite_trimming:

	# Submit bcl2fastq only if no existing results or if user wants to
        # overwrite
        if len(os.listdir(FASTQ_DIRECTORY)) == 0 or \
            args.force_overwrite_all or \
            args.force_overwrite_bcl2fastq:

            print('Starting bcl2fastq...')
            logging.info('bcl2fastq started.')
	    bcl2fastq_command = ('module load modules modules-init modules-gs '
                'bcl2fastq/2.16 fastqc/0.10.1; bcl2fastq --runfolder-dir %s '
                '-o %s --no-lane-splitting --interop-dir %s --sample-sheet %s/Sample_sheet.csv ' % (args.rundir,
                FASTQ_DIRECTORY, FASTQ_DIRECTORY, FASTQ_DIRECTORY))
            f = open(os.path.join(args.outdir, "bcl2fastq_log.txt"), 'w')
            subprocess.check_call(bcl2fastq_command, shell=True, stdout = f,
                stderr=f)
            logging.info('bcl2fastq ended.')

        else:

            print ('fastq directory is not empty, skipping bcl2fastq. Specify '
                '--force_overwrite_all or --force_overwrite_bcl2fastq to '
                'redo.')
            logging.info('bcl2fastq skipped.')

        # Submit barcode corrector only if no existing results or if user wants
        # to overwrite
        if not os.path.exists(bar_out1) or \
            args.force_overwrite_all or \
            args.force_overwrite_barcodecorrect:

            print "Cleaning and fixing barcodes..."
            logging.info('Barcode corrector started.')
            if args.miseq:
                subprocess.check_call('julia %s -f %s -o %s -e %s --miseq' %
                    (BARCODE_CORRECTER, FASTQ_DIRECTORY, OUTPUT_PREFIX,
                    args.maxedit), shell=True)
            else:
                subprocess.check_call('julia %s -f %s -o %s -e %s' %
                    (BARCODE_CORRECTER, FASTQ_DIRECTORY, OUTPUT_PREFIX,
                    args.maxedit), shell=True)
            logging.info('Barcode corrector ended.')

        else:

            print ('Barcodes already fixed, skipping fix barcodes. Specify '
                '--force_overwrite_all or --force_overwrite_barcodecorrect to '
                'redo.')
            logging.info('Barcode corrector skipped.')

        # Submit trimmer only if no existing results or if user wants
        # to overwrite
        if not os.path.exists(trimmer_out1) or \
            args.force_overwrite_all or \
            args.force_overwrite_trimming:

            print "Trimming adapters..."
            logging.info('Trimmomatic started.')
            trimmer_command = ('java -Xmx1G -jar %s PE -threads %s %s %s %s '
                '%s %s %s ILLUMINACLIP:%s'
                '/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true '
                'MINLEN:20' % (TRIMMOMATIC, args.nthreads, bar_out1, bar_out2,
                trimmer_out1, trimmer_un_out1, trimmer_out2, trimmer_un_out2,
                PIPELINE_PATH))
            subprocess.check_call(trimmer_command, shell=True)
            logging.info('Trimmomatic ended.')

        else:

            print ('Sequences already trimmed, skipping trimming. Specify '
                '--force_overwrite_all or --force_overwrite_trimming to redo.')
            logging.info('Trimmomatic skipped.')

        if not args.keep_intermediates:
            print "Cleaning up..."

            # Remove temporary files created during the pipeline.
            clean_command = ('rm %s; rm %s; rm %s; rm %s' %
                (bar_out1, bar_out2, trimmer_un_out1, trimmer_un_out2))
            subprocess.check_call(clean_command, shell=True)

    # Submit bowtie mapping only if no existing results or if user wants
    # to overwrite
    if not os.path.exists(OUTPUT_PREFIX + ".bam") or \
        args.force_overwrite_all or \
        args.force_overwrite_mapping:

        logging.info('Bowtie2 started.')
        print "Starting mapping..."
	subprocess.check_call('module load bowtie2/latest; bowtie2 -3 1 --un-conc-gz %s.unaligned.fq.gz -X 2000 -p %s '
            '-x %s -1 %s -2 %s | samtools view -Sb - > %s.bam' %
            (OUTPUT_PREFIX, args.nthreads, args.genome, trimmer_out1,
            trimmer_out2, OUTPUT_PREFIX), shell=True)
        logging.info('Bowtie2 ended.')
        print "Mapping complete..."

    else:

        print ('Sequences already mapped, skipping mapping. Specify '
            '--force_overwrite_all or --force_overwrite_mapping to redo.')
        logging.info('Bowtie2 skipped.')

    if not os.path.exists(OUTPUT_PREFIX + ".split.q10.sort.bam"):
        print "Filtering low quality reads..."
        logging.info('Quality filter started.')
        subprocess.check_call("samtools view -h -f3 -F12 -q10 %s.bam | grep "
            " -v '\tchrM\t' | samtools sort -T %s.sorttemp -@ %s - -o "
            "%s.split.q10.sort.bam" % (OUTPUT_PREFIX, OUTPUT_PREFIX,
            args.nthreads, OUTPUT_PREFIX), shell=True)
        subprocess.check_call('samtools index %s.split.q10.sort.bam' % OUTPUT_PREFIX,
            shell=True)
        logging.info('Quality filter finished.')
    else:
        print 'Sequences already filtered, skipping.'
        logging.info('Quality filter skipped.')

    print "Complete."
