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
import datetime

# Construct paths to pipeline scripts as constants
PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))

BARCODE_CORRECTER = os.path.join(PIPELINE_PATH, 'src/barcode_correct_scatac.jl')
TRIMMOMATIC = os.path.join(PIPELINE_PATH,
    'Trimmomatic-0.36/trimmomatic-0.36.jar')

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
    parser.add_argument('--keep_intermediates',
        action='store_true', help='Skip clean up steps to keep intermediate '
        'files.')
    args = parser.parse_args()

    # Initialize directories and file prefixes required for the pipeline
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
    FASTQ_DIRECTORY = os.path.join(args.outdir, 'fastqs')
    QC_DIRECTORY = os.path.join(args.outdir, 'qc_info')
    print(FASTQ_DIRECTORY)
    for directory in [FASTQ_DIRECTORY, QC_DIRECTORY]:
        if not os.path.exists(directory):
            os.mkdir(directory)

    bcl_out1 = FASTQ_DIRECTORY + '/Undetermined_S0_R1_001.fastq.gz'
    bcl_out2 = FASTQ_DIRECTORY + '/Undetermined_S0_R2_001.fastq.gz'
    bar_out1 = FASTQ_DIRECTORY + '/Undetermined_S0_R1_001.fastq.gz.out.fq.gz'
    bar_out2 = FASTQ_DIRECTORY + '/Undetermined_S0_R2_001.fastq.gz.out.fq.gz'
    trimmer_out1 = OUTPUT_PREFIX + '.1.trimmed.paired.fastq.gz'
    trimmer_out2 = OUTPUT_PREFIX + '.2.trimmed.paired.fastq.gz'
    trimmer_un_out1 = OUTPUT_PREFIX + '.1.trimmed.unpaired.fastq.gz'
    trimmer_un_out2 = OUTPUT_PREFIX + '.2.trimmed.unpaired.fastq.gz'
    qc_info = QC_DIRECTORY + "/QC_stats.txt"
    qcf = open(qc_info, 'a')
    qcf.write("Runall QC info for " + str(args.rundir) + "\n%s" % datetime.datetime.now())
    qcf.close()
    # Configure logger
    logging.basicConfig(filename= OUTPUT_PREFIX + '.log',format='%(asctime)s '
        '%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('Pipeline started.')

    # Submit bcl2fastq only if no existing results or if user wants to
    # overwrite
    if not os.path.exists(bcl_out1) or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        logging.info('bcl2fastq started.')
	bcl2fastq_command = ('module load modules modules-init modules-gs '
            'bcl2fastq/2.16 fastqc/0.10.1; bcl2fastq --runfolder-dir %s '
            '-o %s --no-lane-splitting --interop-dir %s --sample-sheet %s/Sample_sheet.csv ' % (args.rundir,
            FASTQ_DIRECTORY, FASTQ_DIRECTORY, FASTQ_DIRECTORY))
        f = open(os.path.join(args.outdir, "bcl2fastq_log.txt"), 'w')
        subprocess.check_call(bcl2fastq_command, shell=True, stdout = f,
            stderr=f)

    else:
        logging.info('bcl2fastq skipped.')

    # Submit barcode corrector only if no existing results or if user wants
    # to overwrite
    if not os.path.exists(bar_out1) or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
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
        logging.info('Barcode corrector skipped.')

    # Submit trimmer only if no existing results or if user wants
    # to overwrite
    if not os.path.exists(trimmer_out1) or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        qcf = open(qc_info, 'a')
        qcf.write("\n\nTrimmomatic:\n\n")
        qcf.flush()
        logging.info('Trimmomatic started.')
        trimmer_command = ('java -Xmx1G -jar %s PE -threads %s %s %s %s '
            '%s %s %s ILLUMINACLIP:%s'
            '/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true '
            'MINLEN:20' % (TRIMMOMATIC, args.nthreads, bar_out1, bar_out2,
            trimmer_out1, trimmer_un_out1, trimmer_out2, trimmer_un_out2,
            PIPELINE_PATH))
        subprocess.check_call(trimmer_command, shell=True, stdout=qcf, stderr=qcf)
        logging.info('Trimmomatic ended.')
        qcf.close()
    else:
        logging.info('Trimmomatic skipped.')

    # Submit bowtie mapping only if no existing results or if user wants
    # to overwrite
    if not os.path.exists(OUTPUT_PREFIX + ".bam") or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        qcf = open(qc_info, 'a')
        qcf.write("\n\nBowtie:\n\n")
        qcf.flush()
        logging.info('Bowtie2 started.')
	subprocess.check_call('module load bowtie2/latest; module load samtools/latest; bowtie2 -3 1 --un-conc-gz %s.unaligned.fq.gz -X 2000 -p %s '
            '-x %s -1 %s -2 %s | samtools view -Sb - > %s.bam' %
            (OUTPUT_PREFIX, args.nthreads, args.genome, trimmer_out1,
            trimmer_out2, OUTPUT_PREFIX), shell=True, stderr=qcf, stdout=qcf)
        logging.info('Bowtie2 ended.')
        qcf.close()
    else:
        logging.info('Bowtie2 skipped.')

    if not os.path.exists(OUTPUT_PREFIX + ".split.q10.sort.bam") or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        logging.info('Quality filter started.')
        qcf = open(qc_info, 'a')
        qcf.write("\n\nQualiy control:\n\nTotal mitochondrial reads: ")
        qcf.flush()
        subprocess.call("module load samtools/latest; samtools view %s.bam | grep -c '\tchrM\t' -" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
        subprocess.check_call("module load samtools/latest; samtools view -h -f3 -F12 -q10 %s.bam | grep "
            " -v '\tchrM\t' | samtools sort -T %s.sorttemp -@ %s - -o "
            "%s.split.q10.sort.bam" % (OUTPUT_PREFIX, OUTPUT_PREFIX,
            args.nthreads, OUTPUT_PREFIX), shell=True)
        subprocess.check_call('module load samtools/latest; samtools index %s.split.q10.sort.bam' % OUTPUT_PREFIX,
            shell=True)
        qcf.write("\nTotal reads after QC10: ")
        qcf.flush()
        subprocess.call("module load samtools/latest; samtools view -c -f3 -F12 %s.split.q10.sort.bam" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
        qcf.write("\n\n")
        qcf.close()
        logging.info('Quality filter finished.')
    else:
        logging.info('Quality filter skipped.')

    if not args.keep_intermediates:
    # Remove temporary files created during the pipeline.
        clean_command = ('rm %s; rm %s; rm %s; rm %s; rm %s; rm %s; rm %s.bam;' %
            (bar_out1, bar_out2, trimmer_un_out1, trimmer_un_out2, trimmer_out1, trimmer_out2, OUTPUT_PREFIX ))
        subprocess.check_call(clean_command, shell=True)

    logging.info('Runall complete.')
    qcf.close()
