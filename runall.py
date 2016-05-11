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

# Construct paths to pipeline scripts as constants
PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))

BARCODE_CORRECTER = os.path.join(PIPELINE_PATH, 'barcode_correct_scatac.py')
DEDUPLICATER = os.path.join(PIPELINE_PATH, 'sc_atac_true_dedup.py')
TRIMMOMATIC = os.path.join(PIPELINE_PATH,
    'Trimmomatic-0.36/trimmomatic-0.36.jar')

# Function to make necessary sub-directories in output folder
def initialize_directories(pipeline_root_directory):
    fastqs = os.path.join(args.outdir, 'fastqs')

    for directory in [fastqs]:
        if not os.path.exists(directory):
            os.mkdir(directory)

# main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to convert NextSeq'
        ' BCL files to cleaned and corrected fastq files for scATAC-seq '
        'analysis.')
    parser.add_argument('-R','--rundir', help='Run directory containing BCL '
        'files', dest='rundir')
    parser.add_argument('-O','--outdir', help='Output directory',   
        dest='outdir')
    parser.add_argument('-P','--prefix',help='Output file prefix, otherwise '
        'default(out)', default = "out", dest='prefix')
    parser.add_argument('-E','--maxedit', help='Maximum allowed edit distance '
        '(default = 3)', default=3, dest='maxedit')
    parser.add_argument('-G','--genome', help='Path to genome annotation '
        '(default=/net/trapnell/vol1/genomes/human/hg19/bowtie2/)',
        default='/net/trapnell/vol1/genomes/human/hg19/bowtie2/hg19',
        dest='genome')
    parser.add_argument('-p','--nthreads', help='Number of cores available '
        '(default=1)', default=1, dest='nthreads')
    parser.add_argument('--force_overwrite_all', action='store_true',
        help='Force overwrite of all steps of pipeline regardless of files '
        'already present.')
    parser.add_argument('--force_overwrite_bcl2fastq', action='store_true',
        help='Force overwrite bcl2fastq regardless of files already present.')
    parser.add_argument('--force_overwrite_cleanfastq', action='store_true',
        help='Force overwrite fastq cleaning regardless of files already
        present.')
    parser.add_argument('--force_overwrite_barcodecorrect',
        action='store_true', help='Force overwrite barcode correction '
        'regardless of files already present.')
    parser.add_argument('--force_overwrite_trimming', action='store_true',
        help='Force overwrite trimming regardless of files already present.')
    parser.add_argument('--force_overwrite_mapping', action='store_true',
        help='Force overwrite mapping regardless of files already present.')
    args = parser.parse_args()

    # Initialize directories required for the pipeline
    print 'Initializing directories for pipeline...'
    initialize_directories(args.outdir)

    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
    FASTQ_DIRECTORY = os.path.join(args.outdir, 'fastqs')

    bar_out1 = OUTPUT_PREFIX + 'split.1.fq.gz'
    bar_out2 = OUTPUT_PREFIX + 'split.2.fq.gz'

    # Submit bcl2fastq only if no existing results or user wants to overwrite
    if len(os.listdir(FASTQ_DIRECTORY)) == 0 or
        args.force_overwrite_all or
        args.force_overwrite_bcl2fastq:
        print 'Starting bcl2fastq...'
        bcl2fastq_command = 'module load modules modules-init modules-gs '
            'bcl2fastq/2.16 fastqc/0.10.1; bcl2fastq --runfolder-dir %s -o %s '
            '--ignore-missing-filter' % (args.rundir, FASTQ_DIRECTORY)
        f = open(os.path.join(args.outdir, "bcl2fastq_log.txt"), 'w')
        subprocess.call(bcl2fastq_command, shell=True, stdout = f, stderr=f)
    else:
        print 'fastq directory is not empty, skipping bcl2fastq. Specify '
            '--force_overwrite_all or --force_overwrite_bcl2fastq to redo.'

    print "Cleaning and fixing barcodes..."

    if not os.path.exists(bar_out1) or args.force_overwrite_all or args.force_overwrite_barcodecorrect:
        subprocess.call('python %s -F %s -o %s -E %s -n %s' % (BARCODE_CORRECTER, FASTQ_DIRECTORY, OUTPUT_PREFIX, args.maxedit, args.nthreads), shell=True)
    else:
        print 'Barcodes already fixed, skipping fix barcodes. Specify --force_overwrite_all or --force_overwrite_barcodecorrect to redo.'

    print "Trimming adapters..."

    trimmer_out1 = OUTPUT_PREFIX + '.split.1.trimmed.paired.fastq.gz'
    trimmer_out2 = OUTPUT_PREFIX + '.split.2.trimmed.paired.fastq.gz'
    trimmer_un_out1 = OUTPUT_PREFIX + '.split.1.trimmed.unpaired.fastq.gz'
    trimmer_un_out2 = OUTPUT_PREFIX + '.split.2.trimmed.unpaired.fastq.gz'

    if not os.path.exists(trimmer_out1) or args.force_overwrite_all or args.force_overwrite_trimming:
        trimmer_command = 'java -Xmx1G -jar -threads %s %s PE %s %s %s %s %s %s ILLUMINACLIP:%sTrimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:20' % \
            (args.nthreads, TRIMMOMATIC, bar_out1, bar_out2, trimmer_out1, trimmer_un_out1, trimmer_out2, trimmer_un_out2, PIPELINE_PATH)

        subprocess.call(trimmer_command, shell=True)
    else:
        print 'Sequences already trimmed, skipping trimming. Specify --force_overwrite_all or --force_overwrite_trimming to redo.'

    print "Cleaning up..."

    clean_command = 'rm %s; rm %s; rm %s; rm %s; ' % (bar_out1, bar_out2, trimmer_un_out1, trimmer_un_out2)
    subprocess.call(clean_command, shell=True)

    print "Starting mapping..."
    if not os.path.exists(OUTPUT_PREFIX + ".split.bam") or args.force_overwrite_all or args.force_overwrite_mapping:
        subprocess.call('bowtie2 -X 1000 -p %s -x %s -1 %s -2 %s | samtools view -Sb - > %s .split.bam' % (args.nthreads, args.genome, trimmer_out1, trimmer_out2, OUTPUT_PREFIX), shell=True)
    else:
        print 'Sequences already mapped, skipping mapping. Specify --force_overwrite_all or --force_overwrite_mapping to redo.'

    print "Mapping complete..."
    print "Deduplicating..."

    if not os.path.exists(OUTPUT_PREFIX + ".true.nodups.bam"):
        subprocess.call("samtools view -h -f3 -F12 -q10 %s.split.bam | grep -v '[0-9]'$'\t'chrM | samtools view -u - | samtools sort - %s.split.q10.sort.bam" % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.call('samtools index %s.split.q10.sort.bam' % OUTPUT_PREFIX, shell=True)
        subprocess.call('python %s %s.split.q10.sort.bam %s.true.nodups.bam' % (DEDUPLICATER, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.call('samtools view %s.true.nodups.bam | sort -u -k1,1 | cut -f9 > %s.insertsize.txt' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.call('samtools index %s.true.nodups.bam' % OUTPUT_PREFIX, shell=True)
    else:
        print 'Sequences already deduplicating, skipping.'

    print "Complete."
