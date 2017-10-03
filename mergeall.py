#!/usr/bin/python
import argparse
import os
import subprocess
import gzip
import os.path
import sys
import logging

PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))
DEDUPLICATER = os.path.join(PIPELINE_PATH, 'src/sc_atac_true_dedup.py')
HG19_BLACKLIST = os.path.join(PIPELINE_PATH, 'src/ENCFF001TDO.bed')
PICARD = os.path.join(PIPELINE_PATH, 'picard-tools-1.141/picard.jar')
MCLUST = os.path.join(PIPELINE_PATH, 'src/mclust_call.R')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to aggregate '
        'mapped bams from sequencing runs and run hotspot for scATAC-seq '
        'analysis.')
    parser.add_argument('-B','--bamlist', nargs='*', help='Paths to bams to '
        'combined', dest='bamlist', required=True)
    parser.add_argument('-O','--outdir', help='Output directory',
        dest='outdir', required=True)
    parser.add_argument('-P','--prefix',help='Output file prefix, otherwise '
        'default(out)', default = "out", dest='prefix', required=False)
    parser.add_argument('-C','--barcodes', help='Barcodes combinations allowed'
        ' in a text file.', required=True, dest='barcodes', default="None")
    parser.add_argument('--no_complexity', action='store_true',
        help='Add flag if you would like to skip running picard tools '
        'EstimateLibraryComplexity')
    parser.add_argument('--override_reads_per_cell', default = "mclust",
        dest = 'cell_count_cutoff',
        help='Number of reads per cell to count as valid, calculated using '
        'mclust if left blank')
    parser.add_argument('--force_overwrite_all', action='store_true',
        help='Force overwrite of all steps of pipeline regardless of files '
        'already present.')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
    QCDIR = os.path.join(args.outdir, 'qc_info')

    if not os.path.exists(QCDIR):
        os.mkdir(QCDIR)

    # Configure logger
    logging.basicConfig(filename= OUTPUT_PREFIX + '.log',format='%(asctime)s '
        '%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('Pipeline started.')

    # Make sorted bed file from all bams
    if not os.path.exists(OUTPUT_PREFIX + ".merge.bam") or \
        args.force_overwrite_all:
        if len(args.bamlist) == 1:
            print ('Only 1 bam, skipping merge.')
            mergename = args.bamlist[0]
        else:
            logging.info('Merge started.')
            subprocess.check_call('samtools merge %s.merge.bam %s' % (OUTPUT_PREFIX,
                ' '.join(args.bamlist)), shell=True)
            subprocess.check_call('samtools index %s.merge.bam' % OUTPUT_PREFIX,
                shell=True)
            logging.info('Merge ended.')
            mergename = OUTPUT_PREFIX + '.merge.bam'
    else:
        mergename = OUTPUT_PREFIX + '.merge.bam'
        print ('Bams already merged, skipping.')
        logging.info('Merge skipped.')

    if not args.no_complexity and not os.path.exists(OUTPUT_PREFIX + ".complexity_metrics.txt"):
        subprocess.check_call('java -jar %s EstimateLibraryComplexity '
        'I=%s O=%s.complexity_metrics.txt QUIET=true' % (PICARD,
        mergename, OUTPUT_PREFIX), shell=True)

    if not os.path.exists(OUTPUT_PREFIX + ".true.nodups.bam"):
        logging.info('Deduplication started.')
        subprocess.check_call('python %s %s %s.true.nodups.bam' %
            (DEDUPLICATER, mergename, OUTPUT_PREFIX), shell=True)
        subprocess.check_call('samtools view %s.true.nodups.bam | sort -u -k1,1 | '
            'cut -f9 > %s.insertsize.txt' % (OUTPUT_PREFIX, OUTPUT_PREFIX),
            shell=True)
        subprocess.check_call('samtools index %s.true.nodups.bam' % OUTPUT_PREFIX,
            shell=True)
        logging.info('Deduplication ended.')

    else:
        print ('Sequences already deduplicated, skipping.')
        logging.info('Deduplication skipped.')

# Remove blacklisted regions
    if not os.path.exists(OUTPUT_PREFIX + ".clean.bed") or \
        args.force_overwrite_all:
        logging.info('Read clean up started.')
        subprocess.check_call('''bedtools intersect -bed -a %s.true.nodups.bam  '''
            '''-b %s -v | sed 's/\/[0-9]//' > %s.clean.bed''' %
            (OUTPUT_PREFIX, HG19_BLACKLIST, OUTPUT_PREFIX),
            shell=True)
        subprocess.check_call('''grep -Fwf %s %s.clean.bed > %s.cleant.bed''' % (args.barcodes, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
	subprocess.check_call('mv %s.cleant.bed %s.clean.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        logging.info('Read clean up ended.')
    else:
        print ('Read clean up already done, skipping.')
        logging.info('Read clean up skipped.')

    # Remove low read count cells
    if not os.path.exists(OUTPUT_PREFIX + ".for_macs.bed") or \
        args.force_overwrite_all:
        logging.info('Remove low count cells started.')
        # Count cell reads
        subprocess.check_call("awk '{h[$4]++}; END { for(k in h) print k, h[k] }' "
            "%s.clean.bed > %s.cell_read_counts.txt" % (OUTPUT_PREFIX, OUTPUT_PREFIX),
            shell=True)
        if args.cell_count_cutoff == "mclust":
            # Exclude cells with less than n reads where n is determined by #mclust
            args.cell_count_cutoff = subprocess.call("Rscript --vanilla %s %s" % (MCLUST, OUTPUT_PREFIX),
                shell=True)

        subprocess.check_call("awk '{if ($2 > %s) print $1}' %s.cell_read_counts.txt > high_read_cells.txt"
            % (args.cell_count_cutoff, OUTPUT_PREFIX))
        subprocess.check_call('''grep -Fwf high_read_cells.txt %s.clean.bed | '''
            '''awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $6, $7} > %s.for_macs.bed'''
            % (OUTPUT_PREFIX, OUTPUT_PREFIX))

    if not os.path.exists(OUTPUT_PREFIX + ".clean.bed") or \
        args.force_overwrite_all:
        subprocess.check_call('''module load python/2.7.3; module load '''
            '''numpy/1.8.1; module load setuptools/25.1.1; module load '''
            '''MACS/2.1.0; macs2 callpeak -t %s.for_macs.bed --nomodel '''
            '''--keep-dup all --extsize 200 --shift -100 -f BED -g hs -n '''
            ''''%s_macs --call-summits''' % (OUTPUT_PREFIX, OUTPUT_PREFIX))
# count invalid barcodes

    logging.info('Mergeall Complete.')
