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
QCPLOTS = os.path.join(PIPELINE_PATH, 'src/qcplots.R')

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
    parser.add_argument('-R','--blacklist', help='Path to bed of blacklisted regions to be '
        'removed. Default is HG19 encode blacklist: ENCFF001TDO.bed. For no '
        'blacklist, put None.', required=False, dest='blacklist', default=HG19_BLACKLIST)
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
    parser.add_argument('--keep_intermediates',
        action='store_true', help='Skip clean up steps to keep intermediate '
        'files.')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
    QC_DIRECTORY = os.path.join(args.outdir, 'qc_info')
    MACS_DIRECTORY = os.path.join(args.outdir, 'macs_output')

    for directory in [QC_DIRECTORY, MACS_DIRECTORY]:
        if not os.path.exists(directory):
            os.mkdir(directory)

    qc_info = QC_DIRECTORY + "/QC_stats.txt"
    qcf = open(qc_info, 'a')
    qcf.write("Mergeall QC info for " + str(args.bamlist))
    qcf.close()

    # Configure logger
    logging.basicConfig(filename= OUTPUT_PREFIX + '.log',format='%(asctime)s '
        '%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('Pipeline started.')

    # Make sorted bed file from all bams
    qcf = open(qc_info, 'a')
    if len(args.bamlist) == 1:
        logging.info('Only 1 bam, skipping merge.')
        mergename = args.bamlist[0]
    else:
        if not os.path.exists(OUTPUT_PREFIX + ".merge.bam") or \
            args.force_overwrite_all:
	    args.force_overwrite_all = True
            logging.info('Merge started.')
            subprocess.check_call('samtools merge %s.merge.bam %s' % (OUTPUT_PREFIX,
                ' '.join(args.bamlist)), shell=True)
            subprocess.check_call('samtools index %s.merge.bam' % OUTPUT_PREFIX,
                shell=True)
            logging.info('Merge ended.')
            mergename = OUTPUT_PREFIX + '.merge.bam'
            qcf.write("\n\nTotal reads after merge: ")
            qcf.flush()
            subprocess.call("samtools view -c -f3 -F12 %s" % (mergename), shell=True, stdout=qcf, stderr=qcf)
        else:
            mergename = OUTPUT_PREFIX + '.merge.bam'
            logging.info('Merge skipped.')
    qcf.close()

    if not args.no_complexity and not os.path.exists(QC_DIRECTORY + "/" + args.prefix +  ".complexity_metrics.txt"):
        subprocess.check_call('java -jar %s EstimateLibraryComplexity '
        'I=%s O=%s/%s.complexity_metrics.txt QUIET=true' % (PICARD,
        mergename, QC_DIRECTORY, args.prefix), shell=True)

    if not os.path.exists(OUTPUT_PREFIX + ".true.nodups.bam") or \
	args.force_overwrite_all:
	args.force_overwrite_all = True
        logging.info('Deduplication started.')
        subprocess.check_call('python %s %s %s.true.nodups.bam' %
            (DEDUPLICATER, mergename, OUTPUT_PREFIX), shell=True)
        subprocess.check_call('samtools view %s.true.nodups.bam | sort -u -k1,1 | '
            'cut -f9 > %s/%s.insertsize.txt' % (OUTPUT_PREFIX, QC_DIRECTORY, args.prefix),
            shell=True)
        subprocess.check_call('samtools index %s.true.nodups.bam' % OUTPUT_PREFIX,
            shell=True)
        logging.info('Deduplication ended.')
        qcf = open(qc_info, 'a')
        qcf.write("\nTotal reads after deduplication: ")
        qcf.flush()
        subprocess.call("samtools view -c -f3 -F12 %s.true.nodups.bam" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
        qcf.close()

    else:
        logging.info('Deduplication skipped.')

    subprocess.call("Rscript --vanilla %s %s %s" % (QCPLOTS, QC_DIRECTORY, args.prefix),
         shell=True)

# Remove blacklisted regions
    if not os.path.exists(OUTPUT_PREFIX + ".clean.bed") or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        logging.info('Read clean up started.')
        if args.blacklist == None:
            logging.info('No blacklist provided, skipping.')
            clean_name = OUTPUT_PREFIX + ".true.nodups.bam"
        else:
            subprocess.check_call('''bedtools intersect -bed -a %s.true.nodups.bam  '''
                '''-b %s -v | sed 's/\/[0-9]//' > %s.clean.bed''' %
                (OUTPUT_PREFIX, args.blacklist, OUTPUT_PREFIX),
                shell=True)
            qcf = open(qc_info, 'a')
            qcf.write("\nTotal reads after removing blacklist: ")
            qcf.flush()
            subprocess.call("wc -l %s.clean.bed" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
            clean_name = OUTPUT_PREFIX + ".clean.bed"

# Remove non existant barcodes
        qcf.write("\nTotal reads after removing bad barcodes: ")
        qcf.flush()
        subprocess.check_call('''grep -Fwf %s %s > %s.cleant.bed''' % (args.barcodes, clean_name, OUTPUT_PREFIX), shell=True)
	subprocess.check_call('mv %s.cleant.bed %s.clean.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        logging.info('Read clean up ended.')
        subprocess.call("wc -l %s.clean.bed" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
        qcf.close()
    else:
        logging.info('Read clean up skipped.')

    # Remove low read count cells
    if not os.path.exists(OUTPUT_PREFIX + ".for_macs.bed") or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        logging.info('Remove low count cells started.')
        # Count cell reads
        subprocess.check_call("awk '{h[$4]++}; END { for(k in h) print k, h[k] }' "
            "%s.clean.bed > %s.cell_read_counts.txt" % (OUTPUT_PREFIX, OUTPUT_PREFIX),
            shell=True)
        if args.cell_count_cutoff == "mclust":
            # Exclude cells with less than n reads where n is determined by #mclust
            cutoff = subprocess.Popen("Rscript --vanilla %s %s" % (MCLUST, OUTPUT_PREFIX),
                shell=True, stdout=subprocess.PIPE)
            args.cell_count_cutoff = cutoff.stdout.readline().split(" ")[1].strip()
        qcf = open(qc_info, 'a')
        qcf.write("\nRead per cell cutoff: " + str(args.cell_count_cutoff))
        qcf.write("\n\nTotal reads after removing low read cells: ")
        qcf.flush()
        subprocess.check_call("awk '{if ($2 > %s) print $1}' %s.cell_read_counts.txt > %s.high_read_cells.txt"
            % (args.cell_count_cutoff, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.check_call('''grep -Fwf %s.high_read_cells.txt %s.clean.bed | '''
            '''awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4, $5, $6}' > %s.for_macs.bed'''
            % (OUTPUT_PREFIX, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.call("wc -l %s.for_macs.bed" % (OUTPUT_PREFIX), shell=True, stdout=qcf, stderr=qcf)
        qcf.write("\nTotal cells with more than " + args.cell_count_cutoff + " reads: ")
	qcf.flush()
	subprocess.call("wc -l %s.high_read_cells.txt" % OUTPUT_PREFIX, shell=True, stdout=qcf, stderr=qcf)
        qcf.write("\n\n")
        qcf.close()

    # Call peaks using MACS2
    if not os.path.exists(MACS_DIRECTORY + "/" + args.prefix + "_macs_peaks.narrowPeak") or \
        args.force_overwrite_all:
        args.force_overwrite_all = True
        logging.info('MACS peak calling started.')
        subprocess.check_call('''module load python/2.7.3; module load '''
            '''numpy/1.8.1; module load setuptools/25.1.1; module load '''
            '''MACS/2.1.0; macs2 callpeak -t %s.for_macs.bed --nomodel '''
            '''--keep-dup all --extsize 200 --shift -100 -f BED -g hs -n '''
            '''%s/%s_macs --call-summits''' % (OUTPUT_PREFIX, MACS_DIRECTORY, args.prefix), shell=True)
        logging.info('MACS peak calling ended.')
    else:
        logging.info('MACS skipped.')

    if not os.path.exists(OUTPUT_PREFIX + ".intersect.bed") or \
        args.force_overwrite_all:
        logging.info('Intersect started.')
        args.force_overwrite_all = True
        subprocess.check_call("bedtools intersect -b %s.for_macs.bed -a %s/%s_macs_peaks.narrowPeak -wa -wb >  %s.intersect.bed" % (OUTPUT_PREFIX, MACS_DIRECTORY, args.prefix, OUTPUT_PREFIX), shell=True)
        logging.info('Intersect ended.')

    if not os.path.exists(OUTPUT_PREFIX + ".counts.txt") or \
        args.force_overwrite_all:
	args.force_overwrite_all = True
        logging.info('Count matrix started.')
        subprocess.check_call('''awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $14}' %s.intersect.bed | awk '!x[$0]++' | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $4}' | awk 'BEGiN {OFS = "\t"}; {h[$0]++}; END { for(k in h) print k, h[k] }' | awk 'BEGIN {OFS="\t"}; {print $1 "_" $2 "_" $3, $4, $5}' > %s.counts.txt''' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        logging.info('Count matrix ended.')

    if not args.keep_intermediates:
    # Remove temporary files created during the pipeline.
        clean_command = ('rm %s.for_macs.bed; rm %s.intersect.bed; rm %s.clean.bed; rm %s.merge.bam*; rm %s.true.nodups.bam;' %
        (OUTPUT_PREFIX, OUTPUT_PREFIX, OUTPUT_PREFIX, OUTPUT_PREFIX, OUTPUT_PREFIX))
        subprocess.check_call(clean_command, shell=True)

    logging.info('Mergeall Complete.')
