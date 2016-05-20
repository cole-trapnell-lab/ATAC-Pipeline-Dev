import argparse
import os
import subprocess
import gzip
import os.path
import sys
import logging

PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))
DEDUPLICATER = os.path.join(PIPELINE_PATH, 'sc_atac_true_dedup.py')
PICARD = os.path.join(PIPELINE_PATH, 'picard-tools-1.141/picard.jar')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to aggregate mapped bams from sequencing runs and run hotspot for scATAC-seq analysis.')
    parser.add_argument('-B','--bamlist', nargs='*', help='Paths to bams to combined', dest='bamlist', required=True)
    parser.add_argument('-O','--outdir', help='Output directory', dest='outdir', required=True)
    parser.add_argument('-P','--prefix',help='Output file prefix, otherwise default(out)', default = "out", dest='prefix', required=False)
    parser.add_argument('-C','--barcodes', help='Barcodes combinations allowed in a text file, if none provided, no filtering on possible barcodes.', required=False, dest='barcodes', default="None")
    parser.add_argument('--run_hotspot', action='store_true',
        help='Create tag file and run hotspot with default conditions')
    parser.add_argument('--no_complexity', action='store_true',
        help='Add flag if you would like to skip running picard tools EstimateLibraryComplexity')
    parser.add_argument('--force_overwrite_all', action='store_true',
        help='Force overwrite of all steps of pipeline regardless of files '
        'already present.')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)

    # Configure logger
    logging.basicConfig(filename= OUTPUT_PREFIX + '.log',format='%(asctime)s '
        '%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.info('Pipeline started.')

    # Make sorted bed file from all bams
    if not os.path.exists(OUTPUT_PREFIX + ".merge.bam") or \
        args.force_overwrite_all:
        logging.info('Merge started.')
        subprocess.call('samtools merge %s.merge.bam %s' % (OUTPUT_PREFIX, ' '.join(args.bamlist)), shell=True)
        subprocess.call('samtools index %s.merge.bam' % OUTPUT_PREFIX, shell=True)
        logging.info('Merge ended.')
        if not args.no_complexity:
            subprocess.call('java -jar %s EstimateLibraryComplexity I=%s.merge.bam O=%s.complexity_metrics.txt QUIET=true' % (PICARD, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

    else:
        print 'Bams already merged, skipping.'
        logging.info('Merge skipped.')

    if not os.path.exists(OUTPUT_PREFIX + ".true.nodups.bam"):
        logging.info('Deduplication started.')
        subprocess.call('python %s %s.merge.bam %s.true.nodups.bam' %
            (DEDUPLICATER, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
        subprocess.call('samtools view %s.true.nodups.bam | sort -u -k1,1 | '
            'cut -f9 > %s.insertsize.txt' % (OUTPUT_PREFIX, OUTPUT_PREFIX),
            shell=True)
        subprocess.call('samtools index %s.true.nodups.bam' % OUTPUT_PREFIX,
            shell=True)
        logging.info('Deduplication ended.')

    else:
        print 'Sequences already deduplicating, skipping.'
        logging.info('Deduplication skipped.')

    if not os.path.exists(OUTPUT_PREFIX + ".all.bed") or \
            args.force_overwrite_all:
        subprocess.call('bedtools bamtobed -i %s.true.nodups.bam > %s.all.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

    if not os.path.exists(OUTPUT_PREFIX + ".sort.bed") or \
        args.force_overwrite_all:
        subprocess.call('bedtools sort -i %s.all.bed > %s.sort.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

    if not os.path.exists(OUTPUT_PREFIX + ".bc.bed") or \
        args.force_overwrite_all:

        # Split cell name to only have barcode
        subprocess.call('''awk 'gsub(/(:| )+/,"\t")' %s.sort.bed > %s.bc.bed''' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

        if args.barcodes != "None":
            # Keep only barcodes that are allowed
            subprocess.call("grep -Fwf %s %s.bc.bed > %s.bc_only.bed" % (args.barcodes, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
            after_bc = "%s.bc_only.bed" % OUTPUT_PREFIX
        else:
            after_bc = "%s.bc.bed" % OUTPUT_PREFIX

        # Count cell reads
        subprocess.call("awk '{h[$4]++}; END { for(k in h) print k, h[k] }'  %s > %s.cell_read_counts.txt" % (after_bc, OUTPUT_PREFIX), shell=True)

    if args.run_hotspot and (not os.path.exists(OUTPUT_PREFIX +
        ".hotspot_tags.bed") or args.force_overwrite_all):
        subprocess.call('''awk {'if ($6 == "-") print $1 "\t" $3 "\t" ($3+1) "\t" $4; else print $1 "\t" $2 "\t" ($2+1) "\t" $4'} %s | bedToBam -i - -g %s/human.hg19.genome > %s.hotspot_tags.bam''' % (after_bc, PIPELINE_PATH, OUTPUT_PREFIX), shell=True)
