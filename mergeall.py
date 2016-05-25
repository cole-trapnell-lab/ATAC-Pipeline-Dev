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
        ".hotspot_tags.bam") or args.force_overwrite_all):
        subprocess.call('''awk {'if ($6 == "-") print $1 "\t" $3 "\t" ($3+1) "\t" $4; else print $1 "\t" $2 "\t" ($2+1) "\t" $4'} %s | bedToBam -i - -g %s/human.hg19.genome > %s.hotspot_tags.bam''' % (after_bc, PIPELINE_PATH, OUTPUT_PREFIX), shell=True)
        if not os.path.exists(os.path.join(args.outdir, 'hotspot_calls')):
            os.mkdir(os.path.join(args.outdir, 'hotspot_calls'))

    if args.run_hotspot and (not os.path.join(args.outdir, 'hotspot_calls/', args.prefix, '.hotspot_tags-final'):
        ".hotspot_tags-final") or args.force_overwrite_all):
        hotspot_tokens = "_TAGS_ = %s.hotspot_tags.bam\n_USE_INPUT_ = F\n_GENOME_ = hg19\n_K_ = 36\n_CHROM_FILE_ = %s/hotspot-distr/data/hg19.chromInfo.bed\n_MAPPABLE_FILE_ = %s/hotspot-distr/data/hg19.K36.mappable_only.bed.starch\n_DUPOK_ = T\n_FDRS_ = '0.01'\n_DENS_:\n_OUTDIR_ = %s\n_RANDIR_ = %s\n_OMIT_REGIONS_: %s/hotspot-distr/data/Satellite.hg19.bed\n_CHECK_ = T\n_CHKCHR_ = chrX\n_HOTSPOT_ = %s/hotspot-distr/hotspot-deploy/bin/hotspot\n_CLEAN_ = T\n_PKFIND_BIN_ = %s/hotspot-distr/hotspot-deploy/bin/wavePeaks\n_PKFIND_SMTH_LVL_ = 3\n_SEED_=101\n_THRESH_ = 2\n_WIN_MIN_ = 200\n_WIN_MAX_ = 300\n_WIN_INCR_ = 50\n_BACKGRD_WIN_ = 50000\n_MERGE_DIST_ = 150\n_MINSIZE_ = 10\n" % (OUTPUT_PREFIX, PIPELINE_PATH, PIPELINE_PATH, os.path.join(args.outdir, 'hotspot_calls'), os.path.join(args.outdir, 'hotspot_calls'), PIPELINE_PATH, PIPELINE_PATH )

        hotspot_file = open(os.path.join(args.outdir, "runall.tokens.txt"), 'w')
        hotspot_file.write(hotspot_tokens)
        hotspot_file.close()

        rh = "#! /usr/bin/env bash\nset -e -o pipefail\nscriptTokBin=%s/hotspot-distr/ScriptTokenizer/src/script-tokenizer.py\npipeDir=%s/hotspot-distr/pipeline-scripts\ntokenFile=%s\nscripts='$pipeDir/run_badspot\n$pipeDir/run_make_lib\n$pipeDir/run_wavelet_peak_finding\n$pipeDir/run_10kb_counts\n$pipeDir/run_generate_random_lib\n$pipeDir/run_pass1_hotspot\n$pipeDir/run_pass1_merge_and_thresh_hotspots\n$pipeDir/run_pass2_hotspot\n$pipeDir/run_rescore_hotspot_passes\n$pipeDir/run_spot\n$pipeDir/run_thresh_hot.R\n$pipeDir/run_both-passes_merge_and_thresh_hotspots\n$pipeDir/run_add_peaks_per_hotspot\n$pipeDir/run_final'\n\n$scriptTokBin --clobber --output-dir=args.outdir $tokenFile $scripts\n\nfor script in $scripts\ndo\n\t./$(basename $script).tok\ndone" %(PIPELINE_PATH, PIPELINE_PATH, os.path.join(args.outdir, "runall.tokens.txt"))

        hotspot_run = open(os.path.join(args.outdir, "runhotspot"), 'w')
        hotspot_run.write(rh)
        hotspot_run.close()

        subprocess.call(os.path.join(args.outdir, 'runhotspot'), shell=True)

    if args.run_hotspot and not os.path.exists(OUTPUT_PREFIX + ".intersect_hotspot_counts.bed":
        subprocess.call("bedtools merge -d 150 -c 1 -o count -i %s > %s" % (OUTPUT_PREFIX + ".hotspot_tags-final/" + args.prefix + ".hotspot_tags.fdr0.01.pks.bed", OUTPUT_PREFIX + ".hotspot_merge_pks.bed"), shell=True)

        subprocess.call("bedtools intersect -b %s -a %s -wa -wb > %s" % (after_bc, OUTPUT_PREFIX + ".hotspot_merge_pks.bed", OUTPUT_PREFIX + ".hotspot_intersect.bed"), shell=True)

        subprocess.call('''awk 'BEGIN {OFS="\t"}; {print $1, $2, $3, $8}' %s | awk '{h[$0]++}; END { for(k in h) print k, h[k] }' > %s''' % (OUTPUT_PREFIX + ".hotspot_intersect.bed", OUTPUT_PREFIX + ".intersect_hotspot_counts.bed"), shell=True)
