import argparse
import os
import subprocess
import gzip
import os.path
import sys

PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A program to aggregate mapped bams from sequencing runs and run hotspot for scATAC-seq analysis.')
    parser.add_argument('-B','--bamlist', nargs='*', help='Paths to bams to combined', dest='bamlist')
    parser.add_argument('-O','--outdir', help='Output directory', dest='outdir')
    parser.add_argument('-P','--prefix',help='Output file prefix, otherwise default(out)', default = "out", dest='prefix')
    parser.add_argument('-C','--barcodes', help='Barcodes combinations allowed in a text file.', default=3, dest='barcodes')
    parser.add_argument('--force_overwrite_all', action='store_true',
        help='Force overwrite of all steps of pipeline regardless of files '
        'already present.')
    args = parser.parse_args()
    
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)

    # Make sorted bed file from all bams
    if not os.path.exists(OUTPUT_PREFIX + ".all.bed") or \
        args.force_overwrite_all: 
        for bam in args.bamlist:
            subprocess.call('bedtools bamtobed -i %s >> %s.all.bed' % (bam, OUTPUT_PREFIX), shell=True)

    if not os.path.exists(OUTPUT_PREFIX + ".sort.bed") or \
        args.force_overwrite_all:
        subprocess.call('bedtools sort -i %s.all.bed > %s.sort.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
    
    if not os.path.exists(OUTPUT_PREFIX + ".bc.bed") or \
        args.force_overwrite_all:
        
        # Split cell name to only have barcode
        subprocess.call('''awk 'gsub(/(:| )+/,"\t")' %s.sort.bed > %s.bc.bed''' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

        # Keep only barcodes that are allowed
        subprocess.call("grep -Fwf %s %s.bc.bed > %s.bc_only.bed" % (args.barcodes, OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

        # Count cell reads
        subprocess.call("awk '{h[$4]++}; END { for(k in h) print k, h[k] }'  %s.bc_only.bed > %s.cell_read_counts.txt" % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

    subprocess.call('''awk {'if ($6 == "-") print $1 "\t" $3 "\t" ($3+1) "\t" $4; else print $1 "\t" $2 "\t" ($2+1) "\t" $4'} %s.bc_only.bed | bedToBam -i - -g %s/human.hg19.genome > %s.hotspot_tags.bam''' % (OUTPUT_PREFIX, PIPELINE_PATH, OUTPUT_PREFIX), shell=True)






