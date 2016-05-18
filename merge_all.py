import argparse
import os
import subprocess
import gzip
import os.path
import sys


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A program to aggregate mapped bams from sequencing runs and run hotspot for scATAC-seq analysis.')
	parser.add_argument('-B','--bamlist', help='Paths to bams to combined', dest='bamlist')
	parser.add_argument('-O','--outdir', help='Output directory', dest='outdir')
	parser.add_argument('-P','--prefix',help='Output file prefix, otherwise default(out)', default = "out", dest='prefix')
	parser.add_argument('-C','--barcodes', help='Barcodes combinations allowed in a text file.', default=3, dest='barcodes')

    if not os.path.exists(outdir):
        os.mkdir(outdir)
	OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)

    # Make sorted bed file from all bams
	for bam in bamlist:
		subprocess.call('bedtools bamtobed -i %s >> %s.all.bed' % (bam, OUTPUT_PREFIX), shell=True)

	subprocess.call('bedtools sort -i %s.all.bed > %s.sort.bed' % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)

    # Split cell name to only have barcode
	subprocess.call('''awk 'gsub(/(:| )+/,"\t")''' %s.sort.bed > %s.bc.bed

    # Keep only barcodes that are allowed
    subprocess.call("grep -Fwf %s %s.bc.bed > %s.bc_only.bed" % (barcodes, OUTPUT_PREFIX, OUTPUT_PREFIX)

    #print "Reads kept after matching barcodes: %s, %s%" %()
    # Count cell reads
    subprocess.call("awk '{h[$4]++}; END { for(k in h) print k, h[k] }'  %s.bc_only.bed > %s.cell_read_counts.txt" % (OUTPUT_PREFIX, OUTPUT_PREFIX), shell=True)
