#Python pipeline to convert BCL files to qseq,
#convert qseq to fastq, fix errors in barcodes
#and trim adapters from fastqs.

#Current issues:
# - Depends on perl scripts (yuck!) in Darren's
#   home directory to run properly.

# - Does not currently allow you to modify
#   default parameters of the perl scripts.

# - Hardcodes trimmomatic as located in Darren's
#   home directory and requires cd'ing there to
#   run the program properly.

import argparse
import os
import subprocess

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

parser = argparse.ArgumentParser(description='A program to convert MiSeq BCL files to fastq files for scATAC-seq analysis.')
parser.add_argument('-R','--rundir', help='Run directory containing BCL files',dest='rundir')
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir')
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix')
args = parser.parse_args()

moduler = 'module load modules modules-init modules-gs boost/1.35.0 fftw/3.2.2 OLB/1.9.3'
submitter(moduler)

setter = 'setupBclToQseq.py -b ' + args.rundir + 'Data/Intensities/BaseCalls/ -o ' + args.rundir + 'Data/Intensities/BaseCalls/qseq/ -P .locs --overwrite'
submitter(setter)

os.chdir(args.rundir + 'Data/Intensities/BaseCalls/qseq/')

maker = 'make -j 8'
submitter(maker)

os.mkdir(args.outdir)

fastqer = 'perl /net/shendure/vol1/home/cusanovi/scatac/NCP_MiSeq_qseq_to_fastq.pl ' + args.outdir + args.prefix
submitter(fastqer)

splitter = 'perl /net/shendure/vol1/home/cusanovi/scatac/NCP_fastq_10bpbarcode_split.pl -1 ' + args.outdir + args.prefix + '.1.fq.gz -2 ' + args.outdir + args.prefix + '.2.fq.gz -O ' + args.outdir + args.prefix
submitter(splitter)

os.chdir('/net/shendure/vol1/home/cusanovi/')

trimmer = 'java -Xmx1G -jar /net/shendure/vol1/home/cusanovi/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE '+ args.outdir + args.prefix + '.split.1.fq.gz ' + args.outdir + args.prefix + '.split.2.fq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz ILLUMINACLIP:bin/Trimmomatic-0.32/adapters/NexteraPE-PE.fa:2:30:10:4:true MINLEN:36'
submitter(trimmer)

cleaner = 'rm ' + args.outdir + args.prefix + '.split.1.fq.gz; rm ' + args.outdir + args.prefix + '.split.2.fq.gz; rm ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz; rm ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz;'
submitter(cleaner)
