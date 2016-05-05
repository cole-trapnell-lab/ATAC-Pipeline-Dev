#Python pipeline to convert NextSeq BCL files to fastq,
#fix errors in barcodes and trim adapters from fastqs.

#Cannibalized from Darren with a bunch of stuff stolen from Andrew and compiled/maintained by Hannah


import argparse
import os
import subprocess
import gzip
import os.path
import sys



# Construct paths to pipeline scripts as constants
PIPELINE_PATH = os.path.dirname(os.path.realpath(__file__))


BARCODE_CORRECTER = os.path.join(PIPELINE_PATH, 'barcode_correct_scatac.py')
TRIMMOMATIC = os.path.join(PIPELINE_PATH, 'Trimmomatic-0.36/trimmomatic-0.36.jar')

COMMANDS_TO_ARRAY_JOB_SCRIPT = os.path.join(PIPELINE_PATH, 'commands2arrayjob.sh')

SGE_LOGS_DIRECTORY = 'sge_logs'
SGE_CONFIGURATION = '-b y -shell y -cwd -j y -o %s -e %s' % (SGE_LOGS_DIRECTORY, SGE_LOGS_DIRECTORY)

def initialize_directories(pipeline_root_directory):
	fastqs = os.path.join(args.outdir, 'fastqs')
	sge_logs = os.path.join(pipeline_root_directory, SGE_LOGS_DIRECTORY)

	
	for directory in [sge_logs,fastqs]:
		if not os.path.exists(directory):
			os.mkdir(directory)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A program to convert NextSeq BCL files to cleaned and corrected fastq files for scATAC-seq analysis.')
	parser.add_argument('-R','--rundir', help='Run directory containing BCL files', dest='rundir')
	parser.add_argument('-O','--outdir', help='Output directory', dest='outdir')
	parser.add_argument('-P','--prefix',help='Output file prefix, otherwise default(out)', default = "out", dest='prefix')
	parser.add_argument('-E','--maxedit', help='Maximum allowed edit distance (default = 3)', default=3, dest='maxedit')
	parser.add_argument('--force_overwrite_all', action='store_true', help='Force overwrite of all steps of pipeline regardless of files already present.')
	parser.add_argument('--force_overwrite_bcl2fastq', action='store_true', help='Force overwrite bcl2fastq regardless of files already present.')
	parser.add_argument('--force_overwrite_cleanfastq', action='store_true', help='Force overwrite fastq cleaning regardless of files already present.')
	parser.add_argument('--force_overwrite_barcodecorrect', action='store_true', help='Force overwrite barcode correction regardless of files already present.')
	parser.add_argument('--force_overwrite_trimming', action='store_true', help='Force overwrite trimming regardless of files already present.')
	args = parser.parse_args()

	# Initialize directories required for the pipeline
	print 'Initializing directories for pipeline...'
	initialize_directories(args.outdir)

	OUTPUT_PREFIX = os.path.join(args.outdir, args.prefix)
	FASTQ_DIRECTORY = os.path.join(args.outdir, 'fastqs')

	bar_out1 = OUTPUT_PREFIX + 'split.1.fq.gz'
	bar_out2 = OUTPUT_PREFIX + 'split.2.fq.gz'

	# Submit bcl2fastq only if no existing results or if user wants to overwrite
	if len(os.listdir(FASTQ_DIRECTORY)) == 0 or args.force_overwrite_all or args.force_overwrite_bcl2fastq:
		print 'Starting bcl2fastq...'
		bcl2fastq_command = 'module load modules modules-init modules-gs bcl2fastq/2.16 fastqc/0.10.1; bcl2fastq --runfolder-dir %s -o %s --ignore-missing-filter' % (args.rundir, FASTQ_DIRECTORY)  
		f = open(os.path.join(args.outdir, "bcl2fastq_log.txt"), 'w')
		subprocess.call(bcl2fastq_command, shell=True, stdout = f, stderr=f)
	else:
		print 'fastq directory is not empty, skipping bcl2fastq. Specify --force_overwrite_all or --force_overwrite_bcl2fastq to redo.'

	print "Cleaning and fixing barcodes..."
	print FASTQ_DIRECTORY
	if not os.path.exists(bar_out1) or args.force_overwrite_all or args.force_overwrite_barcodecorrect:
		subprocess.call('python %s -F %s -o %s -E %s -n 10' % (BARCODE_CORRECTER, FASTQ_DIRECTORY, OUTPUT_PREFIX, args.maxedit), shell=True)
	else:
		print 'Barcodes already fixed, skipping fix barcodes. Specify --force_overwrite_all or --force_overwrite_barcodecorrect to redo.'
	
	print "Trimming adapters..."

	trimmer_out1 = OUTPUT_PREFIX + '.split.1.trimmed.paired.fastq.gz'
	trimmer_out2 = OUTPUT_PREFIX + '.split.2.trimmed.paired.fastq.gz'
	trimmer_un_out1 = OUTPUT_PREFIX + '.split.1.trimmed.unpaired.fastq.gz'
	trimmer_un_out2 = OUTPUT_PREFIX + '.split.2.trimmed.unpaired.fastq.gz'
	
	if not os.path.exists(trimmer_out1) or args.force_overwrite_all or args.force_overwrite_trimming:
		trimmer_command = 'java -Xmx1G -jar %s PE %s %s %s %s %s %s ILLUMINACLIP:%sTrimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:20' % \
			(TRIMMOMATIC, BAR_OUT1, BAR_OUT2, trimmer_out1, trimmer_un_out1, trimmer_out2, trimmer_un_out2, PIPELINE_PATH)

		subprocess.call(trimmer_command, shell=True)
	else:
		print 'Sequences already trimmed, skipping trimming. Specify --force_overwrite_all or --force_overwrite_trimming to redo.'

	print "Cleaning up..."

	clean_command = 'rm %s; rm %s; rm %s; rm %s; ' % (BAR_OUT1, BAR_OUT2, trimmer_un_out1, trimmer_un_out2)
	subprocess.call(clean_command, shell=True)




