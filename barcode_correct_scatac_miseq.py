#!/usr/bin/python
##
##	barcode_correct_scatac.py
##
##  Barcode correction for scATAC-seq
##
##  Begun                 2016-05-02 HAP
##  Last modified         2016-05-02 HAP
##

import argparse
import os
import subprocess
import difflib
import gzip
import datetime
import runall
import os.path
import multiprocessing
import glob
from operator import itemgetter
import operator
from itertools import imap, repeat

# Standard nextera barcodes
NEX_I7 = ["ATTACTCG", "TCCGGAGA", "CGCTCATT", "GAGATTCC", "ATTCAGAA",
    "GAATTCGT", "CTGAAGCT", "TAATGCGC", "CGGCTATG", "TCCGCGAA", "TCTCGCGC",
    "AGCGATAG"]

PCR_I7 = ["TCGGATTCGG" ,"GCGGCTGCGG", "AGATTACGTT", "CTAACTAGGT", "CATAGCGACC",
    "CCGCTAAGAG", "ATGGAACGAA", "GCGTTCCGTT", "GGTTATCGAA", "GCATCGTATG",
    "AATACGATAA", "TTCCGTCGAC", "TCCGGCTTAT", "ACCAGGCGCA", "AGAGGAGAAT",
    "GTACTCCTAT", "GCTAACGGAT", "AGTTGAATCA", "TGATTAGGTA", "TCGTAGCATC",
    "TCTTGAGGTT", "AGGTCAGCTT", "TATTAGACTT", "CTCAATTAGT", "TCGCCGCCGG",
    "CCGTATGATT", "AACGCGCAGA", "CTCGTCGTAG", "CTAATTGCGA", "CGCGGCCATA",
    "AATATTACTT", "ATTGGCAGAT", "ATGGCGCCTG", "ATAAGGACTC", "TAGTAAGCCG",
    "ATTATGCAAG", "TTGGCAAGCC", "TTGATTGGCG", "GCATATGAGC", "GAACTCGACT",
    "CTAGCCAGCC", "TGCGACCTCT", "ATTCTTAGCT", "TTGATACGAT", "TATAATAGTT",
    "TTGCCGTAGG", "AGACCATATC", "TTGGTAAGGA", "CAGCTAGCGG", "CTAAGCCTTG",
    "CGTTACCGCT", "GACTGGACCA", "GCAAGACCGT", "TCAATCTCCT", "ATACCTCGAC",
    "TAGAGGCGTT", "TAGGTAACTT", "TTCGAATATT", "TGGACGACTA", "GTAGGCTGCA",
    "GTAGGATAAG", "CGTCGAGCGC", "ACTATTCATT", "TTGCTTAGAT", "CGAATGGAGC",
    "CTATATAGCC", "CTACTAATAA", "TGGTTGCCGT", "TCCTCTGCCG", "GATTCTTGAA",
    "GTAGCAGCTA", "CCTCAGCTCC", "AAGTAGCTCA", "TATTGCTGGA", "CCAGATACGG",
    "AACGAATTCG", "CGCTTATCGT", "AAGTACGCGA", "GATCTTCGCA", "TCTTAGCCTG",
    "TTATTGAGGC", "TTGCGAGCAT", "GCTTGAAGAG", "AGTCCGCTGC", "TAAGTCCTGA",
    "AGTTCTCATG", "CAGACTAAGG", "TCTATCGCTG", "GCGCTATGGT", "CATTATTATT",
    "AGCCGTAGTT", "TGATATTGCG", "ACGGCGTTAA", "GGCTTACTCC", "GCGCGTTCAT",
    "GAGCGCGATG"]

PCR_I5 = ["CTCCATCGAG", "TTGGTAGTCG", "GGCCGTCAAC", "CCTAGACGAG", "TCGTTAGAGC",
    "CGTTCTATCA", "CGGAATCTAA", "ATGACTGATC", "TCAATATCGA", "GTAGACCTGG",
    "TTATGACCAA", "TTGGTCCGTT", "GGTACGTTAA", "CAATGAGTCC", "GATGCAGTTC",
    "CCATCGTTCC", "TTGAGAGAGT", "ACTGAGCGAC", "TGAGGAATCA", "CCTCCGACGG",
    "CATTGACGCT", "TCGTCCTTCG", "TGATACTCAA", "TTCTACCTCA", "TCGTCGGAAC",
    "ATCGAGATGA", "TAGACTAGTC", "GTCGAAGCAG", "AGGCGCTAGG", "AGATGCAACT",
    "AAGCCTACGA", "GTAGGCAATT", "GGAGGCGGCG", "CCAGTACTTG", "GGTCTCGCCG",
    "GGCGGAGGTC", "TAGTTCTAGA", "TTGGAGTTAG", "AGATCTTGGT", "GTAATGATCG",
    "CAGAGAGGTC", "TTAATTAGCC", "CTCTAACTCG", "TACGATCATC", "AGGCGAGAGC",
    "TCAAGATAGT", "TAATTGACCT", "CAGCCGGCTT", "AGAACCGGAG", "GAGATGCATG",
    "GATTACCGGA", "TCGTAACGGT", "TGGCGACGGA", "AGTCATAGCC", "GTCAAGTCCA",
    "ATTCGGAAGT", "GTCGGTAGTT", "AGGACGGACG", "CTCCTGGACC", "TAGCCTCGTT",
    "GGTTGAACGT", "AGGTCCTCGT", "GGAAGTTATA", "TGGTAATCCT", "AAGCTAGGTT",
    "TCCGCGGACT", "TGCGGATAGT", "TGGCAGCTCG", "TGCTACGGTC", "GCGCAATGAC",
    "CTTAATCTTG", "GGAGTTGCGT", "ACTCGTATCA", "GGTAATAATG", "TCCTTATAGA",
    "CCGACTCCAA", "GCCAAGCTTG", "CATATCCTAT", "ACCTACGCCA", "GGAATTCAGT",
    "TGGCGTAGAA", "ATTGCGGCCA", "TTCAGCTTGG", "CCATCTGGCA", "CTTATAAGTT",
    "GATTAGATGA", "TATAGGATCT", "AGCTTATAGG", "GTCTGCAATC", "CGCCTCTTAT",
    "GTTGGATCTT", "GCGATTGCAG", "TGCCAGTTGC", "CTTAGGTATC", "GAGACCTACC",
    "ATTGACCGAG"]

NEX_I5 = ["TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA", "AGGCGAAG",
    "TAATCTTA", "CAGGACGT", "GTACTGAC"]

def hamming(str1, str2):
    return sum(imap(operator.ne, str1, str2))

def split_files((i, R, file_list, fastqpath)):
	subprocess.call("gunzip -c %s | split -l 30000000 -d -a 4 - %s "
        "--filter='gzip > %s'" % (os.path.join(fastqpath, file_list[i]),
        'tempR' + str(R) + str(i) + '.fq', fastqpath + '/$FILE.gz'),
        shell=True)

def closest_match(str1, strlist):
	all_h = [hamming(str1, s2) for s2 in strlist]
	return strlist[min(enumerate(all_h), key=itemgetter(1))[0]]

def reverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

def clean_and_correct((ifile, fastqpath)):
    fileR1 = ifile
    fileR2 = ifile.replace('R1','R2', 1)
    kept = 0
    with gzip.open(os.path.join(fastqpath, fileR1), 'rb') as f:
        with gzip.open(os.path.join(fastqpath, fileR2), 'rb') as r:
            with gzip.open(fileR1 + '.out.fq.gz', 'ab') as o:
                with gzip.open(fileR2 + '.out.fq.gz', 'ab') as g:
                    for tag_line in f:
                        tag_line = tag_line.strip().split()[1].split(':')\
                        [3].replace('+','')
			if len(tag_line) != 36:
                            tag_line2 = next(r)
                            read_line2 = next(r)
                            plus_line2 = next(r)
                            qual_line2 = next(r)
                            continue
                        read_line = next(f)
                        plus_line = next(f)
                        qual_line = next(f)
                        b1 = tag_line[0:8]
                        b2 = tag_line[8:18]
                        b3 = tag_line[18:26]
                        b4 = tag_line[26:36]
                        if b1 in NEX_I7:
                            b1_cor = [b1]
                        else:
                            b1_cor = difflib.get_close_matches(b1, NEX_I7)
                        if b2 in PCR_I7:
                            b2_cor = [b2]
                        else:
                            b2_cor = difflib.get_close_matches(b2, PCR_I7, 1)
                        if b3 in PCR_I5:
                            b3_cor = [b3]
                        else:
                            b3_cor = difflib.get_close_matches(b3, PCR_I5, 1)
                        if b4 in NEX_I5:
                            b4_cor = [b4]
                        else:
                            b4_cor = difflib.get_close_matches(b4, NEX_I5, 1)
                        cor_barcode = ''.join(b1_cor + b2_cor + b3_cor +
                            b4_cor)
                        tag_new = ''.join(b1 + b2 + b3 + b4)
                        edit_dist = hamming(cor_barcode, tag_new)
                        tag_line2 = next(r)
                        read_line2 = next(r)
                        plus_line2 = next(r)
                        qual_line2 = next(r)
                        if edit_dist <= 3:
                            kept += 1
                            content = ('@' + cor_barcode + ':' + str(kept) + '#'
                                + str(edit_dist) + '/1' + '\n' + read_line +
                                plus_line + qual_line)
                            content2 = ('@' + cor_barcode + ':' + str(kept) +
                                '#' + str(edit_dist) + '/1' + '\n' + read_line2
                                + plus_line2 + qual_line2)
                            o.write(content)
                            g.write(content2)
    return kept


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A program to correct '
        'barcodes and report edit distance in scATAC-seq analysis.')
	parser.add_argument('-n','--numthreads', help='Number of threads '
        'available', dest='numthreads', required=True)
	parser.add_argument('-F','--fastqpath', help='Path for fastqs path',
        dest='fastqpath', required=True)
	parser.add_argument('-o','--outpref', help='Output file prefix, otherwise '
        'default(out)', default = "out", dest='outpref')
	parser.add_argument('-E','--maxedit', help='Maximum allowed edit distance '
        '(default = 3)', default=3, dest='maxedit')
	args = parser.parse_args()


	# Open barcode correction log
	log = open(args.outpref + '.barcode_correct.log', 'w')
	log_mes = ('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) +
        ' Starting processing\n')
	log.write(str(log_mes))
	# Open output files
	output1 = args.outpref + 'split.1.fq.gz'
	output2 = args.outpref + 'split.2.fq.gz'

	# Collect fastqs
	fastq_files = [f for f in os.listdir(args.fastqpath) if \
        os.path.isfile(os.path.join(args.fastqpath, f))]
	R1_files = [f for f in fastq_files if 'R1' in f]
	R2_files = [f for f in fastq_files if 'R2' in f]
	log_mes = ('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) +
        ' Starting file split\n')
	log.write(str(log_mes))

	file_names1 = []

	pool = multiprocessing.Pool(processes=int(args.numthreads))
	if not os.path.exists(args.fastqpath + '/tempR10.fq0000.gz'):
		print 'starting file split'
		pool.map(split_files, zip(range(len(R1_files)), repeat(1),
            repeat(R1_files), repeat(args.fastqpath)))
		print 'done file split 1'
	else:
                print 'split1 already exists, moving on'

	file_names1 = glob.glob(args.fastqpath + '/tempR1*[!q].gz')

	if not os.path.exists(args.fastqpath + '/tempR20.fq0000.gz'):
	   	pool.map(split_files, zip(range(len(R2_files)), repeat(2),
            repeat(R2_files), repeat(args.fastqpath)))
		print 'done file split 2'
	else:
		print 'split2 already exists, moving on'

	if not os.path.exists(args.fastqpath + '/tempR10.fq0000.gz.out.fq.gz'):
		kept_list = pool.map(clean_and_correct,
            zip(file_names1,repeat(args.fastqpath)))
		log_mes = "Sequences kept: " + str(sum(kept_list)) + "\n"
		log.write(log_mes)
	else:
		print 'barcode corrected files already exist, moving on'

	outfiles1 = [f + '.out.fq.gz' for f in file_names1]

	for f in outfiles1:
		subprocess.call("cat %s >> %s" % (f, output1), shell=True)
		f2 = f.replace('R1', 'R2', 1)
		subprocess.call("cat %s >> %s" % (f2, output2), shell=True)

	log_mes = ('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) +
        ' Done\n')
	log.write(str(log_mes))
	log.close()
