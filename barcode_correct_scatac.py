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



def ed(string1, string2):
	edist = 0
	for i in range(len(string1)):
		if string1[i] != string2[i]:
			edist += 1
	return(edist)

def reverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])



if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='A program to correct barcodes and report edit distance in scATAC-seq analysis.')
	#parser.add_argument('-F','--forwardin', help='First fastq input file', dest='forwardin', required=True)
	#parser.add_argument('-R','--reversein', help='Second fastq input file', dest='reversein', required=True)
	#parser.add_argument('-O','--outdir', help='Output directory', dest='outdir', required=True)
	#parser.add_argument('-o','--outfile', help='Output file prefix, otherwise default(out)', default = "out", dest='outfile')
	parser.add_argument('-E','--maxedit', help='Maximum allowed edit distance (default = 3)', default=3, dest='maxedit')
	args = parser.parse_args()

	nex_i7 = ["ATTACTCG","TCCGGAGA","CGCTCATT","GAGATTCC","ATTCAGAA","GAATTCGT","CTGAAGCT","TAATGCGC","CGGCTATG","TCCGCGAA","TCTCGCGC","AGCGATAG"]

	pcr_i7 = ["TCGGATTCGG","GCGGCTGCGG","AGATTACGTT","CTAACTAGGT","CATAGCGACC","CCGCTAAGAG","ATGGAACGAA","GCGTTCCGTT","GGTTATCGAA","GCATCGTATG","AATACGATAA","TTCCGTCGAC","TCCGGCTTAT","ACCAGGCGCA","AGAGGAGAAT","GTACTCCTAT","GCTAACGGAT","AGTTGAATCA","TGATTAGGTA","TCGTAGCATC","TCTTGAGGTT","AGGTCAGCTT","TATTAGACTT","CTCAATTAGT","TCGCCGCCGG","CCGTATGATT","AACGCGCAGA","CTCGTCGTAG","CTAATTGCGA","CGCGGCCATA","AATATTACTT","ATTGGCAGAT","ATGGCGCCTG","ATAAGGACTC","TAGTAAGCCG","ATTATGCAAG","TTGGCAAGCC","TTGATTGGCG","GCATATGAGC","GAACTCGACT","CTAGCCAGCC","TGCGACCTCT","ATTCTTAGCT","TTGATACGAT","TATAATAGTT","TTGCCGTAGG","AGACCATATC","TTGGTAAGGA","CAGCTAGCGG","CTAAGCCTTG","CGTTACCGCT","GACTGGACCA","GCAAGACCGT","TCAATCTCCT","ATACCTCGAC","TAGAGGCGTT","TAGGTAACTT","TTCGAATATT","TGGACGACTA","GTAGGCTGCA","GTAGGATAAG","CGTCGAGCGC","ACTATTCATT","TTGCTTAGAT","CGAATGGAGC","CTATATAGCC","CTACTAATAA","TGGTTGCCGT","TCCTCTGCCG","GATTCTTGAA","GTAGCAGCTA","CCTCAGCTCC","AAGTAGCTCA","TATTGCTGGA","CCAGATACGG","AACGAATTCG","CGCTTATCGT","AAGTACGCGA","GATCTTCGCA","TCTTAGCCTG","TTATTGAGGC","TTGCGAGCAT","GCTTGAAGAG","AGTCCGCTGC","TAAGTCCTGA","AGTTCTCATG","CAGACTAAGG","TCTATCGCTG","GCGCTATGGT","CATTATTATT","AGCCGTAGTT","TGATATTGCG","ACGGCGTTAA","GGCTTACTCC","GCGCGTTCAT","GAGCGCGATG"]

	pcr_i5 = ["CTCCATCGAG","TTGGTAGTCG","GGCCGTCAAC","CCTAGACGAG","TCGTTAGAGC","CGTTCTATCA","CGGAATCTAA","ATGACTGATC","TCAATATCGA","GTAGACCTGG","TTATGACCAA","TTGGTCCGTT","GGTACGTTAA","CAATGAGTCC","GATGCAGTTC","CCATCGTTCC","TTGAGAGAGT","ACTGAGCGAC","TGAGGAATCA","CCTCCGACGG","CATTGACGCT","TCGTCCTTCG","TGATACTCAA","TTCTACCTCA","TCGTCGGAAC","ATCGAGATGA","TAGACTAGTC","GTCGAAGCAG","AGGCGCTAGG","AGATGCAACT","AAGCCTACGA","GTAGGCAATT","GGAGGCGGCG","CCAGTACTTG","GGTCTCGCCG","GGCGGAGGTC","TAGTTCTAGA","TTGGAGTTAG","AGATCTTGGT","GTAATGATCG","CAGAGAGGTC","TTAATTAGCC","CTCTAACTCG","TACGATCATC","AGGCGAGAGC","TCAAGATAGT","TAATTGACCT","CAGCCGGCTT","AGAACCGGAG","GAGATGCATG","GATTACCGGA","TCGTAACGGT","TGGCGACGGA","AGTCATAGCC","GTCAAGTCCA","ATTCGGAAGT","GTCGGTAGTT","AGGACGGACG","CTCCTGGACC","TAGCCTCGTT","GGTTGAACGT","AGGTCCTCGT","GGAAGTTATA","TGGTAATCCT","AAGCTAGGTT","TCCGCGGACT","TGCGGATAGT","TGGCAGCTCG","TGCTACGGTC","GCGCAATGAC","CTTAATCTTG","GGAGTTGCGT","ACTCGTATCA","GGTAATAATG","TCCTTATAGA","CCGACTCCAA","GCCAAGCTTG","CATATCCTAT","ACCTACGCCA","GGAATTCAGT","TGGCGTAGAA","ATTGCGGCCA","TTCAGCTTGG","CCATCTGGCA","CTTATAAGTT","GATTAGATGA","TATAGGATCT","AGCTTATAGG","GTCTGCAATC","CGCCTCTTAT","GTTGGATCTT","GCGATTGCAG","TGCCAGTTGC","CTTAGGTATC","GAGACCTACC","ATTGACCGAG"]

	nex_i5 = ["TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA","AGGCGAAG","TAATCTTA","CAGGACGT","GTACTGAC"]


	log = open(os.path.join(runall.OUTPUT_PATH, 'log.txt', 'a'))
	log_mes = print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())) + " Starting processing\n"
	log.write(log_mes)
	
	output1 = os.path.join(runall.OUTPUT_PATH, runall.OUTPUT_PREFIX, 'split.1.fq.gz')
	output2 = os.path.join(runall.OUTPUT_PATH, runall.OUTPUT_PREFIX, 'split.2.fq.gz')

	with gzip.open(runall.BAR_OUT1, 'wb') as o:
		o.write('')
	with gzip.open(runall.BAR_OUT2, 'wb') as g:
		g.write('')


	count = 0
	kept = 0
	with gzip.open(runall.CLEAN_R1, 'rb') as f:
		with gzip.open(runall.CLEAN_R2, 'rb') as r:
			for tag_line in f:
				tag_line = tag_line.strip('@0123456789:\n')
				count += 1
				read_line = next(f)
				plus_line = next(f)
				qual_line = next(f)
				b1 = tag_line[0:8]
				b2 = tag_line[8:18]
				b4 = reverseComplement(tag_line[18:26])
				b3 = reverseComplement(tag_line[26:36])
				b1_cor = difflib.get_close_matches(b1, nex_i7, 1)
				b2_cor = difflib.get_close_matches(b2, pcr_i7, 1)
				b3_cor = difflib.get_close_matches(b3, pcr_i5, 1)
				b4_cor = difflib.get_close_matches(b4, nex_i5, 1)
				cor_barcode = ''.join(b1_cor + b2_cor + b3_cor + b4_cor)
				tag_new = ''.join(b1 + b2 + b3 + b4)
				edit_dist = ed(cor_barcode, tag_new)
				tag_line2 = next(r)
				read_line2 = next(r)
				plus_line2 = next(r)
				qual_line2 = next(r)
				if edit_dist <= args.maxedit:
					kept += 1
					content = '@' + cor_barcode + ':' + str(count) + '#' + str(edit_dist) + '/1' + '\n' + read_line + plus_line + qual_line
					content2 = '@' + cor_barcode + ':' + str(count) + '#' + str(edit_dist) + '/1' + '\n' + read_line2 + plus_line2 + qual_line2
					with gzip.open(runall.BAR_OUT1, 'ab') as o:
						o.write(content)
					with gzip.open(runall.BAR_OUT2, 'ab') as g:
						g.write(content2)


	log_mes = "Sequences processed: " + str(count)  + "\nSequences kept: " + str(kept) + " " + str(kept/float(count)) + " %" 
	log.write(log_mes)

	log_mes = print('{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())) + " Done\n"
	log.write(log_mes)
	log.close()




