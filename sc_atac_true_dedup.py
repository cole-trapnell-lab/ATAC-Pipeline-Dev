import pysam   #/0.8.1
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]

readsin = pysam.AlignmentFile(inbam,"rb")
readsout = pysam.AlignmentFile(outbam,"wb",template=readsin)
refs = readsin.references
for refchrom in refs:
	if refchrom == 'chrM':
		continue
	readdic = {}
	print "Deduplicating " + refchrom + "..."
	for read in readsin.fetch(refchrom):
		readname = read.qname.split(':')[0]
		if 'CTF' in readname or 'AMBIG' in readname:
			continue
		readpos = str(read.pos + 1)
		matepos = str(read.mpos)
		try:
			readdic[readname + readpos + matepos]
		except KeyError:
			readdic[readname + readpos + matepos] = read
			readsout.write(read)

