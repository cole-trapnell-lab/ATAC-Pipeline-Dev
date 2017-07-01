#Single-cell ATAC-seq Read Processing Pipeline
The purpose of this pipeline is to get reads from bcl (directly off the illumina machines) to a peak by cell matrix. Currently, this pipeline takes care of the first step: running bcl2fastq, and cleaning up and correcting barcodes. Upcoming are additions to run map reads, run hotspot, and make an overlap matrix. I am still very actively debugging, so report issues if you have them.

###Installation:
`module load julia/latest`
`julia
Pkg.clone("https://github.com/hpliner/Levenshtein.jl.git")
Pkg.add("GZip")
Pkg.add("ArgParse")`

Currently, this must be run on the lab cluster because it recruits cluster modules.
Requires the following modules, add this to your path:
`module load coreutils/8.24`
`module load pysam/0.8.1`

###Example run
`python runall.py -R /net/shendure/vol9/seq/NEXTSEQ/160311_NS500488_0144_AH3L3MBGXY/ -O /net/trapnell/vol1/hannah/data/ -P Experiment1`

The first argument is the full path to your flowcell data. The second argument is the full path to where you want the output to go. The last argument is optional, and gives each output file an experiment prefix.

For more details on furhter arguments, run `python runall.py --help`

###Output
1. A folder called `fastq` that holds the original fastq output from bcl2fastq along with some stats from illumina
2. bcl2fastq_log.txt which holds the usual output from bcl2fastq - cluster density etc.
3. prefix.split.#.trimmed.paired.fastq.gz which has the trimmed and cleaned fastqs
