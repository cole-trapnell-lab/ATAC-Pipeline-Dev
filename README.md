# sci-ATAC-seq Read Processing Pipeline
The purpose of this pipeline is to get reads from bcl (directly off the illumina machines) to a peak by cell matrix. Currently, this pipeline takes care of the first step: running bcl2fastq, cleaning up and correcting barcodes and mapping.  I am still very actively debugging, so report issues if you have them.

### Installation:
~~~~
module load julia/latest
julia
Pkg.clone("https://github.com/hpliner/Levenshtein.jl.git")
Pkg.add("GZip")
Pkg.add("ArgParse")
~~~~

Currently, this must be run on the lab cluster because it recruits cluster modules.

### Basic usage:

Make a script called `runcall.sh` with the following contents:
~~~~ 
#$ -pe serial 10
#$ -l mfree=10G

module load julia/latest
module load pysam/0.8.1
module load coreutils/8.24
path/to/runall.py -R [Flowcell run directory] -O [Path to output folder] -P [Prefix you want on output files] -p 10 
~~~~ 

To execute, run `qsub runcall.sh`

The first argument of runall.py is the full path to your flowcell data. The second argument is the full path to where you want the output to go. The third argument is optional, and gives each output file an experiment prefix.

(example: `runall.py -R /net/shendure/vol9/seq/170607_NS500488_0394_AHVHL5BGX2 -O /net/trapnell/vol1/hannah/Projects/tests/ -P hifT -p 10`)

For more details on further arguments, run `python runall.py --help`

### Output
1. A folder called `fastq` that holds the original fastq output from bcl2fastq along with some stats from illumina and from barcode correction.
2. bcl2fastq_log.txt which holds the usual output from bcl2fastq - cluster density etc.
