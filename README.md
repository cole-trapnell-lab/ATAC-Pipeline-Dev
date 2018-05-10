# sci-ATAC-seq Read Processing Pipeline
The purpose of this pipeline is to get reads from bcl (directly off the illumina machines) to a peak by cell matrix. Currently, this pipeline takes care of the first step: running bcl2fastq, cleaning up and correcting barcodes and mapping.  I am still very actively debugging, so report issues if you have them.

## Installation:
~~~~
module load julia/latest
julia
Pkg.clone("https://github.com/hpliner/Levenshtein.jl.git")
Pkg.add("GZip")
Pkg.add("ArgParse")
~~~~

Currently, this must be run on the lab cluster because it recruits cluster modules.

## Part 1 - Flowcell to QC-ed Sorted BAM
### Basic usage:

Make a script called `runcall.sh` with the following contents: 
NOTE: These defaults assume that you are mapping to hg19 and that you ran the run on a NextSeq. If not, use `--help` to see alternative flags.
~~~~ 
#$ -pe serial 10
#$ -l mfree=10G

module load julia/latest
module load pysam/0.8.1
module load coreutils/8.24
path/to/runall.py -R [Flowcell run directory] -O [Path to output folder] 
  -P [Prefix you want on output files] -p 10 --keep_intermediates
~~~~ 

To execute, run `qsub runcall.sh`

The first argument of runall.py is the full path to your flowcell data. The second argument is the full path to where you want the output to go. The third argument is optional, and gives each output file an experiment prefix. --keep_intermediates tells the pipeline to keep all the intermediate files which lets you restart where you were when something goes wrong. A good idea to have it if there's a chance something will get messed up.

(example: `runall.py -R /net/shendure/vol9/seq/170607_NS500488_0394_AHVHL5BGX2 -O /net/trapnell/vol1/hannah/Projects/tests/ -P hifT -p 10`)

For more details on further arguments, run `python runall.py --help`

### Output
1. A folder called `fastq` that holds:
  a. The original fastq output from bcl2fastq along with some stats from illumina.
2. bcl2fastq_log.txt which holds the usual output from bcl2fastq - cluster density etc.
3. Prefix.log, that has times when processes started and ended.
4. Prefix.split.q10.sort.bam/bai - an indexed and sorted bam with all mapped reads with quality above 10.
5. A folder called `qc_info` that contains some QC stats and plots
some other stuff

## Part 2 - QC-ed Sorted BAM(s) to sparse matrix
### Basic usage:

Make a script called `mergecall.sh` with the following contents:
~~~~ 
#$ -l mfree=100G

module load julia/latest
module load coreutils/8.24
path/to/mergeall.py -B [list of split.q10.sort.bams] -O [Path to output folder] -P [Prefix you want on output files] -C [path to file with each of the valid barcodes allowed in the experiment] --keep_intermediates
~~~~ 

To execute, run `qsub mergecall.sh`

The first argument of mergeall.py is a list of the full paths to the BAMs created as output of Part 1. The second argument is the full path to where you want the output to go. The third argument is optional, and gives each output file an experiment prefix. The last argument is the path to a text file with all the allowed barcodes for your experiment, one barcode per line (see below).

For more details on further arguments, run `python mergeall.py --help`


### Output
a bunch of stuff

There is a helper function in src called `barcode_helper.R` that will assemble your barcode file for you. To run it do:
`Rscript --vanilla path/to/barcode_helper.R [LIST OF BARCODE COMBOS]` The list of barcodes will be letter plus index number combos that were used in your experiment, and that correspond with the PCR and Tn5 plates. For example, if I used Column 5 and row A for one plate and column 7 and row B for a second plate, I would run `Rscript --vanilla path/to/barcode_helper.R A i5 B i7` 
