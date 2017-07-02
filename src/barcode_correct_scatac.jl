##
##
##    barcode_correct_scatac.jl
##
##  Barcode correction for sci-ATAC-seq
##
##

using Levenshtein
using GZip
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--fastq", "-f"
            help = "Point to fastqs to be corrected"
            required = true
        "--output", "-o"
            help = "Output prefix"
            required = true
        "--miseq"
            help = "Was this a miseq run?"
            action = :store_true
        "--maxedit", "-e"
            help = "Max allowed edit distance"
            arg_type = Int
            default = 3
        "--minrunnerup", "-m"
            help = "Minimum edit distance between first and second best matches. Prevents ambiguous barcode assignments"
            arg_type = Int
            default = 2
    end
    return parse_args(s)
end

# Standard nextera barcodes
const NEX_I7 = ["ATTACTCG", "TCCGGAGA", "CGCTCATT", "GAGATTCC", "ATTCAGAA",
    "GAATTCGT", "CTGAAGCT", "TAATGCGC", "CGGCTATG", "TCCGCGAA", "TCTCGCGC",
    "AGCGATAG"]

const PCR_I7 = ["TCGGATTCGG" ,"GCGGCTGCGG", "AGATTACGTT", "CTAACTAGGT", "CATAGCGACC",
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

const PCR_I5 = ["CTCCATCGAG", "TTGGTAGTCG", "GGCCGTCAAC", "CCTAGACGAG", "TCGTTAGAGC",
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

const NEX_I5 = ["TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA", "AGGCGAAG",
    "TAATCTTA", "CAGGACGT", "GTACTGAC"]

const SEQ_DICT = Dict('A'=>'T','T'=>'A','G'=>'C','C'=>'G','N'=>'N')

function reverseComplement(seq)
    return(map( base -> SEQ_DICT[base], reverse(seq)))
end

const FULL_BARC = Dict{AbstractString,Int64}()
for i in NEX_I7
    for j in PCR_I7
        for k in PCR_I5
            for l in NEX_I5
                currbarc = i * j * k * l
                FULL_BARC[currbarc] = 0
            end
        end
    end
end

function do_correction(bc::AbstractString, possibles, max_edit::Int64, min_runner_up::Int64)
    if bc in possibles
        return(bc)
    end
    winner_ed::Int64 = max_edit + 1
    runner_up_ed::Int64 = max_edit + min_runner_up + 1
    winner::AbstractString = ""
    for barcode in possibles
        current::Int64 = levenshtein(bc, barcode)
        if current <= winner_ed
            runnerup_ed = winner_ed
            winner = barcode
            winner_ed = current
        end
    end

    if winner_ed > max_edit || runner_up_ed - winner_ed < min_runner_up
        return("barcode_fail")
    else
        return(winner)
    end
end

function main()
    perfect = 0
    corrected = 0
    fail = 0
    total = 0
    parsed_args = parse_commandline()
    fastq_files = readdir(parsed_args["fastq"])
    max_edit = parsed_args["maxedit"]
    min_runner_up = parsed_args["minrunnerup"]

    ifile = string(parsed_args["fastq"], "/", AbstractString(fastq_files[find(x -> contains(x, "R1"), fastq_files)][1]))

    fileR1 = ifile
    fileR2 = replace(ifile, "R1","R2")
    o = GZip.gzopen(string(fileR1, ".out.fq.gz"), "w")
    g = GZip.gzopen(string(fileR2, ".out.fq.gz"), "w")
    log = open(string(parsed_args["fastq"], "/barcode_correct.log", ), "w")
    write(log, "failed\tcorrected\tperfect\ttotal\n")
    GZip.gzopen(fileR1, "r") do f
        GZip.gzopen(fileR2, "r") do r
            while !eof(f)
                total += 1
                tag_line::AbstractString = readline(f)
                read_line::AbstractString = readline(f)
                plus_line::AbstractString = readline(f)
                qual_line::AbstractString = readline(f)
                tag_line = strip(replace(split(split(tag_line, " ")[end], ":")[end], "+", ""))
                tag_line2::AbstractString = readline(r)
                read_line2::AbstractString = readline(r)
                plus_line2::AbstractString = readline(r)
                qual_line2::AbstractString = readline(r)
                if length(tag_line) != 36
                    fail += 1
                    continue
                end
                b1::AbstractString = tag_line[1:8]
                b2::AbstractString = tag_line[9:18]
                if(!parsed_args["miseq"])
                    b4 = reverseComplement(tag_line[19:26])
                    b3 = reverseComplement(tag_line[27:36])
                else
                    b3 = tag_line[19:28]
                    b4 = tag_line[29:36]
                end
                ori_barcode::AbstractString = b1 * b2 * b3 * b4
                if haskey(FULL_BARC, ori_barcode)
                    perfect += 1
                    content = string("@", ori_barcode, ":0\n", read_line,
                        plus_line,  qual_line)
                    content2 = string("@", ori_barcode, ":0\n", read_line2,
                        plus_line2, qual_line2)
                    write(o, content)
                    write(g, content2)
                    continue
                end
                b1_cor::AbstractString = do_correction(b1, NEX_I7, max_edit, min_runner_up)
                if b1_cor == "barcode_fail"
                    fail += 1
                    continue
                end
                b2_cor::AbstractString = do_correction(b2, PCR_I7, max_edit, min_runner_up)
                if b2_cor == "barcode_fail"
                    fail += 1
                    continue
                end
                b3_cor::AbstractString = do_correction(b3, PCR_I5, max_edit, min_runner_up)
                if b3_cor == "barcode_fail"
                    fail += 1
                    continue
                end
                b4_cor::AbstractString = do_correction(b4, NEX_I5, max_edit, min_runner_up)
                if b4_cor == "barcode_fail"
                    fail += 1
                    continue
                end
                corrected += 1
                cor_barcode::AbstractString = b1_cor * b2_cor * b3_cor * b4_cor
                edit_dist = levenshtein(cor_barcode, ori_barcode)
                content::AbstractString = string("@", cor_barcode, ":",
                    string(edit_dist), "\n", read_line,
                    plus_line,  qual_line)
                content2::AbstractString = string("@", cor_barcode, ":",
                    string(edit_dist), "\n", read_line2,
                    plus_line2, qual_line2)
                write(o, content)
                write(g, content2)
            end
        end
    end
    close(o)
    close(g)
    write(log, string(fail, "\t", corrected, "\t", perfect, "\t", total, "\n"))
    close(log)
end


main()
