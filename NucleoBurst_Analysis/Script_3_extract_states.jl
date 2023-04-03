###############################################################################
# Extract spacer configuration from sequencing data
# -----------------------------------------------------------------------------
# Author:  Thomas E. Gorochowski <thomas.gorochowski@bristol.ac.uk>
# Licence: MIT
###############################################################################

using CSV
using DataFrames
using Query
using Glob
using Printf
using BioSequences
using BioAlignments
using FASTX

###############################################################################
# PARAMETERS – Should be optimised for particular use case
###############################################################################

# Minimum alignment score for start and end of memory register
# Previously 300, 100
global MIN_ALN_SCORE_REGION = 20
global MIN_ALN_SCORE_SPACER = 0

# Alignemnt parameters (these have been tweaked for 25 bp barcodes)
#global ALN_MATCH = 5
#global ALN_MISMATCH = -4
#global ALN_GAP_OPEN = -4
#global ALN_GAP_EXTEND = -2

# Parameters from Matt's paper (need to adjust the thresholds)
global ALN_MATCH = 2
global ALN_MISMATCH = -3
global ALN_GAP_OPEN = -5
global ALN_GAP_EXTEND = -2


global BASE_FRACTION_CALL_INSERT = 0.5

###############################################################################
# FUNCTIONS TO DO THE HEAVY LIFTING
###############################################################################

function load_spacers(filename)
    ref_seqs = Dict{String,LongDNASeq}()
    reader = FASTA.Reader(open(filename, "r"))
    record = FASTA.Record()
    read_count = 0
    while !eof(reader)
        read!(reader, record)
        ref_seqs[FASTA.identifier(record)] = FASTA.sequence(record)
        read_count = read_count + 1
    end
    close(reader)
    @info string("Processed ", filename, " containing ", read_count, " sequences")
    return ref_seqs
end

function ref_alignment_range(aln)
    start_pos = -1
    end_pos = -1
    aln_anchors = aln.a.aln.anchors
    first_match = true
    for a in aln_anchors
        if a.op == OP_SEQ_MATCH
            if first_match == true
                start_pos = a.refpos
                first_match = false
                end_pos = a.refpos
            else
                end_pos = a.refpos
            end
        end
    end
    return (start_pos, end_pos)
end

function check_length_fix_orientation(test_seq, start_seq, end_seq)
    # Alignment type used and scoring model
    problem = SemiGlobalAlignment()
    scoremodel = AffineGapScoreModel(
                  match=ALN_MATCH,
                  mismatch=ALN_MISMATCH,
                  gap_open=ALN_GAP_OPEN,
                  gap_extend=ALN_GAP_EXTEND
                )
    test_seq_rc = reverse_complement(test_seq)
    alignment_start = pairalign(problem, start_seq, test_seq, scoremodel)
    alignment_start_rc = pairalign(problem, start_seq, test_seq_rc, scoremodel)
    alignment_end = pairalign(problem, end_seq, test_seq, scoremodel)
    alignment_end_rc = pairalign(problem, end_seq, test_seq_rc, scoremodel)
    score_start = score(alignment_start)
    score_start_rc = score(alignment_start_rc)
    score_end = score(alignment_end)
    score_end_rc = score(alignment_end_rc)
    # Test to see what start and end are most likely and return sequence
    if score_start > score_start_rc && score_end > score_end_rc &&
        score_start > MIN_ALN_SCORE_REGION && score_end > MIN_ALN_SCORE_REGION
        (s_idx0, s_idx1) = ref_alignment_range(alignment(alignment_start))
        (e_idx0, e_idx1) = ref_alignment_range(alignment(alignment_end))
        return test_seq[s_idx0:e_idx1]
    elseif score_start_rc > score_start && score_end_rc > score_end &&
        score_start_rc > MIN_ALN_SCORE_REGION && score_end_rc > MIN_ALN_SCORE_REGION
        (s_idx0, s_idx1) = ref_alignment_range(alignment(alignment_start_rc))
        (e_idx0, e_idx1) = ref_alignment_range(alignment(alignment_end_rc))
        return test_seq_rc[s_idx0:e_idx1]
    else
        return nothing
    end
end

function find_spacer_list(seq, ref_seqs)
    spacers = []
    # Alignment type used and scoring model
    problem = SemiGlobalAlignment()
    scoremodel = AffineGapScoreModel(
                  match=ALN_MATCH,
                  mismatch=ALN_MISMATCH,
                  gap_open=ALN_GAP_OPEN,
                  gap_extend=ALN_GAP_EXTEND
                )
    seq_len = length(seq)
    seq_rc = reverse_complement(seq)
    for (ref_key, ref_seq) in ref_seqs
        # Check reference is for a spacer
        if ref_key != "start" && ref_key != "end"
            spacer_aln = pairalign(problem, ref_seq, seq, scoremodel)
            spacer_aln_rc = pairalign(problem, ref_seq, seq_rc, scoremodel)
            spacer_score = score(spacer_aln)
            spacer_score_rc = score(spacer_aln_rc)
            # Check alignments meet minimum score
            if spacer_score > spacer_score_rc && spacer_score > MIN_ALN_SCORE_SPACER
                (s_idx0, s_idx1) = ref_alignment_range(alignment(spacer_aln))
                push!(spacers, [ref_key, s_idx0])
            elseif spacer_score_rc > spacer_score && spacer_score_rc > MIN_ALN_SCORE_SPACER
                (s_idx0, s_idx1) = ref_alignment_range(alignment(spacer_aln_rc))
                push!(spacers, [string(ref_key, "_r"), seq_len-s_idx1])
            end
        end
    end
    # Sort the found spacers on their position
    sort!(spacers, lt=(x,y)->isless(x[2], y[2]))
    # Generate string for configuration

    #=
    # Only return results where all spacers are present
    if length(spacers) == length(ref_seqs) - 2
        spacer_string = ""
        for el in spacers
            if spacer_string == ""
                spacer_string = el[1]
            else
                spacer_string = string(spacer_string, " ", el[1])
            end
        end
        if spacers == []
            return "None"
        else
            return spacer_string
        end
    else
        return "None"
    end
    =#

    spacer_string = ""
    for el in spacers
        if spacer_string == ""
            spacer_string = el[1]
        else
            spacer_string = string(spacer_string, " ", el[1])
        end
    end
    if spacers == []
        return "None"
    else
        return spacer_string
    end
end

function find_spacers(fastq_filename, ref_seqs)
    @info string("Started processing ", fastq_filename, "...\n")
    results = Dict{String, Int}()
    try
        reader = FASTQ.Reader(open(fastq_filename, "r"))
        record = FASTQ.Record()
        read_count = 0
        processed_read_count = 0
        while !eof(reader)
            # Comment out if processing full data
            if read_count >= 10000
                break
            end
            if mod(read_count, 1000) == 0
                print(".")
            end
            read!(reader, record)
            fixed_seq = check_length_fix_orientation(FASTQ.sequence(record), ref_seqs["start"], ref_seqs["end"])
            if fixed_seq != nothing
                # Search for order and orientation of spacers and add to list
                cur_key = find_spacer_list(fixed_seq, ref_seqs)
                if haskey(results, cur_key)
                    results[cur_key] += 1
                else
                    results[cur_key] = 1
                end
                processed_read_count = processed_read_count + 1
            end
            read_count = read_count + 1
        end
        close(reader)
        @info string("Processed ", fastq_filename, " containing ", read_count, " reads (",
                     processed_read_count, " full length memory registers)\n")
    catch y
        @info string("Error in ", fastq_filename, "\n")
        println(y)
    end
    return results
end

function write_results_to_file(results, filename)
    open(filename, "w") do file
        write(file, "spacer_config,count\n")
        for (key, val) in results
            write(file, string(key, ",", val, "\n"))
        end
    end
end

###############################################################################
# RUN THE ANALYSIS
###############################################################################

# Load the metadata for the run
metadata = DataFrame(CSV.File("metadata.csv", header=1, delim=","))

# Create mapping from barcode to sample name using the metadata for run 1
bc_map1 = @from i in metadata begin
    @where i.run == 1
    @select i.barcode=>i.name
    @collect Dict
end

# Generate the input and output file names
read_filenames = []
output_filenames = []

for (bc, name) in bc_map1
    read_filename = "fastq/barcode$(@sprintf("%02d", bc)).fastq"
    output_filename = "output/$(name).csv"
    push!(read_filenames, read_filename)
    push!(output_filenames, output_filename)
end

# Load the spacer sequences to search for
spacer_refs_filename = "spacers.fasta"
ref_seqs = load_spacers(spacer_refs_filename)

# Extract the spacer configurations from each sample
Threads.@threads for idx in 1:length(read_filenames)
    read_filename = read_filenames[idx]
    output_filename = output_filenames[idx]
    results = find_spacers(read_filename, ref_seqs)
    write_results_to_file(results, output_filename)
end
