###############################################################################
# Summarise spacer configurations
# -----------------------------------------------------------------------------
# Author:  Thomas E. Gorochowski <thomas.gorochowski@bristol.ac.uk>
# Licence: MIT
###############################################################################

using CSV
using DataFrames
using Query
using Glob
using Printf
using DelimitedFiles

###############################################################################
# RUN THE ANALYSIS
###############################################################################

int_names = ["Int2", "Int3", "Int4", "Int5", "Int7", "Int8", "Int9", "Int10",
             "Int11", "Int12", "Int13"] 
int_order = Dict(zip(int_names, collect(1:11)))
summary = Dict()

fs = glob("output/*.csv")
for f in fs
    filename = split(f, "/")[end]
    cur_key = split(filename, ".")[1]
    df = DataFrame(CSV.File(f, header=1, delim=","))
    summary[cur_key] = zeros(length(int_names), 2)
    for row in eachrow(df)
        split_config = split(row.spacer_config, " ")
        for s in split_config
            if s != "None"
                if contains(s, "_")
                    cur_int = split(s, "_")[1]
                    summary[cur_key][int_order[cur_int], 1] += row.count
                    summary[cur_key][int_order[cur_int], 2] += row.count
                else
                    cur_int = s
                    summary[cur_key][int_order[cur_int], 2] += row.count
                end
            end
        end
    end
end

header_flip = string.(int_names, "_flipped")
header_total = string.(int_names, "_total")
header = ["sample"]
append!(header, header_flip)
append!(header, header_total)

open("output/summary.tsv", "w") do io
    out = join(header, '\t')
    write(io, "$out\n")
    for k in keys(summary)
        row_data = append!(string.([k]), string.(summary[k][:,1]))
        append!(row_data, string.(summary[k][:,2]))    
        out = join(row_data, '\t')
        write(io, "$out\n")
    end
end
