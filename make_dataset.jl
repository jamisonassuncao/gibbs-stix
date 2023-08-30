using JSON
using DataFrames

include("functions.jl")

filename = "stx11.json"  # Replace with the path to your JSON file
data = open(filename, "r") do file
    read(file, String)
end
data = JSON.parse(data)
data = DataFrame(data)

output_filename = "stx11_data.json"
open(output_filename, "w") do io
    write(io, "[")
    for (index_row, row) in enumerate(eachrow(data))
        # count components
        comps = split(row.cmp, r"[()]")
        pop!(comps) # pop what is after the `)`
        for c in 1:2:size(comps)[1]
            p = findfirst(x -> x == comps[c], cmp)
            if isnothing(p)
                error("Invalid chemical component")
            end
            n_cmp[p] = parse(Float64, comps[c+1])
        end
        cmp_dict = Dict(zip(cmp, n_cmp))

        phase = Phase(
            row.id,
            row.fml,
            cmp_dict,
            row.G0,
            row.S0,
            row.V0,
            row.c1,
            row.c2,
            row.c3,
            row.c4,
            row.c5,
            row.c6,
            row.c7
        )

        write(io, JSON.json(phase))
        if index_row != nrow(data)
            write(io, ",\n")
        end

    end

    write(io, "]")
end