function main()
    data = Vector{Vector{String}}()
    files = readdir(join=true)          # files names in current dir
    for file in files
        
        if (!occursin(".tex", file) 
            && !occursin(".jl", file) 
            && !occursin(".md", file) 
            && !occursin(".toml", file) 
            && !occursin(".git", file) 
            && !occursin("changelog", file)
            && !occursin("out", file))
            
            
            file = open(file, "r")              # open file
            
            for line in eachline(file)
                line = split(line, limit=2)     # split in first space
                line[2] = lstrip(line[2])       # strip leading spaces
                line[2] = rstrip(line[2])       # strip trailing spaces
                push!(data, line)               # append <line> to <data>
            end
            close(file)                         # close file
            # println(data)
        end
    end
    # println(file)
    return data
end

out = main();