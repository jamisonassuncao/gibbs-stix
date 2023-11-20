PARAMS = [
    
    # "Z, Formula units in unit cell",
    # "Formula mass (g/mol)",
    # "T_0 (K)",
    "F_0 (kJ/mol)",
    "Atoms in formula unit",
    "V_0 (cm^3/mol)",
    "K_0 (GPa)",
    "K_0 prime",
    # "K_0K_0 prime prime (zero for third order)",
    "Theta_0 (K)", 
    # "Debye Acoustic Branch 2",
    # "Debye Acoustic Branch 3",
    # "Sin Acoustic Branch 1",
    # "Sin Acoustic Branch 2",
    # "Sin Acoustic Branch 3",
    # "Einstein Oscillator 1",
    # "Weight of Einstein Oscillator 1",
    # "Einstein Oscillator 2",
    # "Weight of Einstein Oscillator 2",
    # "Einstein Oscillator 3",
    # "Weight of Einstein Oscillator 3",
    # "Einstein Oscillator 4",
    # "Weight of Einstein Oscillator 4",
    # "Upper Limit of Optic Continuum",
    # "Lower Limit of Optic Continuum",
    "gamma_0",
    "q_0",
    # "beta (electronic contribution)",
    # "gammael_0 (electronic contribution)",
    # "q_2A_2 (anharmonic contribution)",
    # "High Temperature Approximation (1=yes)",
    # "Birch-Murnaghan (0) or Vinet (1)",
    # "Einstein (0) or Debye (1)",
    # "Zero Point Pressure (1=yes)",
    # "Ambient Shear Modulus",
    # "Pressure Derivative",
    # "Temperature Derivative",
    # "Critical Temperature (K)",
    # "Critical Entropy (J/mol/K)",
    # "Critical Volume (cm^3/mol)",
    # "Van Laar size parameter (-)",
    # "C12 prime",
    # "C44 prime"
    ]

function main()
    data = Vector{Vector{String}}()
    files = readdir(join=true)          # files names in current dir
    sfile = open("data.json", "w")
    print(sfile, "[")
    for file in files
        if (!occursin(".tex", file) 
            && !occursin(".jl", file) 
            && !occursin(".md", file) 
            && !occursin(".toml", file) 
            && !occursin(".git", file) 
            && !occursin(".json", file) 
            && !occursin("changelog", file)
            && !occursin("out", file)
            && !occursin("phase", file))
            
            bname = basename(file)
            println("cur file is ", bname)
            file = open(file, "r")              # open file
            
            count = 0
            for line in eachline(file)
                line = split(line, limit=2)     # split in first space
                line[2] = lstrip(line[2])       # strip leading spaces
                line[2] = rstrip(line[2])       # strip trailing spaces
                push!(data, line)               # append <line> to <data>
                # println(sfile, line)
                if count == 0
                    printstr = "{\"id\":\"" * bname * "\",\"fml\":\"" * line[1] * "\","
                    print(sfile, printstr)
                else 
                    if line[2] in PARAMS
                        if line[2] == "F_0 (kJ/mol)"
                            value = round(1_000 * parse(Float64, line[1]), digits=5)
                            printstr = "\"F0\":" * "$value,"
                        elseif line[2] == "Atoms in formula unit"
                            printstr = "\"n\":-" * line[1] * ","
                        elseif line[2] == "V_0 (cm^3/mol)"
                            value = round(0.1 * parse(Float64, line[1]), digits=5)
                            printstr = "\"V0\":" * "$value,"
                        elseif line[2] == "K_0 (GPa)"
                            printstr = "\"K0\":" * line[1] * ","
                        elseif line[2] == "K_0 prime"
                            printstr = "\"Kp\":" * line[1] * ","
                        elseif line[2] == "Theta_0 (K)"
                            printstr = "\"θ0\":" * line[1] * ","
                        elseif line[2] == "gamma_0"
                            printstr = "\"γ0\":" * line[1] * ","
                        elseif line[2] == "q_0"
                            printstr = "\"q0\":" * line[1] * ","
                        # elseif line[2] == ""
                            # printstr = "\"ηS0\":" * line[1] * ","
                            # printstr = "\"cme\":" * line[1] * ","
                        end

                        print(sfile, printstr)
                    end
                end
                count += 1
                
            end
            print(sfile, "},\n")
            print(count-1, " files were processed\n")
            close(file)                         # close file
        end
    end
    print(sfile, "]")
    close(sfile)
    return nothing
end

out = main();