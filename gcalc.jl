using JSON

include("functions.jl")

function activity(my_comp, model)
    actvty = 0
    n_species = size(model["species"])[1]
    n_sites = model["sites"]
    n_components = 1

    for i in 1:n_species
        println("specie: ", model["species"][i])
        # G = gibbs(model["species"][i], t=1000.0, p=1000.0)
        for k in 1:n_sites
            for j in 1:n_components
                println(j, k, i)
            end
        end
    end

    return actvty
end


function main()
    # load endmembers
    species = JSON.parsefile("stx11_data.json")
    # set components
    my_comp = Dict("MGO" => 1.0, "SIO2" => 1.0, "FEO" => 1.0)
    # load solution model
    modelname = "olivine"
    model = read_model(modelname, "stx11_solution.json")

    temperature = 1000.0 # K
    pressure = 1000.0 # bar 

    n_species = size(model["species"])[1]
    G = zeros(n_species)

    i = 1
    for specie in species
        if specie["id"] in model["species"]
            G[i] = gibbs(specie, 1000.0, 1000.0)
            i += 1
        end
    end
    # actv = R * temperature * activity(my_comp, model);
    return i
end



function foo(x::Int64, y::Int64)
    z = x + y
    return z
end

function foo(x::Float64, y::Float64)
    z = x + 2*y
    return z
end

function foo(x::_T, y::_S)  where {_T<:Number,_S<:Number}
    z = x + 3*y
    return z
end