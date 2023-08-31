using JSON

include("functions.jl")

function activity(specie, model)
    actvty = 1.0
    
    n_species = length(model["species"])
    n_sites = model["sites"]
    
    Sik = 0.0 # sum of the stochiometric coefficient (sijk) of component j on site k in species i
    logNk = 0.0
    
    sijk = values(specie["cmp"])

    for _ in n_sites
        Njk = 0.0
        Sik += sum(sijk) 
        for _ in n_species
            Njk += 1
        end
        logNk += log(sum(Njk))
    end

    activity 

    # for i in 1:n_species
    #     println("specie: ", model["species"][i])
    #     # G = gibbs(model["species"][i], t=1000.0, p=1000.0)
    #     for k in 1:n_sites
    #         for j in 1:n_components
    #             println(j, k, i)
    #         end
    #     end
    # end

    return actvty
end

function make_comp(comp)
    # ["SIO2", "MGO", "FEO", "CAO", "AL2O3", "NA2O"]
    sc = size(cmp)
    my_comp = zeros(sc)
    for key in keys(comp)
        p = findfirst(x -> x == key, cmp)
        my_comp[p] = comp[key]
    end
    println(my_comp)
    return my_comp
end


function main()
    # load endmembers
    species = JSON.parsefile("stx11_data.json")
    # set components
    my_comp = make_comp(Dict("SIO2" => 1.0, "MGO" => 10.0, "FEO" => 100.0))
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
            G[i] = gibbs(specie, temperature, pressure)
            a = R * temperature * activity(specie, model)
            println("a: ", a)
            i += 1
        end
    end
    # actv = R * temperature * activity(my_comp, model);
    return nothing
end

@time test = main()

