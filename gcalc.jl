using CairoMakie

include("functions.jl")

function main()

    # user input
    pressure, temperature = 1_000.0, 1_000.0
    message("start", [pressure, temperature])

    # model_names = ["C2/c pyroxene"]
    # species_fractions = [[0.5, 0.5]] # C2/c pyroxene

    # model_names = ["olivine"]
    # species_fractions = [[0.5, 0.5]] # olivine

    # model_names = ["perovskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # peroviskite

    model_names = ["spinel"]
    species_fractions = [[1.0, 0.0]] # 100% spinel
    # species_fractions = [[0.5, 0.5]] # spinel
    
    # read datasets
    data = read_data("data/stx11_data.json")
    models = read_models("data/stx11_solution.json", data, model_names)

    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, species_fractions) 
    # gibbs = span_gcalc(10, pressure, temperature, models)

    return gibbs
end

# @code_warntype main()
aux = main();
