using CairoMakie

include("functions.jl")

function main()
    message("start")

    # user input
    pressure, temperature = 1_000.0, 1_000.0

    # model_name = ["olivine"]
    # comp = make_comp(Dict("FEO" => 1.0, "MGO" => 1.0)) # olivine
    # species_fractions = [[0.5, 0.5]] # olivine

    # model_names = ["perovskite"]
    # comp = make_comp(Dict("FEO" => 1.0, "MGO" => 1.0)) # peroviskite
    # species_fractions = [[0.3, 0.3, 0.4]] # peroviskite

    model_names = ["spinel"]
    comp = make_comp(Dict("AL2O3" => 2.0, "FEO" => 1.0, "MGO" => 1.0)) # spinel
    species_fractions = [[0.5, 0.5]] # spinel

    # read datasets
    data = read_data("data/stx11_data.json")
    models = read_models("data/stx11_solution.json", data, model_names)

    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, comp, models, species_fractions) 
    # gibbs = span_gcalc(10, pressure, temperature, comp, models)

    return gibbs
end

# @code_warntype main()
aux = main();
