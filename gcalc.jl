using CairoMakie

include("functions.jl")

function main()
    message("start")

    # user input
    pressure, temperature = 1_000.0, 1_000.0
    comp = make_comp(Dict("SIO2" => 1.0, "FEO" => 1.0, "MGO" => 1.0))
    model_names = ["olivine"]#, "spinel"]
    amounts = [[0.5, 0.5]]#, [0.2, 0.8]]

    # read datasets
    data = read_data("data/stx11_data.json")
    models = read_models("data/stx11_solution.json", data, model_names)

    # calculate gibbs free energy
    # gibbs = gcalc(pressure, temperature, comp, models, amounts)
    gibbs = span_gcalc(10, pressure, temperature, comp, models)

    return gibbs
end

# @code_warntype main()
aux = main();
