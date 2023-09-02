include("functions.jl")

function main()
    message("start")

    # user input
    pressure, temperature = 1_000.0, 1_000.0
    comp = make_comp(Dict("SIO2" => 1.0, "FEO" => 100.0, "MGO" => 10.0))
    model_names = ["olivine", "spinel"]
    amounts = [[0.1, 0.9], [0.2, 0.8]]

    # read datasets
    data = read_data("stx11_data.json")
    models = read_models("stx11_solution.json", data, model_names)

    # calculate gibbs free energy 
    gibbs = gcalc(pressure, temperature, comp, models, amounts)

    message("line")
    return models
end

# @code_warntype main()
# @time main()
aux = main();
