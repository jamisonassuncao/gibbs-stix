using CairoMakie

include("functions.jl")

function main()

    # user input
    pressure, temperature = 1_000.0, 1_000.0 # K
    # pressure, temperature = 300_000.0, 2_000.0 # bar
    message("start", [pressure, temperature])

    ############################################################################
    # benchmarked
    ############################################################################

    # model_names = ["spinel"]
    # endmember_fractions = [[0.5, 0.5]] # spinel hercynite

    # model_names = ["olivine"]
    # endmember_fractions = [[0.3, 0.7]] # fayalite fosterite

    # model_names = ["wadsleyite"]
    # endmember_fractions = [[0.5, 0.5]] # mg-wadsleyite fe-wadsleyite

    # model_names = ["ringwoodite"]
    # endmember_fractions = [[0.5, 0.5]] # mg-ringwoodite fe-ringwoodite

    # model_names = ["hp-clinopyroxene"] # energy_interaction @ burnman = 0
    # endmember_fractions = [[0.5, 0.5]] # hp-clinoestantite hp-clinoferrosilite

    # model_names = ["ca-perovskite"]
    # endmember_fractions = [[1.0]] # ca-peroviskite 

    # model_names = ["coesite"]
    # endmember_fractions = [[1.0]] # coesite

    # model_names = ["seifertite"]
    # endmember_fractions = [[1.0]] # seifertite

    # model_names = ["kyanite"]
    # endmember_fractions = [[1.0]] # kyanite

    # model_names = ["nepheline"]
    # endmember_fractions = [[1.0]] # nepheline

    # model_names = ["orthopyroxene"]
    # endmember_fractions = [[0.25, 0.25, 0.25, 0.25]] # enstatite ferrosilite mg-tschermak ortho-diopside

    ############################################################################
    # failed benchmark
    ############################################################################

    # For plagioclase, the difference between my code and burnman is the value 
    # of R*T*config, this might have something to do with the failed benchmark.
    # model_names = ["plagioclase"] 
    # endmember_fractions = [[0.5, 0.5]] # anorthite albite

    # model_names = ["quartz"]
    # endmember_fractions = [[1.0]] # quartz

    # model_names = ["stishovite"] # diff. < 54
    # endmember_fractions = [[1.0]] # stishovite

    model_names = ["perovskite"] # not working for new approach
    # endmember_fractions = [[0.3, 0.3, 0.4]] # mg-peroviskite, fe-peroviskite, al-peroviskite
    endmember_fractions = [Dict("mg-perovskite" => 0.3, "fe-perovskite"=> 0.3, "al-perovskite" => 0.4)]

    # model_names = ["post-perovskite"]
    # endmember_fractions = [[0.3, 0.3, 0.4]] # mg-post-peroviskite fe-post-peroviskite al-post-peroviskite

    # model_names = ["magnesio-wustite"]
    # endmember_fractions = [[0.5, 0.5]] # periclase wustite

    # model_names = ["ca-ferrite"]
    # endmember_fractions = [[0.3, 0.3, 0.4]] # mg-ca-ferrite fe-ca-ferrite na-ca-ferrite

    # model_names = ["clinopyroxene"] # energy_interaction @ burnman + van laar + zeroed 3rd site + asymmetric
    # endmember_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # diopside hedenbergite clinoestantite ca-tschermak jadeite

    # model_names = ["akimotoite"]
    # endmember_fractions = [[0.3, 0.3, 0.4]] # mg-akimotoite fe-akimotoite corundum

    # model_names = ["garnet"] # asymmetric
    # endmember_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # pyrope almandine grossular mg-majorite jd-majorite 

   ############################################################################
    # run
    ############################################################################

    # read datasets
    data = read_data("data/stx11_data.json")
    models = read_models("data/stx11_solution.json", data, model_names, endmember_fractions)
    # println(models)
    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, endmember_fractions)
    # gibbs = span_gcalc(10, pressure, temperature, models)

    return gibbs, "Si,Ca,Al,Fe,Mg,Na";
end

# @code_warntype main()
aux, ele = main();
