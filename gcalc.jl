using CairoMakie

include("functions.jl")

function main()

    # user input
    pressure, temperature = 1_000.0, 1_000.0
    message("start", [pressure, temperature])

    ############################################################################
    # benchmarked
    ############################################################################

    # model_names = ["plagioclase"] # zeroed 2nd site
    # species_fractions = [[0.5, 0.5]] # anorthite albite

    model_names = ["spinel"]
    species_fractions = [[0.5, 0.5]] # spinel hercynite

    # model_names = ["olivine"]
    # species_fractions = [[0.5, 0.5]] # fosterite fayalite

    # model_names = ["wadsleyite"]
    # species_fractions = [[0.5, 0.5]] # mg-wadsleyite fe-wadsleyite

    # model_names = ["ringwoodite"]
    # species_fractions = [[0.5, 0.5]] # mg-ringwoodite fe-ringwoodite

    # model_names = ["hp-clinopyroxene"] # energy_interaction @ burnman = 0
    # species_fractions = [[0.5, 0.5]] # hp-clinoestantite hp-clinoferrosilite

    # model_names = ["ca-perovskite"]
    # species_fractions = [[1.0]] # ca-peroviskite 

    # model_names = ["coesite"]
    # species_fractions = [[1.0]] # coesite

    # model_names = ["seifertite"]
    # species_fractions = [[1.0]] # seifertite

    ############################################################################
    # failed benchmark
    ############################################################################

    # model_names = ["quartz"]
    # species_fractions = [[1.0]] # quartz

    # model_names = ["stishovite"]
    # species_fractions = [[1.0]] # stishovite

    ############################################################################
    # not benchmarked
    ############################################################################

    # model_names = ["orthopyroxene"] # energy_interaction @ burnman
    # species_fractions = [[0.25, 0.25, 0.25, 0.25]] # enstatite ferrosilite mg-tschermak ortho-diopside

    # model_names = ["clinopyroxene"] # energy_interaction @ burnman + van laar + zeroed 3rd site
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # diopside hedenbergite clinoestantite ca-tschermak jadeite

    # model_names = ["akimotoite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-akimotoite fe-akimotoite corundum

    # model_names = ["garnet"]
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # pyrope almandine grossular mg-majorite jd-majorite 

    # model_names = ["perovskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-peroviskite, fe-peroviskite, fperov

    








    # model_names = ["C2/c pyroxene"]
    # species_fractions = [[0.5, 0.5]] # C2/c, FC2/c

    # model_names = ["magnesio-wuestite"]
    # species_fractions = [[0.5, 0.5]] # per wus

    

    

    # model_names = ["perovskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # peroviskite

    

    

    





    # model_names = ["clinopyroxene"]
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # jadeite, diopside, hedenbergite, clinoestantite, ca-tschermaks

    # model_names = ["akimotoite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # corundum, mg-akimotoite, fe-akimotoite

    # model_names = ["garnet"]
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # garnet, almandine, mg-majorite, pyrope, jd-majorite

    # model_names = ["post-peroviskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # post-peroviskite, fe-post-peroviskite
    
    # read datasets
    data = read_data("data/stx11_data.json")
    models = read_models("data/stx11_solution.json", data, model_names)

    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, species_fractions)
    # gibbs = span_gcalc(10, pressure, temperature, models)

    return gibbs
end

# @code_warntype main()
aux = main()
