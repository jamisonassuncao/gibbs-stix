using CairoMakie

include("functions.jl")

function main()

    # user input
<<<<<<< HEAD
    pressure, temperature = 1_000.0, 1_000.0 # K
    # pressure, temperature = 300_000.0, 2_000.0 # bar
=======
    pressure, temperature = 1_000.0, 1_000.0
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
    message("start", [pressure, temperature])

    ############################################################################
    # benchmarked
    ############################################################################

<<<<<<< HEAD
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
=======
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

    # model_names = ["kyanite"]
    # species_fractions = [[1.0]] # kyanite

    # model_names = ["nepheline"]
    # species_fractions = [[1.0]] # nepheline

    ############################################################################
    # benchmarked with different criteria for sites
    ############################################################################

    # model_names = ["orthopyroxene"]
    # species_fractions = [[0.25, 0.25, 0.25, 0.25]] # enstatite ferrosilite mg-tschermak ortho-diopside
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

    ############################################################################
    # failed benchmark
    ############################################################################

    # For plagioclase, the difference between my code and burnman is the value 
    # of R*T*config, this might have something to do with the failed benchmark.
    # model_names = ["plagioclase"] 
<<<<<<< HEAD
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
=======
    # species_fractions = [[0.5, 0.5]] # anorthite albite

    # model_names = ["quartz"]
    # species_fractions = [[1.0]] # quartz

    # model_names = ["stishovite"] # diff. < 54
    # species_fractions = [[1.0]] # stishovite

    # model_names = ["perovskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-peroviskite, fe-peroviskite, al-peroviskite

    # model_names = ["post-perovskite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-post-peroviskite fe-post-peroviskite al-post-peroviskite

    # model_names = ["magnesio-wustite"]
    # species_fractions = [[0.5, 0.5]] # periclase wustite

    # model_names = ["ca-ferrite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-ca-ferrite fe-ca-ferrite na-ca-ferrite

    # model_names = ["clinopyroxene"] # energy_interaction @ burnman + van laar + zeroed 3rd site + asymmetric
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # diopside hedenbergite clinoestantite ca-tschermak jadeite

    # model_names = ["akimotoite"]
    # species_fractions = [[0.3, 0.3, 0.4]] # mg-akimotoite fe-akimotoite corundum

    # model_names = ["garnet"] # asymmetric
    # species_fractions = [[0.2, 0.2, 0.2, 0.2, 0.2]] # pyrope almandine grossular mg-majorite jd-majorite 
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

   ############################################################################
    # run
    ############################################################################

    # read datasets
    data = read_data("data/stx11_data.json")
<<<<<<< HEAD
    models = read_models("data/stx11_solution.json", data, model_names, endmember_fractions)
    # println(models)
    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, endmember_fractions)
    # gibbs = span_gcalc(10, pressure, temperature, models)

    return gibbs, "Si,Ca,Al,Fe,Mg,Na";
end

# @code_warntype main()
aux, ele = main();
=======
    models = read_models("data/stx11_solution.json", data, model_names)

    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, species_fractions)
    # gibbs = span_gcalc(10, pressure, temperature, models)

    return gibbs
end

# @code_warntype main()
aux = main()
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
