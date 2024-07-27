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

    model_names = ["spinel"]
    endmember_fractions = [Dict("spinel" => 0.5, "hercynite" => 0.5)] 

    model_names = ["olivine"]
    endmember_fractions = [Dict("fayalite" => 0.3, "forsterite" => 0.7)]

    model_names = ["wadsleyite"]
    endmember_fractions = [Dict("mg-wadsleyite" => 0.5, "fe-wadsleyite" => 0.5)] 

    model_names = ["ringwoodite"]
    endmember_fractions = [Dict("mg-ringwoodite" => 0.5, "fe-ringwoodite" => 0.5)]

    model_names = ["hp-clinopyroxene"] # energy_interaction @ burnman = 0
    endmember_fractions = [Dict("hp-clinoenstatite" => 0.5, "hp-clinoferrosilite" => 0.5)] 

    model_names = ["ca-perovskite"]
    endmember_fractions = [Dict("ca-perovskite" => 1.0)] 

    model_names = ["coesite"]
    endmember_fractions = [Dict("coesite" => 1.0)]

    model_names = ["seifertite"]
    endmember_fractions = [Dict("seifertite" => 1.0)]

    model_names = ["kyanite"]
    endmember_fractions = [Dict("kyanite" => 1.0)] 

    model_names = ["nepheline"]
    endmember_fractions = [Dict("nepheline" =>  1.0)]

    model_names = ["orthopyroxene"]
    endmember_fractions = [Dict("enstatite" => 0.25, "ferrosilite" => 0.25, "mg-tschermak" => 0.25, "ortho-diopside" => 0.25)]

    ############################################################################
    # failed benchmark
    ############################################################################

    model_names = ["plagioclase"] 
    endmember_fractions = [Dict("anorthite" => 0.5, "albite" => 0.5)]

    model_names = ["quartz"]
    endmember_fractions = [Dict("quartz" => 1.0)]

    model_names = ["stishovite"] # diff. < 54
    endmember_fractions = [Dict("stishovite" =>  1.0)]

    # model_names = ["perovskite"] 
    # endmember_fractions = [Dict("mg-perovskite" => 0.3, "fe-perovskite"=> 0.3, "al-perovskite" => 0.4)]

    # model_names = ["post-perovskite"]
    # endmember_fractions = [Dict("mg-post-perovskite" => 0.3, "fe-post-perovskite" => 0.3, "al-post-perovskite" => 0.4)]

    # model_names = ["magnesio-wustite"]
    # endmember_fractions = [Dict("periclase" => 0.5, "wustite" => 0.5)]

    # model_names = ["ca-ferrite"]
    # endmember_fractions = [Dict("mg-ca-ferrite" => 0.3, "fe-ca-ferrite" => 0.3, "na-ca-ferrite" => 0.4)]

    # model_names = ["clinopyroxene"] # energy_interaction @ burnman + van laar + zeroed 3rd site + asymmetric
    # endmember_fractions = [Dict("diopside" => 0.2, "hedenbergite" => 0.2, "clinoenstatite" => 0.2, "ca-tschermak" => 0.2, "jadeite" => 0.2)]

    # model_names = ["akimotoite"]
    # endmember_fractions = [Dict("mg-akimotoite" => 0.3, "fe-akimotoite" => 0.3, "corundum" => 0.4)]

    # model_names = ["garnet"] # asymmetric
    # endmember_fractions = [Dict("pyrope" => 0.2, "almandine" => 0.2, "grossular" => 0.2, "mg-majorite" => 0.2, "jd-majorite" => 0.2)] 

   ############################################################################
    # run
    ############################################################################

    # read datasets
    data = read_data("data/stx11_data.json")
    # println(data)
    models = read_models("data/stx11_solution.json", data, model_names, endmember_fractions)
    print_endmembers(models)
    # print_endmember("")
    # calculate gibbs free energy
    gibbs = gcalc(pressure, temperature, models, endmember_fractions)

    return gibbs, "Si,Ca,Al,Fe,Mg,Na";
end

# @code_warntype main()
aux, ele = main();
