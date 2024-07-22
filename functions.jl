using DataFrames, JSON3, Printf
include("types.jl")

global R = 8.31446261815324
<<<<<<< HEAD
# global COMP = ["SIO2", "CAO", "AL2O3", "FEO", "MGO", "NA2O"]
global COMP = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
=======
global COMP = ["SIO2", "CAO", "AL2O3", "FEO", "MGO", "NA2O"]
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

function read_data(fname::String)
    return JSON3.read(fname, Vector{Phase}) |> DataFrame
end

function restructure(s::DataFrame, m::Float64, v::Vector{Vector{Float64}})
    return (id = s.id[1], fml = s.fml[1], F0 = s.F0[1], n = s.n[1], V0 = s.V0[1], K0 = s.K0[1], Kp = s.Kp[1], Θ0 = s.Θ0[1], γ0 = s.γ0[1], q0 = s.q0[1], ηS0 = s.ηS0[1], cme = s.cme[1], molar_fraction = m, sites_cmp = v)
end

<<<<<<< HEAD
function read_models(fname::String, data::DataFrame, model_names::Vector{String}, endmember_fractions::Vector{Dict{String, Float64}})
=======
function read_models(fname::String, data::DataFrame, model_names::Vector{String})
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
    read_models = JSON3.read(fname, Vector{ModelJSON})                          # read the json
    models = Vector{Model}()                                                    # create a vector of Models
    
    for model in read_models                                                    # for each model in read_models
        aux_data = DataFrame()                                                  # create an empty dataframe
        aux_model = DataFrame()
        model_number = 1
        
        if model.name in model_names
            for (cont, em) in enumerate(model.endmembers)                                          # put each endmember in a Dict to be added to Model
                p = findfirst(x -> x == em[1], data.id)
                aux_fraction = endmember_fractions[model_number][em[1]]
                aux_model = restructure(DataFrame(data[p, :]), aux_fraction, em[2])
                push!(aux_data, aux_model)
            end

            aux_data, n_sites, site_multiplicities = set_multiplicity(aux_data)
            push!(models, Model(model.name, n_sites, site_multiplicities, aux_data, model.margules, model.van_laar))
            
        else
            continue
        end

        model_number += 1

    end
    return models
end

<<<<<<< HEAD
function set_multiplicity(data::DataFrame)

    n_endmembers = size(data)[1]                                                # get the number of endmembers
    n_sites = size(data.sites_cmp[1])[1]                                        # get the number of sites
    aux_max = zeros(n_endmembers, n_sites)                                      # matrix to store maximum values
    multiplicity = ones(n_sites)                                                # vactor to store multiplicities

    for j in 1:n_sites
        for i in 1:n_endmembers
            aux_max[i, j] = maximum(data.sites_cmp[i][j])
        end
    end

    for j in 1:n_sites
        for i in 1:n_endmembers
            multiplicity[j] = maximum(aux_max[:, j])
            data.sites_cmp[i][j] ./= multiplicity[j]
        end
    end

    return data, n_sites, multiplicity
end

=======
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
function message(str::String, arg::Vector{Float64}=[0.0])
    max_dist = 80
    dist = 3

    if str == "line"
        println(repeat("=", max_dist))
    elseif str == "start"
        message("line")
        @printf("Starting computation for \n")
        @printf(" * P: %.2f bar", arg[1])
        @printf(", T: %.2f K\n", arg[2])
        message("line")
    elseif str == "PT"
        @printf(" * P: %.2f bar", arg[1])
        @printf(", T: %.2f K\n", arg[2])
    elseif str == "gibbs"
        @printf(" * gibbs: \t%15.2f\n", arg[1])
    elseif str == "config"
        @printf(" * R*T*config: \t%15.2f\n", arg[1])
    elseif str == "excess"
        @printf(" * excess: \t%15.2f\n", arg[1])
    elseif str == "μi"
        @printf(" * μ: \t\t%15.2f\n", arg[1])
    elseif str == "μ"
        @printf(" * μ: \t\t%15.2f\n", arg[1])
    else
        println(repeat("=", dist), str, repeat("=", max_dist - dist - length(str)))
    end
end

<<<<<<< HEAD
# function make_comp(comp::Dict{String, Float64})
#     sc = size(COMP)
#     my_comp = zeros(sc)
#     for key in keys(comp)
#         p = findfirst(x -> x == key, COMP)
#         my_comp[p] = comp[key]
#     end
#     return my_comp
# end
=======
function make_comp(comp::Dict{String, Float64})
    sc = size(COMP)
    my_comp = zeros(sc)
    for key in keys(comp)
        p = findfirst(x -> x == key, COMP)
        my_comp[p] = comp[key]
    end
    return my_comp
end
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

"""
    plg(t)

This function evaluates the Debye integral:
    int((ln(1-exp(-t))*t^2),t=0..t)

# Arguments
- `t::Float64`: Temperature.

# Returns
- `plg::Float64`: Value of the integral
"""
function plg(t)
    p0 = exp(-t)
    p1 = 1.0
    p2 = t * t
    p3 = 2.0 * t
    nopt50 = 4.7637509757004091e-015
    # c 45/Pi
    plg = -2.1646464674222763831

    i = 1
    while i < 100000

        p4 = Float64(i)
        p1 = p1 * p0
        dinc = (p2 + (p3 + 2.0 / p4) / p4) * p1 / p4 / p4
        plg = plg + dinc

        if (abs(dinc / (1.0 + abs(plg))) < nopt50)
            return plg
        end

        i += 1
    end
    return plg
end

"""
    gibbs(t, p, phase)

This function calculates the Gibbs energy of the `phase`

# Arguments
- `phase::Phase`: Phase object containing the phase data.
- `t::Float64`: Temperature value in K.
- `p::Float64`: Pressure value in bar.

# Returns
- `G::Float64`: Gibbs energy value.
"""
function calc_gibbs(phase::DataFrameRow{DataFrame, DataFrames.Index}, p::Float64=1000.0, t::Float64=1000.0)
    
    tr = 300.0

    v0 = -phase.V0
    nr9 = -9.0 * phase.n * R
    nr9t0 = nr9 * tr
    c1 = -9.0 * phase.V0 * phase.K0
    c2 = phase.Kp / 2.0 - 2.0
    c3 = 3.0 * c1 * c2
    aii = 6.0 * phase.γ0
    aiikk2 = 0.5 * aii * (-2.0 + 6.0 * phase.γ0 - 3.0 * phase.q0)
    aii2 = 3 * phase.γ0

    # aiikk  = thermo[16]

    r23 = 2.0 / 3.0
    r59 = 5.0 / 9.0

    iopt21 = 100
    nopt51 = 1.8590370495272219e-13

    t1 = phase.Θ0 / t
    t2 = t / tr
    nr9t = nr9 * t
    # initial guess for volume:
    #   taylor(diff(FTH,v),v=v0,1)
    #   JADC Feb 26, 2008.
    #   the dfth0 could be loaded as a constant.
    tht = t1
    tht0 = tht * t2

    k00 = phase.K0
    k0p = phase.Kp
    γ0 = phase.γ0

    dfth = nr9t * γ0 / v0 * (3.0 * plg(tht) / tht^3 - log(1.0 - exp(-tht)))
    dfth0 = nr9t0 * γ0 / v0 * (3.0 * plg(tht0) / tht0^3 - log(1.0 - exp(-tht0)))
    #  taylor(diff(FC,v),v=v0,2)
    #       v = (k00-dfth+dfth0-p)/k00*v0
    #  taylor(diff(FC,v),v=v0,3)
    root = k00 * ((2.0 + 2.0 * k0p) * (p + dfth - dfth0) + k00)

    v = 0.0
    if (root > 0.0)
        v = ((2.0 + k0p) - sqrt(root) / k00) * v0 / (1.0 + k0p)
        if (v < v0 / 10.0) || (v > v0 * 10.0)
            v = v0
        end
    else
        v = v0
    end

    itic = 0
    ibad = 4
    bad = true

    while (itic < 100)

        itic += 1
        # println("itic: ", itic)
        #  f, and derivatives
        v23 = (v0 / v)^r23
        f = 0.5 * v23 - 0.5
        df = -v23 / v / 3.0
        d2f = r59 * v23 / v^2
        # cold part derivatives
        dfc = (c3 * f + c1) * f * df
        d2fc = (2.0 * c3 * f + c1) * df^2 + (c3 * f + c1) * f * d2f
        # debye T/T (tht)
        z = 1.0 + (aii + aiikk2 * f) * f

        if (z < 0.0 || v / v0 > 100.0 || v / v0 < 1e-2)
            println("ERROR z or v/v0")
        end
        # println(z)

        root = sqrt(z)

        tht = t1 * root
        tht0 = tht * t / tr
        # tht derivatives
        a2f = aii2 + aiikk2 * f
        da = a2f / root
        dtht = t1 * da * df
        d2tht = t1 * ((aiikk2 / root - a2f^2 / z^1.5) * df^2 + da * d2f)
        # println(d2tht)

        dtht0 = dtht * t2
        d2tht0 = d2tht * t2
        # polylog functions:
        fpoly = 3.0 * plg(tht) / tht^3
        fpoly0 = 3.0 * plg(tht0) / tht0^3
        # thermal part derivatives:
        etht = exp(-tht)

        if (1.0 - etht < 0.0)
            println("ERROR 1-etht")
        end

        letht = log(1.0 - etht)

        dfth = (letht - fpoly) * nr9t * dtht / tht
        d2fth = ((4.0 * dtht^2 / tht - d2tht) * (fpoly - letht) + dtht^2 * etht / (1.0 - etht)) * nr9t / tht

        etht0 = exp(-tht0)

        if (1.0 - etht0 < 0.0)
            println("ERROR 1-tht0")
        end

        letht0 = log(1.0 - etht0)

        dfth0 = (letht0 - fpoly0) * nr9t0 * dtht0 / tht0
        d2fth0 = ((4.0 * dtht0^2 / tht0 - d2tht0) * (fpoly0 - letht0) + dtht0^2 * etht0 / (1.0 - etht0)) * nr9t0 / tht0

        f1 = -dfc - dfth + dfth0 - p
        df1 = -d2fc - d2fth + d2fth0
        dv = f1 / df1

        if (v - dv < 0.0)
            dv = v / 2.0
        end

        v -= dv
<<<<<<< HEAD
        
=======
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b

        if (itic > iopt21 || abs(f1) > 1e40)
            # allow bad result
            if (abs(f1 / p) < 0.0)
                ibad = 5
                println("ERROR abs(f1/p)")
            end
        elseif (abs(dv / (1.0 + v)) < nopt51)
            bad = false
            break
        end
    end

    # if (bad)
    # if we get here, failed to converge
    #       if (izap.le.iopt(1)) then

    #       msg = 'STXGJI/'//names(id)

    #       call conwrn (ibad,msg)

    #       izap = izap + 1

    #       if (izap.eq.iopt(1)) call warn (49,r,93,'STXGJI')

    #       end if

    #       if (ibad.eq.4) then 
    # # c                                 destabilize the phase
    #       gstxgi  = 1d2*p
    #       badend(id) = .true.
    #       return
    #       end if
    # end if

    # get helmoltz energy:
    f = 0.5 * (v0 / v)^r23 - 0.5
    z = 1.0 + (aii + aiikk2 * f) * f
    root = sqrt(z)
    # final estimate for tht
    tht = t1 * root
    tht0 = tht * t2
    # helmholtz energy
    a = phase.F0 + c1 * f^2 * (0.5 + c2 * f) + nr9 * (t / tht^3 * plg(tht) - tr / tht0^3 * plg(tht0))
    G = a + p * v - t * phase.cme
    return G
end

<<<<<<< HEAD
function no_nan_sum(vec::Vector{Float64})
    return sum([x for x in vec if !isnan(x)])
end

function calc_config(phase::DataFrameRow{DataFrame, DataFrames.Index}, index::Int64, model::Model)
    # initialization
    config::Vector{Float64} = []
    n_endmembers            = size(model.endmembers)[1]
    total_molar_fraction    = sum(model.endmembers.molar_fraction)
    n_sites                 = model.sites
    n_species               = size(model.endmembers.sites_cmp[1][1])[1]

    for site in 1:n_sites
        
        molar_fraction::Vector{Float64} = []
        aux_config::Vector{Float64}     = []

        for specie in 1:n_species
            push!(molar_fraction, 0.0)
            
            for endmember in 1:n_endmembers
                
                if model.endmembers.sites_cmp[endmember][site][specie] != 0
                    # println(model.endmembers.id[endmember])
                    # println("endmember_fractions: ", model.endmembers.molar_fraction[endmember])
                    molar_fraction[end] += model.endmembers.molar_fraction[endmember]
                end                
            end
            
            # println("model.site_multiplicities[site]: ", model.site_multiplicities[site])
            aux_config = model.site_multiplicities[site] .* molar_fraction .* log.(molar_fraction ./ total_molar_fraction)
            # println("aux_config: ", aux_config)
            
        end
        # println("molar_fraction: ", molar_fraction)
        # println(aux_config)
        push!(config, no_nan_sum(aux_config))
        # println(repeat("-", 80))
    end
    
    # println(sum(config))
    # println(config)
    return sum(config)
end

function eye(i::Int64, j::Int64)
    return i == j ? 1.0 : 0.0
end

function calc_excess(i::Int64, model::Model)
    
    excess = 0.0
    n_endmembers = size(model.endmembers)[1]
    W = [value for value in values(model.margules)]
    v = model.van_laar

    asymmetric = false
    if any(value -> value != 1, v)
        asymmetric = true
        # print("asymmetric")
    end

    if asymmetric
        sum_v = 0.0
        for i in 1:n_endmembers
            sum_v += model.endmembers.molar_fraction[i] * v[i]
        end

        mat_phi = zeros(n_endmembers)
        for i in 1:n_endmembers
            mat_phi[i] = (model.endmembers.molar_fraction[i] * v[i]) / sum_v
        end
    end

    excess = 0.0
    it = 1
    for j in 1:n_endmembers-1
        for k in j+1:n_endmembers
            if asymmetric
                excess += (eye(i,j) - mat_phi[j]) * (eye(i,k) - mat_phi[k]) * (W[it] * 2.0 * v[i] / (v[j] + v[k]))
            else
                excess += (eye(i,j) - model.endmembers.molar_fraction[j]) * (eye(i,k) - model.endmembers.molar_fraction[k]) * W[it]
            end
            it += 1
        end
    end

    return excess
end

<<<<<<< HEAD
function gcalc(pressure::Float64, temperature::Float64, models::Vector{Model}, endmember_fractions::Vector{Dict{String, Float64}})
=======
function gcalc(pressure::Float64, temperature::Float64, models::Vector{Model}, endmembers_fractions::Vector{Vector{Float64}})
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
    
    μ = Vector{Float64}()
    
    for (m, model) in enumerate(models)
        gi = Vector{Float64}() 
        ai = Vector{Float64}()    
        xi = Vector{Float64}()
        μi = Vector{Float64}()
<<<<<<< HEAD
        # n_endmembers = size(model.endmembers)[1]
        # for i in 1:n_endmembers
        #     phase = model.endmembers[i, :]
        #     println(phase)
        #     title = " " * string(endmember_fractions[m][i] * 100.0) * " % of `" * phase.id * "` (" * phase.fml * ") "
        #     # message(title);
        #     g = calc_gibbs(phase, pressure, temperature)
        #     # message("gibbs", [g])
        #     push!(gi, g)
        #     println("ef: ", endmember_fractions[m])
        #     a = R * temperature * calc_config(phase, i, model, endmember_fractions[m]) 
        #     # message("config", [a])
        #     push!(ai, a)
        #     e = calc_excess(i, model, endmember_fractions[m])
        #     # message("excess", [e])
        #     push!(xi, e)
        #     # message("μi", [g - a - e]) # message("μi", [g - (a/3.0) - e])
        #     push!(μi, (g - a - e)) # push!(μi, (g - (a/3.0) - e))
        # end
        for (i, (key, value)) in enumerate(endmember_fractions[m])
            phase = model.endmembers[i, :]
            title = " " * string(endmember_fractions[m][key] * 100.0) * " % of `" * phase.id * "` (" * phase.fml * ") "
            message(title);
            g = calc_gibbs(phase, pressure, temperature)
            message("gibbs", [g])
            push!(gi, g)
            # println("ef: ", endmember_fractions[m])
            a = R * temperature * calc_config(phase, i, model) 
            message("config", [a])
            push!(ai, a)
            e = calc_excess(i, model)
            message("excess", [e])
            push!(xi, e)
            # message("μi", [g - a - e]) # message("μi", [g - (a/3.0) - e])
            push!(μi, (g - a - e)) # push!(μi, (g - (a/3.0) - e))
        end
        
        push!(μ, sum(μi .* model.endmembers.molar_fraction[m]))
        message("line")
        message("μ", [sum(μ)])
        message("line")    
    end

    return μ
end

function span_gcalc(n_span::Int64, pressure::Float64, temperature::Float64, models::Vector{Model})
    
    g = Vector{Float64}()
    x = [x for x in 0:n_span] ./ n_span

    for a in 0:n_span
        v = a / n_span
<<<<<<< HEAD
        endmember_fractions = [[v, 1-v]]
        push!(g, sum(gcalc(pressure, temperature, models, endmember_fractions)))
=======
        endmembers_fractions = [[v, 1-v]]
        push!(g, sum(gcalc(pressure, temperature, models, endmembers_fractions)))
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
    end

    fig = Figure()
    Axis(
        fig[1, 1], 
        title="Gibss free energy", 
        xticks = 0.0:0.1:1.0, 
        xlabel = "nᵢ", 
        ylabel = "𝒢 (kJ)",
        ytickformat = "{:.2f}"
        )
    xlims!(0, 1)
    lines!(x, g / 1.0e3)
    display(fig)
    
    return g
<<<<<<< HEAD
end

# function jprint(number::Float64, sep::Char=',', dig::Int=2)
#     # format_str = Printf.Format("%." * string(dig) * "f")
#     # str_number = Printf.format(format_str, number)
#     str_number = format(number, commas=true, precision=dig, separator=sep)
#     return str_number
# end
=======
end
>>>>>>> 810b2137dfaab5d6d63e10c27f821e69d2f5324b
