using DataFrames, JSON3, Printf
include("types.jl")

global R = 8.31446261815324
global COMP = ["SIO2", "CAO", "AL2O3", "FEO", "MGO", "NA2O"]

function read_data(fname::String)
    return JSON3.read(fname, Vector{Phase}) |> DataFrame
end

function restructure(s::DataFrame, v::Vector{Vector{Float64}})
    return (id = s.id[1], fml = s.fml[1], F0 = s.F0[1], n = s.n[1], V0 = s.V0[1], K0 = s.K0[1], Kp = s.Kp[1], Θ0 = s.Θ0[1], γ0 = s.γ0[1], q0 = s.q0[1], ηS0 = s.ηS0[1], cme = s.cme[1], sites_cmp = v)
end

function read_models(fname::String, data::DataFrame, model_names::Vector{String})
    read_models = JSON3.read(fname, Vector{ModelJSON})                          # read the json
    models = Vector{Model}()                                                    # create a vector of Models
    
    for model in read_models                                                    # for each model in read_models
        aux_data = DataFrame()                                                  # create an empty dataframe
        aux_model = DataFrame()
        if model.name in model_names                                            
            for em in model.endmembers
                p = findfirst(x -> x == em[1], data.id)
                aux_model = restructure(DataFrame(data[p, :]), em[2])
                push!(aux_data, aux_model)
            end
            push!(models, Model(model.name, model.sites, aux_data, model.margules, model.van_laar))
        else
            continue
        end
    end
    return models
end

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

function make_comp(comp::Dict{String, Float64})
    sc = size(COMP)
    my_comp = zeros(sc)
    for key in keys(comp)
        p = findfirst(x -> x == key, COMP)
        my_comp[p] = comp[key]
    end
    return my_comp
end

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

function calc_config(phase::DataFrameRow{DataFrame, DataFrames.Index}, index::Int64, model::Model, endmember_fractions::Vector{Float64})
    # initialization
    act = 0.0
    n_endmembers = size(model.endmembers)[1]
    n_sites = model.sites
    n_species = size(model.endmembers.sites_cmp[1][1])[1]
    sijk = model.endmembers.sites_cmp

    for k in 1:n_sites

        sijk1 = sijk[index][k] 

        Sik = sum(sijk1)

        Nk = 0.0
        for c in 1:n_species
            Njk = 0.0
            for i in 1:n_endmembers
                sijk2 = sijk[i][k] 
                Njk += sijk2[c] .* endmember_fractions[i]
            end
            Nk += Njk
        end
        a1 = (Nk != 0) ? Sik * log(Nk) : 0.0                                    # Should I do this if Nk == 0?
        # a1 = Sik * log(Nk)

        a2 = 0.0
        for c in 1:n_species
            Njk = 0.0
            for i in 1:n_endmembers
                sijk2 = sijk[i][k]
                Njk += sijk2[c] .* endmember_fractions[i] #.* aux[k]
            end
            a2 += (Njk != 0) ? sijk1[c] * log(Njk) : 0.0
        end
        act += (a1 - a2)
    end
    return act
end

function eye(i::Int64, j::Int64)
    return i == j ? 1.0 : 0.0
end

function calc_excess(i::Int64, model::Model, endmember_fractions::Vector{Float64})
    
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
            sum_v += endmember_fractions[i] * v[i]
        end

        mat_phi = zeros(n_endmembers)
        for i in 1:n_endmembers
            mat_phi[i] = (endmember_fractions[i] * v[i]) / sum_v
        end
    end

    excess = 0.0
    it = 1
    for j in 1:n_endmembers-1
        for k in j+1:n_endmembers
            if asymmetric
                excess += (eye(i,j) - mat_phi[j]) * (eye(i,k) - mat_phi[k]) * (W[it] * 2.0 * v[i] / (v[j] + v[k]))
            else
                excess += (eye(i,j) - endmember_fractions[j]) * (eye(i,k) - endmember_fractions[k]) * W[it]
            end
            it += 1
        end
    end

    return excess
end

function gcalc(pressure::Float64, temperature::Float64, models::Vector{Model}, endmember_fractions::Vector{Vector{Float64}})
    
    μ = Vector{Float64}()
    
    for (m, model) in enumerate(models)
        gi = Vector{Float64}() 
        ai = Vector{Float64}()    
        xi = Vector{Float64}()
        μi = Vector{Float64}()
        n_endmembers = size(model.endmembers)[1]
        for i in 1:n_endmembers
            phase = model.endmembers[i, :]
            title = " " * string(endmember_fractions[m][i] * 100.0) * " % of `" * phase.id * "` (" * phase.fml * ") "
            message(title);
            g = calc_gibbs(phase, pressure, temperature)
            message("gibbs", [g])
            push!(gi, g)
            a = R * temperature * calc_config(phase, i, model, endmember_fractions[m]) 
            message("config", [a])
            push!(ai, a)
            e = calc_excess(i, model, endmember_fractions[m])
            message("excess", [e])
            push!(xi, e)
            message("μi", [g - a - e]) # message("μi", [g - (a/3.0) - e])
            push!(μi, (g - a - e)) # push!(μi, (g - (a/3.0) - e))
        end
        
        push!(μ, sum(μi .* endmember_fractions[m]))
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
        endmember_fractions = [[v, 1-v]]
        push!(g, sum(gcalc(pressure, temperature, models, endmember_fractions)))
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
end

# function jprint(number::Float64, sep::Char=',', dig::Int=2)
#     # format_str = Printf.Format("%." * string(dig) * "f")
#     # str_number = Printf.format(format_str, number)
#     str_number = format(number, commas=true, precision=dig, separator=sep)
#     return str_number
# end