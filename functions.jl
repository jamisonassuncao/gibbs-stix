using DataFrames, JSON3, Printf

global R = 8.31446261815324
global COMP = ["CAO", "AL2O3", "FEO", "MGO", "NA2O"]

struct Phase
    id::String                  # Name id
    fml::String                 # Chemical formula
    cmp::Dict{String, Float64}  # Chemical components (SiO2, MgO, FeO, CaO, Al2O3, Na2O)
    F0::Float64                 # Helmoltz energy (F0, J/mol)
    n::Float64                  # negative of the number of atoms per formula unit (-n)
    V0::Float64                 # negative of the volume (-V0)
    K0::Float64                 # c1: isothermal bulk modulus (K0, bar)
    Kp::Float64                 # c2: pressure derivative of the isothermal bulk modulus (K')
    Θ0::Float64                 # c3: Debye Temperature (Θ0, K)
    γ0::Float64                 # c4: Gruneisen thermal parameter (γ0)
    q0::Float64                 # c5: Mie-Gruneisen exponent (q0)
    ηS0::Float64                # c6: Shear strain derivative of the tensorial Gruneisen parameter (ηS0)
    cme::Float64                # c7: Configurational (and magnetic) entropy (J/mol/K)
end

struct ModelJSON
    name::String
    species::Vector{String}
    margules::Dict{String, Float64}
    sites::Int64
end

struct Model
    name::String
    species::DataFrame
    margules::Dict{String, Float64}
    sites::Int64
end

function read_data(fname::String)
    return JSON3.read(fname, Vector{Phase}) |> DataFrame
end

function read_models(fname::String, data::DataFrame, model_names::Vector{String})
    read_models = JSON3.read(fname, Vector{ModelJSON})
    models = Vector{Model}()
    
    for model in read_models
        aux_data = DataFrame()
        if model.name in model_names
            for i in 1:length(model.species)
                p = findfirst(x -> x == model.species[i], data.id)
                push!(aux_data, data[p, :])
            end
            push!(models, Model(model.name, aux_data, model.margules, model.sites))
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
        println("Starting...")
    elseif str == "PT"
        @printf(" * P: %.2f bar", arg[1])
        @printf(", T: %.2f K\n", arg[2])
    elseif str == "gibbs"
        @printf(" * gibbs: \t%.2f\n", arg[1])
    elseif str == "activity"
        @printf(" * activity: \t%.2f\n", arg[1])
    elseif str == "excess"
        @printf(" * excess: \t%.2f\n", arg[1])
    elseif str == "μ"
        @printf(" * μ: \t\t%.2f\n", arg[1])
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
    # println("F: ", a)
    G = a + p * v - t * phase.cme
    return G
end

function calc_activity(phase::DataFrameRow{DataFrame, DataFrames.Index}, comp::Vector{Float64}, model::Model, species_fractions::Vector{Float64})
    # initialization
    act = 0.0
    n_species = size(model.species)[1]
    n_sites = model.sites
    n_comps = length(model.species.cmp[1])

    aux = [1.0/8.0, 1.0/4.0]

    for k in 1:n_sites

        sijk1 = [value for value in values(phase.cmp)]
        sijk1 .*= aux[k]

        Sik = sum(sijk1)

        Nk = 0.0
        for c in 1:n_comps
            Njk = 0.0
            for i in 1:n_species
                sijk2 = [value for value in values(model.species.cmp[i])] 
                sijk2 .*= aux[k]
                Njk += sijk2[c] .* species_fractions[i]
            end
            Nk += Njk
        end

        a1 = Sik * log(Nk)

        a2 = 0.0
        for c in 1:n_comps
            Njk = 0.0
            for i in 1:n_species
                sijk2 = [value for value in values(model.species.cmp[i])] 
                sijk2 .*= aux[k]
                Njk += sijk2[c] .* species_fractions[i]
            end
            a2 += (Njk != 0) ? sijk1[c] * log(Njk) : 0
        end
        # println(a1, ", ", a2, " = ", a1 -a2)
        act += (a1 - a2)
    end

    return act
end

function eye(i::Int64, j::Int64)
    return i == j ? 1.0 : 0.0
end

function calc_excess(i::Int64, phase::DataFrameRow{DataFrame, DataFrames.Index}, model::Model, species_fractions::Vector{Float64})
    
    excess = 0.0
    n_species = size(model.species)[1]
    W = [value for value in values(model.margules)]
    # println(model.margules)

    # sum_v = 0.0
    # for i in 1:n_species
    #     sum_v += species_fractions[i] * v[i]
    # end

    # for i in 1:n_species
    #     mat_phi[i] = (species_fractions[i] * v[i]) / sum_v
    # end

    # for i in 1:n_species
    excess = 0.0
    it = 1
    for j in 1:n_species-1
        for k in j+1:n_species
            excess += (eye(i,j) - species_fractions[j]) * (eye(i,k) - species_fractions[k]) * W[it]
            # excess -= (eye(i,j) - mat_phi[j]) * (eye(i,k) - mat_phi[k]) * (W[it] * 2.0 * v[i] / (v[j] + v[k]))
            it += 1
        end
    end
    # end

    return excess
end

function gcalc(pressure::Float64, temperature::Float64, comp::Vector{Float64}, models::Vector{Model}, species_fractions::Vector{Vector{Float64}})
    
    μ = Vector{Float64}()

    for (m, model) in enumerate(models)
        gi = Vector{Float64}() 
        ai = Vector{Float64}()    
        xi = Vector{Float64}()
        μi = Vector{Float64}()
        n_species = size(model.species)[1]
        for i in 1:n_species
            phase = model.species[i, :]
            title = " " * string(species_fractions[m][i] * 100.0) * " % of `" * phase.id * "` [" * phase.fml * "] "
            message(title);
            g = calc_gibbs(phase, pressure, temperature)
            message("PT", [pressure, temperature])
            message("gibbs", [g])
            push!(gi, g)
            a = R * temperature * calc_activity(phase, comp, model, species_fractions[m]) 
            message("activity", [a])
            push!(ai, a)
            e = calc_excess(i, phase, model, species_fractions[m])
            message("excess", [e])
            push!(xi, e)
            message("μ", [g - a - e])
            push!(μi, (g - a - e))
            
        end
        
        push!(μ, sum(μi .* species_fractions[m]))
        message("line")
        @printf(" * μ : \t\t%10.2f\n", sum(μ))
        message("line")
    end

    return μ
end

function span_gcalc(n_span::Int64, pressure::Float64, temperature::Float64, comp::Vector{Float64}, models::Vector{Model})
    
    g = Vector{Float64}()
    x = [x for x in 0:n_span] ./ n_span

    for a in 0:n_span
        v = a / n_span
        species_fractions = [[v, 1-v]]
        push!(g, sum(gcalc(pressure, temperature, comp, models, species_fractions)))
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