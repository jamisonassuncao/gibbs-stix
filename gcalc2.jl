global R = 8.31446261815324

struct Phase
    id::String      # Name id
    F0::Float64     # Helmoltz energy (F0, J/mol)
    n::Float64      # negative of the number of atoms per formula unit (-n)
    V0::Float64     # negative of the volume (-V0)
    K0::Float64     # c1: isothermal bulk modulus (K0, bar)
    Kp::Float64     # c2: pressure derivative of the isothermal bulk modulus (K')
    Θ0::Float64     # c3: Debye Temperature (Θ0, K)
    γ0::Float64     # c4: Gruneisen thermal parameter (γ0)
    q0::Float64     # c5: Mie-Gruneisen exponent (q0)
    ηS0::Float64    # c6: Shear strain derivative of the tensorial Gruneisen parameter (ηS0)
    cme::Float64    # c7: Configurational (and magnetic) entropy (J/mol/K)
end
# cp(pr,t) = c1 + c2*t + c3/(t*t) + c4*t*t + c5/t**(1/2) + c6/t + c7 /t**3

seif    = Phase("seif", -794335.4, -3, -1.367000, 3275843.0, 4.01553, 1140.772, 1.3746600, 2.83517, 4.971078, 0.0)
coe     = Phase( "coe", -855068.5, -3, -2.065700, 1135856.0, 4.00000, 852.4267, 0.3915700, 1.00000, 2.397930, 0.0)
q       = Phase(   "q", -858853.4, -3, -2.367003, 495474.30, 4.33155, 816.3307, -0.296e-2, 1.00000, 2.364690, 0.0)
st      = Phase(  "st", -818984.6, -3, -1.401700, 3143352.0, 3.75122, 1107.824, 1.3746600, 2.83517, 4.609040, 0.0)

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
    gcalc(t, p, phase)

This function calculates the Gibbs energy of the `phase`

# Arguments
- `t::Float64`: Temperature value in K.
- `p::Float64`: Pressure value in bar.
- `phase::Phase`: Phase object containing the phase data.

# Returns
- `G::Float64`: Gibbs energy value.
"""
function gcalc(t=1000.0, p=1000.0, phase=q)
    println("Calculating Gibbs energy for `", phase.id, "`")
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
    # # c                                 destabilize the phase.
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
    println("G: ", G)
    println(repeat("=", 40))
    return G
end

function Σαμ()

end

G = gcalc(1000.0, 1000.0, seif);
G = gcalc(1000.0, 1000.0, coe);
G = gcalc(1000.0, 1000.0, q);
G = gcalc(1000.0, 1000.0, st);