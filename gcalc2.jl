seif = [-794335.4000000000,    # Helmoltz energy (F0, J/mol)
        -3.000000000000000,    # negative of the number of atoms per formula unit (-n)
        -1.367000000000000,    # negative of the volume (-V0)
        3275843.0000000000,    # c1: isothermal bulk modulus (K0, bar)
        4.0155300000000000,    # c2: pressure derivative of the isothermal bulk modulus (K')
        1140.7720000000000,    # c3: Debye Temperature (Θ0, K)
        1.3746600000000000,    # c4: Gruneisen thermal parameter (γ0)
        2.8351700000000000,    # c5: Mie-Gruneisen exponent (q0)
        4.9710780000000000,    # c6: Shear strain derivative of the tensorial Gruneisen parameter (ηS0)
        0.0000000000000000]    # c7: Configurational (and magnetic) entropy (J/mol/K)

coe = [-855068.5, -3, -2.065700, 1135856.0, 4.00000, 852.4267, 0.3915700, 1.00000, 2.39793, 0.0]
q   = [-858853.4, -3, -2.367003, 495474.30, 4.33155, 816.3307, -0.296e-2, 1.00000, 2.36469, 0.0]
st  = [-818984.6, -3, -1.401700, 3143352.0, 3.75122, 1107.824, 1.3746600, 2.83517, 4.60904, 0.0]
# cp(pr,t) = c1 + c2*t + c3/(t*t) + c4*t*t + c5/t**(1/2) + c6/t + c7 /t**3
global R = 8.31446261815324
println("Starting...")

function plg(t)
    # -----------------------------------------------------------------------
    # evaluates debye integral: int((ln(1-exp(-t))*t^2),t=0..t)
    # -----------------------------------------------------------------------
    p0 = exp(-t)
    p1 = 1.0
    p2 = t*t
    p3 = 2.0*t
    nopt50 = 4.7637509757004091e-015
    # c 45/Pi
    plg = -2.1646464674222763831

    i = 1
    while i < 100000

        p4 = Float64(i)
        p1 = p1 * p0
        dinc = (p2 + (p3 + 2.0/p4)/p4)*p1/p4/p4
        plg = plg + dinc

        if (abs(dinc/(1.0+abs(plg))) < nopt50) 
            return plg
        end

        i += 1
    end 
    return plg
end

function gcalc(t=1000.0, p=1000.0, thermo=seif)

    tr = 300.0                                          # K

    v0     = -thermo[3]
    nr9    = -9.0*thermo[2]*R                           #thermo[11]
    nr9t0  = nr9*tr                                     #thermo[20]
    c1     = -9.0*thermo[3]*thermo[4]                   #thermo[12]
    c2     = thermo[5]/2.0-2.0                          #thermo[13]
    c3     = 3.0*c1*c2                                  #thermo[14]
    aii    = 6.0*thermo[7]                              #thermo[15]
    aiikk2 = 0.5*aii*(-2.0+6.0*thermo[7]-3.0*thermo[8]) #thermo[18]
    aii2   = 3*thermo[7]                                #thermo[19]

    # aiikk  = thermo[16]

    r23 = 2.0/3.0
    r59 = 5.0/9.0

    iopt21 = 100
    nopt51 = 1.8590370495272219e-13

    t1     = thermo[6]/t
    t2     = t/tr
    nr9t   = nr9*t
    # initial guess for volume:
    #   taylor(diff(FTH,v),v=v0,1)
    #   JADC Feb 26, 2008.
    #   the dfth0 could be loaded as a constant.
    tht    = t1
    tht0   = tht*t2

    k00    = thermo[4]
    k0p    = thermo[5]
    gamma0 = thermo[7]

    dfth   = nr9t*gamma0/v0*(3.0*plg(tht)/tht^3 - log(1.0-exp(-tht)))
    dfth0  = nr9t0*gamma0/v0*(3.0*plg(tht0)/tht0^3 - log(1.0-exp(-tht0)))
    #  taylor(diff(FC,v),v=v0,2)
    #       v = (k00-dfth+dfth0-p)/k00*v0
    #  taylor(diff(FC,v),v=v0,3)
    root = k00*((2.0+2.0*k0p)*(p+dfth-dfth0) + k00)

    v = 0.0
    if (root > 0.0) 
        v = ((2.0+k0p)-sqrt(root)/k00)*v0/(1.0+k0p)
        if (v < v0/10.0) || (v > v0*10.0) v = v0 end
    else
        v = v0
    end 

    itic = 0
    ibad = 4
    bad = true

    while (itic < 100)

        itic += 1
        println("itic: ", itic)
        #  f, and derivatives
        v23 = (v0/v)^r23
        f = 0.5*v23 - 0.5
        df = -v23/v/3.0
        d2f = r59*v23/v^2
        # cold part derivatives
        dfc = (c3*f+c1)*f*df
        d2fc = (2.0*c3*f+c1)*df^2+(c3*f+c1)*f*d2f
        # debye T/T (tht)
        z = 1.0+(aii+aiikk2*f)*f

        if (z < 0.0 || v/v0 > 100.0 || v/v0 < 1e-2) println("ERROR z or v/v0") end
        # println(z)

        root = sqrt(z)

        tht   = t1*root
        tht0  =  tht*t/tr
        # tht derivatives
        a2f   = aii2+aiikk2*f
        da    = a2f/root
        dtht  = t1*da*df
        d2tht = t1*((aiikk2/root-a2f^2/z^1.5)*df^2 + da*d2f)
        # println(d2tht)

        dtht0 = dtht*t2
        d2tht0 = d2tht*t2
        # polylog functions:
        fpoly   = 3.0*plg(tht)/tht^3
        fpoly0  = 3.0*plg(tht0)/tht0^3
        # thermal part derivatives:
        etht  = exp(-tht )

        if (1.0-etht < 0.0) println("ERROR 1-etht") end

        letht = log(1.0-etht)

        dfth = (letht-fpoly)*nr9t*dtht/tht
        d2fth = ((4.0*dtht^2/tht-d2tht)*(fpoly-letht) + dtht^2*etht/(1.0-etht))*nr9t/tht

        etht0 = exp(-tht0)

        if (1.0-etht0 < 0.0) println("ERROR 1-tht0") end

        letht0 = log(1.0-etht0)

        dfth0 = (letht0-fpoly0)*nr9t0*dtht0/tht0
        d2fth0 = ((4.0*dtht0^2/tht0-d2tht0)*(fpoly0-letht0) + dtht0^2*etht0/(1.0-etht0))*nr9t0/tht0

        f1  = -dfc - dfth + dfth0 - p
        df1 = -d2fc - d2fth + d2fth0
        dv = f1/df1

        if (v - dv < 0.0) dv = v/2.0 end

        # println("v: ", v)
        v -= dv

        if (itic > iopt21 || abs(f1) > 1e40) 
        # allow bad result
            if (abs(f1/p) < 0.0) 
                ibad = 5
                println("ERROR abs(f1/p)")
            end
        elseif (abs(dv/(1.0+v)) < nopt51)
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
    # everything is ok, now get helmoltz energy:
    f = 0.5*(v0/v)^r23 - 0.5
    # println(f)
    z = 1.0+(aii+aiikk2*f)*f
    root = sqrt(z)
    # final estimate for tht
    tht   = t1*root
    tht0  = tht*t2
    # helmholtz enery
    a = thermo[1] + c1*f^2*(0.5 + c2*f) + nr9*(t/tht^3*plg(tht) -tr/tht0^3*plg(tht0))
    println("ℱ: ", a)
    gstxgi = a + p*v - t*thermo[10]
    println("𝒢: ", gstxgi)
    # # z = (theta/theta0)^2
    # gamma = (2.0*f+1.0)*(aii+aiikk*f)/6.0/z
    # etas = - gamma - thermo[17]/z*(2.0*f + 1.0)^2
    # # thermal energy/V, based on previous v estimate
    # if (gamma < 0.0)
    #     ethv = (dfth0-dfth)/gamma
    # else
    #     ethv = 0.0
    # end
    # # adiabatic shear modulus
    # smu = (1.0 + 2.0*f)^(2.5)*(emod[1] + f*(thermo[21] + thermo[22]*f)) - etas*ethv
end

gcalc(1000.0, 1000.0, seif)
gcalc(1000.0, 1000.0, coe)
gcalc(1000.0, 1000.0, q)
gcalc(1000.0, 1000.0, st)
println("Finished!")