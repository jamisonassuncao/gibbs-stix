function plg(t)
    # 29th order series expansion of polylog terms in sixtrude's EoS about
    # t = 0, good to 0.0001 accuracy to t = 5.
    # f := int((ln(1-exp(-t))*t^2),t=0..TD);
    # gg := convert(series(f,TD=0,29),polynom);
    # see equation 4 in Stixrude 1990

    c2 = 0.1111111111111111111

    plg = t^3.0 * ((log(t) / 3.0 - c2) - t / 8.0 + t^2.0 / 120.0 \
                    -t^4 / 20160.0 + t^6 / 1632960.0 - t^8 / 106444800.0 \
                    +t^10 / 6227020800.0 - 0.2935661189E-11 * t^12 \
                    +0.562291451E-13 * t^14 - 0.1115026413E-14 * t^16 \
                    +0.2271444989E-16 * t^18 - 0.4727975432E-18 * t^20 \
                    +0.1001636878E-19 * t^22 - 0.2153466772E-21 * t^24)

    return plg
end

function print1000(v1, v2, str1)
    println("**warning** failed to converge at ")
    println("T= ", round(v1, digits=2), " K ")
    println("P= ", round(v2, digits=1), " bar")
    println("Using Stixrude EoS.")
    println("Phase ", str1, " will be destabilized.")
end

function print1010()
    println("#" ^ 80)
    println("# **warnign** will not be repeated for future ")
    println("# instances of this problem")
    println("#" ^ 80)
end

function check_error(IZAP, NAME, T, P)
    if IZAP < 100
        println("**error**")
        print1000(T, P, NAME)
    elseif IZAP == 100 
        print1010()
    end 
    IZAP += 1
    𝒢 = 1.0e30
    return IZAP, 𝒢
end

function gcalc(T=1000.0, P=1000.0, 
    # T0=300.0,
    # H0R=-860.11803, 
    # NN0=1.0, 
    # VAA=5.76, 
    # K0=614.2537, 
    # K01=197.8011, 
    # V0R=22.421050, 
    # γ0=-0.03958, 
    # q0=1.0,
    # θ0=884.20481)
    T0=300.0,
    H0R=-794335.40, 
    NN0=3.0, 
    VAA=4.97107, 
    K0=3275843.0, 
    K01=4.015530, 
    V0R=1.367, 
    γ0=-1.374660, 
    q0=2.835170,
    θ0=1140.77199)
    """
    Volume function from Stixrude and Lithgow-Bertelloni 2011

    # Arguments                                                     -> .f90 arguments
    'H0R': Helmoltz energy (F0, J/mol)                              -> H0R
    'NN0': number of atoms per formula unit (n)                     -> NAT & NN0
    'VAA': Configurations (and magnetic) entropy (J/mol/K)          -> VAA
    'K0' : isothermal bulk modulus (K0, kbar)                       -> VPA & K0
    'K01': pressure derivative of the isothermal bulk modulus (K')  -> VPB & K01
    'V0R': reference volume (V0)                                    -> V0R
    'γ0' : Gruneisen thermal paramenter (γ0)                        -> D2
    'q0' : Mie-Gruneisen exponent (q0)                              -> D3
    'θ0' : Debye Temperature (θ0, K)                                -> D1
    """

    # The value of F0 is:   -794335.40000000002     
    # The value of n is:   -3.0000000000000000     
    # The value of v0 is:   -1.3670000000000000     
    # The value of k0 is:    3275843.0000000000     
    # The value of k0p is:    4.0155300000000000     
    # The value of theta0 is:    1140.7719999999999     
    # The value of gamma0 is:    1.3746600000000000     
    # The value of q is:    2.8351700000000002     
    # The value of etaS0 is:    4.9710799999999997     
    # The value of Smag is:    0.0000000000000000

    IZAP = 0

    # initial parameters
    R = 8.31446261815324 # Gas constant (J/(mol K))
    NR9 = 9.0 * NN0 * R # 9nR
    C1 = 9.0 * K0 * V0R # 9K₀V₀
    C2 = K01 / 2.0 - 2.0 # K'/2 - 2
    C3 = 3.0 * C1 * C2
    C5 = 3.0 * γ0 # 3γ₀
    C7 = C5 * (-2.0 + 6.0 * γ0 - 3.0 * q0)
    QM1 = q0 - 1.0
    DT = T - T0

    # 'Newton–Raphson'-method approximation:
    #   initial guess for volume uses a 4 order expansion on the polylog terms 
    #   before differentiation and a second order series on the resulting 
    #   derivative.
    VOLUM = V0R * (3.0 * ((2.0 - q0) * DT * γ0 * NR9 - 3.0 * P * V0R) + C1) / (C1 - QM1 * DT * γ0 * NR9)
    # VOLUM = V0R * (3.0*((2.0-q0)*DT*γ0*NR9 - 3.0*P*V0R) + C1) / 
    # (C1 - QM1 * DT * γ0 * NR9)
    if (VOLUM < 0.0) VOLUM = V0R end

    # run 'Newton-Raphson'-method approximation
    ITIC = 0
    FL = 1.0e9 # (bar)
    TOL = 1.0e-6 * P # (bar)

    ROOT_ER = 0
    ETHT_ER = 0
    VOLU_ER = 0
    ITIC_ER = 0
    FL_ER   = 0

    while (abs(FL) > TOL)
        ITIC += 1
        VQ = (VOLUM / V0R)^q0
        GVQ = γ0 * VQ
        # Birch-Murnaghan finite strain
        V23 = (V0R / VOLUM)^(2.0/3.0)           # -> eq. 16 @ https://www.perplex.ethz.ch/perplex_thermodynamic_data_file.html
        F = 0.5 * V23 - 0.5                     # -> eq. 16 @ https://www.perplex.ethz.ch/perplex_thermodynamic_data_file.html
        DF = -V23 / VOLUM / 3.0
        D2F = 5.0 / 9.0 * V23 / VOLUM^2.0
        # Debye T/T 
        ROOT = sqrt(1.0 + C7 * F^2.0 + 2.0 * C5 * F) # eq. 21 @ https://www.perplex.ethz.ch/perplex_thermodynamic_data_file.html
        if (ROOT^2.0 < 0.0)
            ROOT_ER += 1
            IZAP, 𝒢 = check_error(IZAP, :ROOT, T, P)
        end

        THT = θ0 / T * ROOT
        THT0 = THT * T / T0

        ETHT = exp(-THT)
        if (1-ETHT < 0.0)
            ETHT_ER += 1
            IZAP, 𝒢 = check_error(IZAP, :ETHT, T, P)
        end
        ETHT0 = exp(-THT0)
        LTHT = log(1.0 - ETHT)
        LTHT0 = log(1.0 - ETHT0)
        # diff(theta/T,v)
        DD2 = (C7 * F + C5) / ROOT
        DD1 = θ0 * DF * DD2
        DTHT = DD1 / T
        DTHT0 = DD1 / T0
        # diff(theta/T,v,v)
        DD3 = θ0 * (DF^2.0 * (C7 - C5^2.0) / ROOT^3.0 + D2F * DD2)
        D2THT = DD3 / T
        D2THT0 = DD3 / T0

        PLGG = plg(THT)
        PLGG0 = plg(THT0)
        DG = THT^2.0 * LTHT * DTHT
        DG0 = THT0^2.0 * LTHT0 * DTHT0
        D2G = ((2.0 * LTHT + THT * ETHT / (1.0 - ETHT)) * DTHT^2 \
               +THT * LTHT * D2THT) * THT
        D2G0 = ((2.0 * LTHT0 + THT0 * ETHT0 / (1.0 - ETHT0)) * DTHT0^2 \
                +THT0 * LTHT0 * D2THT0) * THT0

        DFC = (C3 * F + C1) * F * DF
        D2FC = (2.0 * C3 * F + C1) * DF^2 + (C3 * F + C1) * F * D2F

        DFT = NR9 * T / THT^3 * (DG - 3.0 / THT * PLGG * DTHT)
        DFT0 = NR9 * T0 / THT0^3 * (DG0 - 3.0 / THT0 * PLGG0 * DTHT0)

        D2FT = NR9 * T / THT^3 * (3.0 / THT * (DTHT * (4.0 / THT * PLGG * DTHT \
                                                       -2.0 * DG) - PLGG * D2THT) + D2G)
        D2FT0 = NR9 * T0 / THT0^3 * (3.0 / THT0 * (DTHT0 * (4.0 / THT0 * PLGG0 * DTHT0 \
                                                            -2.0 * DG0) - PLGG0 * D2THT0) + D2G0)

        FL  = - DFC - DFT + DFT0 - P   # P(V) bar
        # println(ITIC, ": ", FL)

        DFL = - D2FC - D2FT + D2FT0    # dP(V)/dV
        #'Newton-Raphson'-method
        VOLUM -= FL / DFL
        if (VOLUM <= 0.0 || ITIC > 40 || abs(FL) > 1.0e40)
            if (VOLUM <= 0) VOLU_ER += 1 end
            if (ITIC > 0) ITIC_ER += 1 end
            if (FL <=0) FL_ER += 1 end
            IZAP, 𝒢 = check_error(IZAP, :VOLUM, T, P)
        end
    end

    VQ = (VOLUM / V0R)^q0
    # Volume (Birch-Murnaghan finite strain)
    F = 0.5 * (V0R / VOLUM)^(2.0/3.0) - 0.5
    ROOT = sqrt(1.0 + C7 * F^2.0 + 2.0 * C5 * F)
    THT = θ0 / T * ROOT
    THT0 = THT * T / T0

    # Helmholtz energy
    ℱ = H0R + C1 * F^2.0 * (0.5 + C2 * F) + NR9 * (T / THT^3 * plg(THT) - T0 / THT0^3.0 * plg(THT0))
    # ℱ = H0R + 
    # C1 * F^2.0 * (0.5 + C2 * F) + 
    # NR9 * (T / THT^3 * plg(THT) - 
    # T0 / THT0^3.0 * plg(THT0))
    # Legendre transformation
    𝒢 = ℱ + P * VOLUM - T * VAA
    println("# ROOT errors: ", ROOT_ER)
    println("# ETHT errors: ", ETHT_ER)
    println("# VOLU errors: ", VOLU_ER)
    println("# ITIC errors: ", ITIC_ER)
    println("# FL errors: ", FL_ER)

    println("ℱ: ", ℱ)
    println("𝒢: ", 𝒢)

    return 𝒢
end

println("Starting...")
𝒢 = gcalc()

println("Done!")

# if (VO2) 
#     FVVOL = V0R * exp(VAA * (T - T0) + VAB * (TT - TT0) / 2.0)
#     if (VB > 0.0) 
#         FGVOL=(FVVOL / VB) * (1.0 - exp(- VB * (P - P0)))
#         GR = GR + FGVOL
#         VOLUM = FVVOL * DEXP(- VB * (P - P0))
#         FVVOL = VOLUM - V0R
#     else
#         FGVOL = FVVOL * (P - P0)
#         GR = GR + FGVOL
#         VOLUM = FVVOL
#         FVVOL = VOLUM - V0R
#     end
# end