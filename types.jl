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
    endmembers::Vector{String}
    margules::Dict{String, Float64}
    sites::Int64
end

struct Model
    name::String
    endmembers::DataFrame
    margules::Dict{String, Float64}
    sites::Int64
end