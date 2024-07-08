T::Float64 = 1000.0 # K
R::Float64 = -8.31446261815324

m::Vector{Float64} = [1.0, 1.0, 1.0, 1.0, 1.0] # multiplicity of each element
C::Vector{Vector{Float64}} = [[1, 0, 0, 0, 1], [0, 1, 0, 0, 1], [0, 0, 1, 1, 0]]

x = m * C'