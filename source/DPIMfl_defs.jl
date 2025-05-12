
#const eps0 = 1

mutable struct ST6
  mat::Int64
  nodes::Vector{Int64}
end

mutable struct SB3
  nodes::Vector{Int64}
  press::Float64
  press_damping::Float64
end

# define a new type
mutable struct Snode
 coor::Vector{Float64}
 dof::Vector{Int64}
 u::Vector{Float64}
end

mutable struct Sinfo
 mesh_file::String 
 NN::Int64
 NE::Int64
 NL::Int64
 neq::Int64
 #nK::Int64
 nA::Int64
 nMat::Int64
 #
 nmm::Int64    # number of master modes
 Lmm::Vector{Int64}  # list of mm
 LambdaSym::Vector{Int64}  # list of mm
 nza::Int64    # autonomous
 #nzna::Int64   # nonautonomous
 #nrom:: Int64
 Ffreq::Int64
 #Fmodes::Vector{Int64}
 #Fmult::Vector{Float64}
 alpha::Float64
 beta::Float64
 tol::Float64
 neig::Int64   # number of computed modes
 style::Char
 max_order::Int64
 max_orderNA::Int64
 Jordan::Bool # Whether or not to use Jordan blocks
 tolHopf::Float64
 symbolic_Jordan::Bool
 Hopf_pairs::Matrix{Int64}

 Sinfo() = new()
end

