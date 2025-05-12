
#
# sets input
#

input_file="cantilever.jl"

using Arpack
using ExtendableSparse
using LinearAlgebra
using MAT
using MATLAB
using Pardiso
using SparseArrays

include("DPIMfl_defs.jl")
include("DPIMfl_globalprocedure.jl")
include("DPIMfl_readgmsh.jl")
include("DPIMfl_outgmsh.jl")
include("DPIMfl_elemental.jl")
include("DPIMfl_analysis.jl")
include("DPIMfl_prova.jl")
include("DPIMfl_newmark.jl")
include("DPIMfl_eig.jl")
include("DPIMfl_param_struct.jl")
include("DPIMfl_dpim.jl")
include("DPIMfl_dpim_FEM.jl")
include("DPIMfl_dpim_realification.jl")
include("DPIMfl_dpim_output.jl")
include("quadrature.jl")
include("shape_functions.jl")
include("DPIMfl_matcont.jl")
include("DPIMfl_eigenvalue_trajectories.jl")

output_folder = "output/output_2modes_o9_after_far_MassDamped_graph"
P = globalprocedure(input_file, output_folder);

# Time-history analysis
# newmark(input_file)

# Eigenvalue trajectories
# P_min = 0
# P_max = 20
# n_P = 1000
# output_folder = "output/eigenvalues_StiffnessDamped"
# eigenvalue_trajectories(input_file, output_folder, P_min, P_max, n_P)

println("done!")




