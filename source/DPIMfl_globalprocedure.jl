function globalprocedure(input_file,output_folder)


fname = "../input/"*input_file
include(fname)
info.mesh_file = file

T6, B3, nodes = readgmsh!()

##############################
# static solution and matrix assemblage
##############################

#prova(nodes,T6,B3)

K, M, F0 = analysis(nodes,T6,B3)
if info.beta != 0
    K0 = damping_stiffness(nodes,T6,B3) # For stiffness damping
else
    K0 = nothing
end

##############################
# eigenvalue computation
##############################

Lambda, VR, VL, VRp = eigAB(K,M,F0,nodes,K0)

##############################
# launch DPIM
##############################

P = dpim(K,M,Lambda,VR,VL,VRp,nodes,T6,B3,K0) # computes parametrization
realification!(P)  # performs realification

##################################################
#  OUTPUT ON FILE
##################################################

isdir("./" * output_folder) || mkdir("./" * output_folder)
output(P, output_folder)  # output

##################################################
#  output for matcont
##################################################

matcont(nodes,M,VR,P,output_folder)

return P
end