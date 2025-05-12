
function prova(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3}) 

neq=0
# equation numbering 
for e in 1:info.NE           
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      for d in 1:2
        if nodes[n].dof[d]==0   
          neq=neq+1             
          nodes[n].dof[d]=neq   
        end  
      end
    end
  end   
end

info.neq=neq


###################################
# follower forces
###################################

# allocation of elements arrays
Xe=zeros(Float64,3,2)
dofUe=zeros(Int64,6,1)

FeT=zeros(Float64,6,1)
Fe0=zeros(Float64,6,1)
Fe1=zeros(Float64,6,1)
Fe2=zeros(Float64,6,1)
KUUe=zeros(Float64,6,6)
KUUe2=zeros(Float64,6,6)

#for e in 1:info.NL                     

e=1

Ue=rand(Float64,3,2)
press=rand(Float64)

for k in 1:3                                 
  n=B3[e].nodes[k]
  Xe[k,:]=nodes[n].coor              
  dofUe[(k-1)*2+1:k*2]=nodes[n].dof
  Ue[k,:]=nodes[n].u
end

B3_fl!(Fe0,KUUe,Xe,Ue,press) 

DUe=rand(Float64,3,2)
Dpress=rand(Float64)

B3_fl!(FeT,KUUe2,Xe,Ue+DUe,press+Dpress) 
println("T ",FeT)

println("0 ",Fe0)
vUe=[DUe[1,:]' DUe[2,:]' DUe[3,:]']'
#println("0 ",vUe)
#println(size(vUe))
Fe1=-KUUe*vUe+Fe0/press*Dpress
println("1 ",Fe1)
#println("T ",Fe0+Fe1)

B3R_quad_fl!(Fe2,Xe,Dpress,Dpress,DUe,DUe) 
println("2 ",Fe2)
println("T ",Fe0+Fe1+Fe2)

#end

end   



function B3R_quad_fl!(Fe::Matrix{Float64},Xe::Matrix{Float64},
                     p1e::Float64,p2e::Float64,U1e::Matrix{Float64},U2e::Matrix{Float64})

fill!(Fe,0.0)

qr=quadrature_points(Val(:LINE3gp))
for (w,a) in qr
  D=[a-.5  a+.5 -2*a]                
  FU1=D*U1e;
  JNU1=[-FU1[2], FU1[1]]                 
  FU2=D*U2e;
  JNU2=[-FU2[2], FU2[1]]                 
  NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2]   # shape functions
  N=[NL[1] 0 NL[2] 0 NL[3] 0;
     0 NL[1] 0 NL[2] 0 NL[3]]   
  Fe[:]+=transpose(N)*(JNU1*p2e+JNU2*p1e)*w/2  # quadratic
end

end


