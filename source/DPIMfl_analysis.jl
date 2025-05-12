
function analysis(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3}) 

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

sol=zeros(Float64,neq)
F=zeros(Float64,neq)
Kfl=zeros(Float64,neq)

residuum=100
iter=0
while residuum>1e-10
#while iter<1

iter+=1

global K   # to make it visible out of the while loop

# faster and safer to reallocate at every iter, as the number of nozero coeffs may vary from iter to iter
K=SparseMatrixLNK(neq,neq)
if iter >0 
  fill!(F,0.0)
  fill!(Kfl,0.0)
end

###################################
# mechanics
###################################

println("Assembling mechanics")

# allocation of elements arrays
Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12)
Ue=zeros(Float64,12)
Fe=zeros(Float64,12)
Ke= zeros(Float64,(12,12))

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
      Ue[(k-1)*2+1:k*2]=nodes[n].u
    end
#    T6_Ke!(Ke,Xe,material[mat])
    T6_KeNL!(Ke,Fe,Xe,Ue,material[mat])
    for i = 1:12
      dofi=dofe[i]
      if dofi>0
        F[dofi]+=Fe[i]
        for j = 1:12
#          F[dofi]-=Ke[i,j]*Ue[j]
          dofj=dofe[j]
          if dofj>0 
            K[dofi,dofj]+=Ke[i,j]
          end
        end
      end
    end
  end
end

###################################
# follower forces
###################################

# allocation of elements arrays
Xe=zeros(Float64,3,2)
dofUe=zeros(Int64,6,1)
Ue=zeros(Float64,3,2)

Fe=zeros(Float64,6,1)
KUUe=zeros(Float64,6,6)

for e in 1:info.NL                     
  for k in 1:3                                 
    n=B3[e].nodes[k]
    Xe[k,:]=nodes[n].coor              
    dofUe[(k-1)*2+1:k*2]=nodes[n].dof
    Ue[k,:]=nodes[n].u
  end

  press=B3[e].press
  B3_fl!(Fe,KUUe,Xe,Ue,press) 

  for i in 1:6
    dofUi=dofUe[i]
    if dofUi>0
      F[dofUi]+=Fe[i]
      Kfl[dofUi]+=Fe[i]/press   # for press
      
      for j in 1:6
        dofUj=dofUe[j]
        if dofUj>0
          K[dofUi,dofUj]+=KUUe[i,j]  
        end
      end
    end
  end    

end


# Solution phase
println("Solving system")
K=SparseMatrixCSC(K)  # convert to other sparse format
sol = K\F
#ps = MKLPardisoSolver()
#solve!(ps,sol,K,F)

# fill dofs
for n in 1:info.NN 
  for d in 1:2       
    dofU=nodes[n].dof[d]
    if dofU>0                      
      nodes[n].u[d]+=sol[dofU]      
    end 
  end
end
 
residuum=norm(sol)
println("Residuum: ",residuum)

end # while residuum

##############################
# mass assembler
##############################

M=SparseMatrixLNK(neq,neq)

# allocation of elements arrays
Xe=zeros(Float64,(6,2))
dofe=zeros(Int64,12,1)
Me= zeros(Float64,(12,12))

for e in 1:info.NE     
  mat=T6[e].mat
  if mat>0 
    for k in 1:6
      n=T6[e].nodes[k]   
      Xe[k,:]=nodes[n].coor
      dofe[(k-1)*2+1:k*2]=nodes[n].dof
    end
    T6_Me!(Me,Xe,material[mat])
    for i = 1:12
      dofi=dofe[i]
      if dofi>0
        for j = 1:12
          dofj=dofe[j]
          if dofj>0 
            M[dofi,dofj]+=Me[i,j]
          end
        end
      end
    end
  end
end

M=SparseMatrixCSC(M)

##############################
# GMSH output of static solution
##############################

outgmsh(info.NN,nodes,T6,B3,"_u0")
#  mycommand = `../gmsh.exe post.msh`
#  run(mycommand)

return K,M,Kfl

end   


# M Utt + K U + Kfl p

function damping_stiffness(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3}) 

  neq = info.neq
  
  sol=zeros(Float64,neq)
  F=zeros(Float64,neq)
  U = zeros(Float64,neq)
  
  residuum=100
  iter=0
  while residuum>1e-10
  #while iter<1
  
  iter+=1
  
  global K   # to make it visible out of the while loop
  
  # faster and safer to reallocate at every iter, as the number of nozero coeffs may vary from iter to iter
  K=SparseMatrixLNK(neq,neq)
  if iter >0 
    fill!(F,0.0)
  end
  
  ###################################
  # mechanics
  ###################################
  
  println("Assembling mechanics")
  
  # allocation of elements arrays
  Xe=zeros(Float64,6,2)
  dofe=zeros(Int64,12)
  Fe=zeros(Float64,12)
  Ke= zeros(Float64,(12,12))
  
  for e in 1:info.NE
    Ue=zeros(Float64,12)   
    mat=T6[e].mat
    if mat>0 
      for k in 1:6
        n=T6[e].nodes[k]   
        Xe[k,:]=nodes[n].coor
        dofe[(k-1)*2+1:k*2]=nodes[n].dof
        for d in 1:2       
          dofU=nodes[n].dof[d]
          if dofU>0                 
            Ue[(k-1)*2+d]=U[dofU]    
          end 
        end
      end
  #    T6_Ke!(Ke,Xe,material[mat])
      T6_KeNL!(Ke,Fe,Xe,Ue,material[mat])
      for i = 1:12
        dofi=dofe[i]
        if dofi>0
          F[dofi]+=Fe[i]
          for j = 1:12
  #          F[dofi]-=Ke[i,j]*Ue[j]
            dofj=dofe[j]
            if dofj>0 
              K[dofi,dofj]+=Ke[i,j]
            end
          end
        end
      end
    end
  end 
  
  ###################################
  # follower forces
  ###################################
  
  # allocation of elements arrays
  Xe=zeros(Float64,3,2)
  dofUe=zeros(Int64,6,1)
  
  Fe=zeros(Float64,6,1)
  KUUe=zeros(Float64,6,6)
  
  for e in 1:info.NL
    Ue=zeros(Float64,3,2)                
    for k in 1:3                                 
      n=B3[e].nodes[k]
      Xe[k,:]=nodes[n].coor              
      dofUe[(k-1)*2+1:k*2]=nodes[n].dof
      for d in 1:2       
        dofU=nodes[n].dof[d]
        if dofU>0                      
          Ue[k,d]=U[dofU]    
        end 
      end
    end
  
    press=B3[e].press_damping
    B3_fl!(Fe,KUUe,Xe,Ue,press) 
  
    for i in 1:6
      dofUi=dofUe[i]
      if dofUi>0
        F[dofUi]+=Fe[i]
        for j in 1:6
          dofUj=dofUe[j]
          if dofUj>0
            K[dofUi,dofUj]+=KUUe[i,j]  
          end
        end
      end
    end    
  
  end
  
  
  # Solution phase
  println("Solving system")
  K=SparseMatrixCSC(K)  # convert to other sparse format
  sol = K\F
  # println("F")
  # println(Fe)
  # println(minimum(F))
  # println("sol")
  # println(maximum(sol))
  # println(minimum(sol))
  #ps = MKLPardisoSolver()
  #solve!(ps,sol,K,F)
  
  # fill dofs                
  U += sol     
   
  residuum=norm(sol)
  println("Residuum: ",residuum)
  
  end # while residuum
  
  return K
  
  end   