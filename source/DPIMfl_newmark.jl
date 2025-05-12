
function newmark(input_file)

fname = "../input/"*input_file
include(fname)
info.mesh_file=file

T6,B3,nodes = readgmsh!()

##############################
# static solution and matrix assemblage
##############################

K,M,Kfl=analysis(nodes,T6,B3)
if info.beta != 0
  K0 = damping_stiffness(nodes,T6,B3) # For stiffness damping
else
  K0 = K
end
#D,VR,VL=eigAB(K,M,Kfl,nodes)

neq=info.neq
u = zeros(Float64,neq,1)
v = zeros(Float64,neq,1)
upred = zeros(Float64,neq,1)
vpred = zeros(Float64,neq,1)
a=zeros(Float64,neq,1)
da=zeros(Float64,neq,1)

# fill dofs
for n in 1:info.NN 
  for d in 1:2       
    dofU=nodes[n].dof[d]
    if dofU>0                      
      u[dofU]=nodes[n].u[d]
    end 
  end
end

#uneq=0
## equation numbering 
#for e in 1:info.NE           
#  mat=T6[e].mat
#  if mat>0 
#    for k in 1:6
#      n=T6[e].nodes[k]   
#      for d in 1:2
#        if nodes[n].dof[d]==0   
#          uneq=uneq+1             
#          nodes[n].dof[d]=uneq   
#        end  
#      end
#    end
#  end   
#end

#neq=uneq
#info.neq=neq
#M = mass_assembly(nodes,T6)

dofnode=nodes[2].dof[2]

# Internal forces
F = zeros(Float64,neq,1)

dpress = 0.5
npress=1
press=7.21373583+dpress
dpress=0.1
vmax=zeros(Float64,npress)
vmin=zeros(Float64,npress)
vpress=zeros(Float64,npress)

Q=5
omega=0.68
# alpha=omega/Q
# betamod = 0
alpha = 0
betamod = 1/(omega*Q)
period=2*pi/omega
moltQ=10
nstepd=150
te = moltQ*Q*period

for ipress=1:npress

t=0
#fill!(u,0.0)  
fill!(v,0.0)  
fill!(a,0.0)  

vpress[ipress]=press

beta=.25
gamma=0.5
 
dt = period/nstepd
nstep=floor(Int,te/dt)+1
hist=zeros(Float64,nstep,1)
for i=1:nstep

  t+=dt

#  if t<10*period
#   alphamod=0.00001*alpha
#  else
   alphamod=alpha
#  end

  upred[:]=u+dt*v+dt^2*(0.5-beta)*a
  vpred[:]=v+dt*(1-gamma)*a

  residuum=100
  iter=0

  while residuum>1e-8
    
    iter+=1
    u[:]=upred+dt^2*beta*a
    v[:]=vpred+dt*gamma*a
    K=SparseMatrixLNK(neq,neq)
    compute_residual!(F,K,u,nodes,T6,B3,press)

    if t<period
     F[dofnode]+=0.01
    end
   
    K=SparseMatrixCSC(K)
    F[:]-=M*(a+alphamod*v)+K0*betamod*v
    #  da[:]=(M*(1+alphamod*gamma*dt)+K*beta*dt^2)\F
    #ps = MKLPardisoSolver()
    #solve!(ps,da,M*(1+alphamod*gamma*dt)+K*beta*(dt^2),F)
    da = (M*(1+alphamod*gamma*dt)+K*beta*(dt^2)+K0*betamod*gamma*dt)\F
    a[:] += da
    residuum=norm(da)
    #    println(residuum)

  end  

  hist[i]=u[dofnode]
#  println(t," ",t-period*floor(Int,t/period)," ",hist[i])
  check=t-period*floor(Int,(t+.0001*dt)/period)
#  println(check," ",floor(Int,t/period)," ",period," ",iter)
  if abs(check)<0.1*dt
    umax=maximum(hist[i-nstepd+1:i])
    umin=minimum(hist[i-nstepd+1:i])    
    vmax[ipress]=umax
    vmin[ipress]=umin
    println("MAX,MIN,MED,N    ",umax," ",umin," ",(umax-umin)/2," ",floor(Int,(t+.0001*dt)/period))
    file = matopen("./output/hist.mat","w")
    write(file, "hist",hist[1:i])
    close(file)
    file = matopen("./output/ampl.mat","w")
    write(file, "vpress",vpress)
    write(file, "vmax",vmax)
    write(file, "vmin",vmin)
    write(file, "vampl",(vmax-vmin)/2)
    close(file)
  end 

end  # end iterations

press+=dpress
println(press)

end # end iomega

end



function compute_residual!(F::Matrix{Float64},K::SparseMatrixLNK{Float64},u::Matrix{Float64},
                           nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},press::Float64)

  fill!(F,0.0)
  fill!(K,0.0)
  
  # allocation of elements arrays
  Xe=zeros(Float64,6,2)
  dofe=zeros(Int64,12)
  Ue=zeros(Float64,12)
  Fe=zeros(Float64,12)
  Ke= zeros(Float64,12,12)
  
  for e in 1:info.NE     
    mat=T6[e].mat
    if mat>0
      fill!(Ue,0.0) 
      for k in 1:6
        n=T6[e].nodes[k]   
        Xe[k,:]=nodes[n].coor
        dofe[2*k-1:2*k]=nodes[n].dof
        for d=1:2
         dofu=dofe[2*(k-1)+d]
         if dofu>0
          Ue[2*(k-1)+d]=u[dofu]
         end
        end  
      end
      T6_KeNL!(Ke,Fe,Xe,Ue,material[mat])
      for i = 1:12
        dofi=dofe[i]
        if dofi>0
          F[dofi]+=Fe[i]
          for j = 1:12
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
  Ke=zeros(Float64,6,6)
  
  for e in 1:info.NL                     
    fill!(Ue,0.0) 
    for k in 1:3                                 
      n=B3[e].nodes[k]
      Xe[k,:]=nodes[n].coor              
      dofUe[2*k-1:2*k]=nodes[n].dof
      for d=1:2
        dofu=dofUe[2*(k-1)+d]
        if dofu>0
         Ue[k,d]=u[dofu]
        end
      end  
    end
  
    B3_fl!(Fe,Ke,Xe,Ue,press) 
 #    B3_Ftot!(Fe,Ke,Xe,Ue,omegat) 
  
    for i in 1:6
      dofUi=dofUe[i]
      if dofUi>0
        F[dofUi]+=Fe[i]
        for j in 1:6
          dofUj=dofUe[j]
          if dofUj>0
            K[dofUi,dofUj]+=Ke[i,j] 
          end  
        end  
      end
    end    
  
  end
end





function mass_assembly(nodes::Vector{Snode},T6::Vector{ST6})

neq=info.neq  
M=SparseMatrixLNK(neq,neq)

# allocation of elements arrays
Xe=zeros(Float64,6,2)
dofe=zeros(Int64,12,1)
Me= zeros(Float64,12,12)

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
return M
end

