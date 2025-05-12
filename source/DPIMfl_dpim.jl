function dpim(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},
              Lambda::Matrix{ComplexF64},VR::Matrix{ComplexF64},VL::Matrix{ComplexF64},VRp::Vector{Float64},
              nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},K0=nothing)

  info.nmm=length(info.Lmm)   # master modes
  nmm=info.nmm

  uneq=info.neq
  nA=2*uneq
  info.nA=nA
  info.nza=2*nmm         
 
  info.nMat=nA+info.nza  # dim of system to be solved

  println("Init Parametrisation")
  P=initParametrisation!(info)

  println("Init System")
  Rhs = Array{ComplexF64}(undef,info.nMat)
  Sol = Array{ComplexF64}(undef,info.nMat)
  Mat=spzeros(ComplexF64,info.nMat,info.nMat)
  
  println("Order 1")
  # Changing eigenvectors ordering
  perm = [if (i%2!=0) i÷2+1 else i÷2+nmm end for i in 1:2*nmm]
  Lambda = Lambda[perm,perm]
  VR = VR[:,perm]
  VL = VL[:,perm]

  # Order 1 parametrisation
  P[1].f[1:2*nmm,1:2*nmm] = Lambda
  P[1].W[1:2*uneq,1:2*nmm] = VR
  B = spzeros(Float64,nA,nA)
  B[1:uneq,1:uneq] = M
  B[uneq+1:2*uneq,uneq+1:2*uneq] = M
  BY = B*VR
  XTB = transpose(VL)*B

# places autovett=0 in last position
  P[1].W[uneq+1:2*uneq,info.nza+1] = VRp[:]
  P[1].W[2*uneq+1,info.nza+1] = 1

  Lambda=Vector{ComplexF64}(undef,info.nza+1)
  for i in 1:info.nza
    Lambda[i]=P[1].f[i,i] 
  end
  Lambda[info.nza+1]=0.0

# no nonautonomous
  
  println("Higher orders")
  for p = 2:info.max_order
    println("Order $p")
    fillrhs_quad!(nodes,T6,B3,P,p)
    fillrhs_cub!(nodes,T6,B3,P,p)
    fillWf!(P,p)
    for i in 1:P[p].m # for every alpha vector
      corresp=P[p].corresp[i]
      if corresp>0
        fillWfnonaut!(P,p,P[p].Av[i],i)
        homological!(Sol,Rhs,Mat,Lambda,P[p],i,K,M,BY,XTB,K0)
        P[p].analysed[i]=1
      elseif corresp<0
        P[p].W[:,i]=conj(P[p].W[:,-corresp])
        P[p].f[1:nmm,i]=conj(P[p].f[nmm+1:2*nmm,-corresp])
        P[p].f[nmm+1:2*nmm,i]=conj(P[p].f[1:nmm,-corresp])
        P[p].analysed[i]=1
      end  
    end  

  end

 return  P

end



function homological!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                      Lambda::Vector{ComplexF64},P::Parametrisation,Apos::Int64,
                      K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},
                      BY::Matrix{ComplexF64},XTB::Matrix{ComplexF64},K0)

             
  uneq=info.neq
  nA=info.nA
  nza=info.nza

  σ=dot(P.Av[Apos],Lambda[1:nza+1])
#  resonant_modes = zeros(Bool,info.nza)
#  fill!(resonant_modes,false)
#  for i = 1:info.nza
#    λ = Lambda[i]
##    if (abs(σ-λ)/abs(λ)<=info.tol)  # does not work is one eigen is 0
#    if (abs(σ-λ)<=info.tol)
#      resonant_modes[i] = true
#    end
#  end

LambdaSym=info.LambdaSym
  σ_sym = dot(P.Av[Apos],LambdaSym)
  resonant_modes = zeros(Bool,info.nza+1)
  if info.style == 'c'
    fill!(resonant_modes,false)
    for i = 1:info.nza+1
      λ_sym = LambdaSym[i]
      if ( σ_sym == λ_sym )
        resonant_modes[i] = true
      end
    end
  else
    fill!(resonant_modes,true)
  end
  println(P.Av[Apos],"  ",resonant_modes)

  fill!(Rhs,0.0)
  Rhs[1:nA]=P.R[:,Apos]
  
  Rhs[1:uneq]-=M*P.Wf[1:uneq,Apos]
  Rhs[uneq+1:2*uneq]-=M*P.Wf[uneq+1:2*uneq,Apos]

  fill!(Mat.nzval,0.0)
  Mat[1:uneq,1:uneq] = σ*M+info.alpha*M
  Mat[uneq+1:2*uneq,uneq+1:2*uneq] = σ*M
  if K0 !== nothing
    Mat[1:uneq,1:uneq] += info.beta*K0 
  end

  Mat[uneq+1:2*uneq,1:uneq] = -M
  Mat[1:uneq,uneq+1:2*uneq] = K
  for j = 1:info.nza
    if resonant_modes[j]  # if resonant
      Mat[1:nA,nA+j]=BY[:,j]
      Mat[nA+j,1:nA]=XTB[j,:]
    else
      Mat[nA+j,nA+j]=1   
    end  
  end  

  #  println("Solving system2")
  #ps = MKLPardisoSolver()
  #solve!(ps,Sol,Mat,Rhs)
  Sol = Mat\Rhs

  P.W[1:nA,Apos]=Sol[1:nA]
  P.f[1:info.nza,Apos]=Sol[nA+1:nA+info.nza]

end


function fillrhs_quad!(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},
                       P::Vector{Parametrisation},p::Int64)

  nA=info.nA
  uneq=info.neq
  for p1 in 1:p-1 
    p2=p-p1
    for k1 in 1:P[p1].m, k2 in 1:P[p2].m 
      Av=P[p1].Av[k1]+P[p2].Av[k2]   
      pos=findfirst(x->x==Av,P[p].Av)
      if P[p].corresp[pos]>0
        W1 = P[p1].W[uneq+1:nA,k1] # second part is U+psi, first is V
        W2 = P[p2].W[uneq+1:nA,k2]
        press1 = P[p1].W[nA+1,k1]
        press2 = P[p2].W[nA+1,k2]
        assembly_quad!(P[p],pos,W1,W2,nodes,T6)
        assembly_cub0!(P[p],pos,W1,W2,nodes,T6)
        assembly_quad_fl!(P[p],pos,W1,W2,press1,press2,nodes,B3)
      end 
    end
  end  

end


function fillrhs_cub!(nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},
                      P::Vector{Parametrisation},p::Int64)
  
  nA=info.nA
  uneq=info.neq
  for p1 in 1:p-2, p2 in 1:p-1-p1  #, p3 in 1:p-2
    p3=p-p1-p2
    for k1 in 1:P[p1].m, k2 in 1:P[p2].m, k3 in 1:P[p3].m  
      Av=P[p1].Av[k1]+P[p2].Av[k2]+P[p3].Av[k3]   
      pos=findfirst(x->x==Av,P[p].Av)
      if P[p].corresp[pos]>0
        W1 = P[p1].W[uneq+1:nA,k1]
        W2 = P[p2].W[uneq+1:nA,k2]
        W3 = P[p3].W[uneq+1:nA,k3]
        assembly_cub!(P[p],pos,W1,W2,W3,nodes,T6)
      end
    end
  end 

end


function fillWf!(P::Vector{Parametrisation},p::Int64)

  for p1 in 2:p-1, p2 in 2:p-1
    if (p1+p2)==p+1  
      for i in 1:P[p1].m
        Av1=P[p1].Av[i][:]
        for j in 1:P[p2].m 
          Av2=P[p2].Av[j][:]      
          for s in 1:info.nza+1
            if Av1[s]>0
              Av=Av1+Av2
              Av[s]-=1
              pos=findfirst(x->x==Av,P[p].Av)
#              P[p].Wf[1:neq,pos]+=Av1[s]*P[p1].W[1:neq,i]*P[p2].f[s,j]
#              P[p].Wf[neq+1:2*neq,pos]+=Av1[s]*P[p1].W[neq+1:2*neq,i]*P[p2].f[s,j]
               P[p].Wf[:,pos]+=Av1[s]*P[p1].W[:,i]*P[p2].f[s,j]
            end  
          end    
        end 
      end
    end
  end

end

function fillWfnonaut!(P::Vector{Parametrisation},p::Int64,Av::Vector{Int64},ind_rhs::Int64)
  # if ind_rhs == 71 || ind_rhs == 36 || ind_rhs == 16 || ind_rhs == 6 || ind_rhs == 2 || ind_rhs == 1
  #   println("Entered ", ind_rhs)
  # end
  for s in 1:info.nza-1
    for j in s+1:info.nza
      if Av[j]>0
        fs_j = P[1].f[s,j]
        # if abs(fs_j)>10^(-8)    # only fills if the reduced dyn is nzero
          Av_W = Av[:]
          Av_W[s] += 1
          Av_W[j] -= 1
          ind_W = findfirst(x->x==Av_W,P[p].Av)
          P[p].Wf[:,ind_rhs] += P[p].W[:,ind_W]*fs_j*Av_W[s]
          # if ind_rhs == 71 || ind_rhs == 36 || ind_rhs == 16 || ind_rhs == 6 || ind_rhs == 2 || ind_rhs == 1
          #   println("Av = ", Av) 
          #   println("Av_W = ", Av_W) 
          #   println("ind_W = ", ind_W)
          #   println("W = ", P[p].W[180,ind_W])
          #   println("fs_j = ", fs_j)
          #   println("Av_W[s] = ", Av_W[s])
          #   println("Wf[:,ind_rhs] = ", P[p].W[180,ind_W]*fs_j*Av_W[s])
          # end
        # end
      end
    end 
  end
end
    



