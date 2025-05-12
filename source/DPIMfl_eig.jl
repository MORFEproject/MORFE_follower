

function eigAB(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},F0::Vector{Float64},nodes::Vector{Snode},K0)

println("Computing eigenvalues")

uneq=info.neq
nA=2*uneq
neig=2*info.neig # adds zero one

B=spzeros(Float64,nA,nA)
B[1:uneq,1:uneq] = M
B[uneq+1:2*uneq,uneq+1:2*uneq] = M

A=spzeros(Float64,nA,nA)
A[1:uneq,1:uneq] = -info.alpha*M 
A[uneq+1:2*uneq,1:uneq] = M
A[1:uneq,uneq+1:2*uneq] = -K
if K0 !== nothing
  A[1:uneq,1:uneq] -= info.beta*K0 
end

mat"""
[$VR,$DR] = eigs($A,$B,$neig,'SM','Tolerance',10^-17);
$DR=diag($DR);
%$VR=$VR*10;
$VR=$VR*5;
[$VL,$DL] = eigs(transpose($A),transpose($B),$neig,'SM','Tolerance',10^-17);
$DL=diag($DL);
"""

println(DR)
println(DL)

# normalisation
for i = 1:neig
#  for j = 1:neig   
    cc = transpose(VL[:,i])*B*VR[:,i]  # compl conj
#    cd = transpose(VL[:,i])*A*VR[:,j]
#    println(i,",",j)
#    println(cc)
#    println(cd/DL[i])
#  end   
  for j = 1:nA   # arbitrary: scales only VR  (seems to coincide with old approach)
    VR[j,i] /= sqrt(cc)
    VL[j,i] /= sqrt(cc)
  end
  outgmsheig(info.NN,nodes,real(VR[uneq+1:2*uneq,i]),"_eig"*string(i))  # prints displacements
end

# define new eigenvectors to form Jordan blocks

# identify where similar eigenvalues are
# supposing the maximum multiplicity is two
if info.Jordan
  if info.symbolic_Jordan
    Hopf_pairs = info.Hopf_pairs
  else
    Hopf_pairs = Matrix{Int64}(undef, 0, 2)
    for i in eachindex(DR)
      for j = i+1:length(DR)
        if abs(DR[i] - DR[j]) < info.tolHopf
          Hopf_pairs = [Hopf_pairs; [i j]]
          break
        end
      end
    end
    info.Hopf_pairs = Hopf_pairs
  end

  dD = DR[Hopf_pairs[:,1]] - DR[Hopf_pairs[:,2]]
  gamma = 0*dD
  for i in eachindex(dD)
    gamma[i] = VR[:,Hopf_pairs[i,1]]'*VR[:,Hopf_pairs[i,2]]/(VR[:,Hopf_pairs[i,1]]'*VR[:,Hopf_pairs[i,1]])
  end

  mu = (1.0+0.0*im)*sparse(I,neig,neig)
  for i = 1:length(dD)
    mu[Hopf_pairs[i,1],Hopf_pairs[i,2]] = 1/dD[i];
    mu[Hopf_pairs[i,2],Hopf_pairs[i,2]] = -1/gamma[i]/dD[i];
  end

  VL = VL*transpose(inv(transpose(VL)*B*VR))
  nu = transpose(inv(transpose(VL)*B*VR*mu));
  Y = VR*mu;
  X = VL*nu;
end

# check they respect wanted properties:
# if Jordan
#   println("X.'*B*Y : ")
#   println(  transpose(X)*B*Y  ) # should be identity

#   println("X.'*A*Y : ")
#   println(  transpose(X)*A*Y  ) # should have ones out of diagonal
# end
# Lambda = transpose(X)*A*Y

VRp=zeros(Float64,uneq)
#ps = MKLPardisoSolver()
#solve!(ps,VRp,K,F0)
VRp = K\F0
outgmsheig(info.NN,nodes,VRp,"_eig"*string(neig+1))  # prints displacements

if info.Jordan
  return transpose(X)*A*Y,Y,X,VRp
else
  return diagm(DR),VR,VL,VRp
end

end


function eig(K::SparseMatrixCSC{Float64},M::SparseMatrixCSC{Float64},nodes::Vector{Snode})

  println("Computing eigenvalues")
  
  #mat"""
  #[$V,$D] = eigs($K,$M,5,'SM');
  #$D=diag($D);
  #"""
 
  D,V = eigs(K,M,nev=info.neig,which=:SM)

  println(D.^0.5)
  D = real(D)
  V = real(V)
  
  # normalisation
  for i = 1:info.neig
    c = transpose(V[:,i])*M*V[:,i]
    for j = 1:info.neq
      V[j,i] /= sqrt(c)
    end
    outgmsheig(info.NN,nodes,V[:,i],"_eig"*string(i))
  end
  
  return D,V
  
  end