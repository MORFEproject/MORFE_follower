function output(P::Vector{Parametrisation}, folder::String)


# writes complex parametrization on file
ofile = open("./" * folder * "/outCFULL.txt","w")

uneq=info.neq

write(ofile,"\nParametrization f\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,string(P[p].Av[i])*"\n")
    for j in 1:info.nza
      write(ofile,string(P[p].f[j,i])*"\n")
    end  
  end
end

write(ofile,"\nParametrization W\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,"\n"*string(P[p].Av[i])*"\n")
    for j in uneq+1:2*uneq+1
      write(ofile,string(P[p].W[j,i])*"\n")
    end  
  end
end
close(ofile)               

# writes real parametrization on file
ofile = open("./" * folder * "/outRFULL.txt","w")
write(ofile,"\nParametrization fr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,string(P[p].Av[i])*"\n")
    for j in 1:info.nmm
      write(ofile,string(2*real(P[p].fr[j,i]))*"\n")
    end  
    for j in 1:info.nmm
      write(ofile,string(2*imag(P[p].fr[j,i]))*"\n")
    end  
    write(ofile,string(P[p].fr[2*info.nmm+1,i])*"\n")
  end
end

write(ofile,"\nParametrization Wr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:P[p].m # for every alpha vector
    write(ofile,"\n"*string(P[p].Av[i])*"\n")
    for j in uneq+1:2*uneq+1
      write(ofile,string(real(P[p].Wr[j,i]))*"\n")
    end  
  end
end
close(ofile)               

write_rdyn(info,P,folder)

end



function write_rdyn(info::Sinfo,P::Vector{Parametrisation},folder::String)

  # additional paramter has no dynamics
  rdyn = ["" for i in 1:info.nza]
  for i = 1:info.nza
    rdyn[i] = "a"*string(i)*"' = "
  end

  for p in 1:info.max_order
    for c = 1:P[p].m
      Av=P[p].Av[c]   
      monomial = ""
      for d = 1:info.nza+1        
        if (Av[d]!=0)
          monomial *= "*a"*string(d)*"^"*string(Av[d])
        end
      end
      for j in 1:info.nmm
        rcoeff=2*real(P[p].fr[j,c])
#        
#        icoeff=2*imag(P[p].fr[j,c]) # ORIGINAL
        icoeff=-2*imag(P[p].fr[j,c])
#
        if abs(rcoeff)>1e-20
          rdyn[j] *= " + "*string(rcoeff)*monomial
        end
        if abs(icoeff)>1e-20
          rdyn[j+info.nmm] *= " + "*string(icoeff)*monomial
        end
      end
#      # odd variable
#      rdyn[2*info.nmm+1] *= " + "*string(real(P[p].fr[2*info.nmm+1,c]))*monomial
    end
  end

  ofile = open("./" * folder * "/equations.txt","w")
  for i = 1:info.nza
    write(ofile,rdyn[i]*";\n")
  end
  close(ofile)  

end

