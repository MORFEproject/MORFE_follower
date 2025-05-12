using DelimitedFiles

function eigenvalue_trajectories(input_file, output_file, mu_min, mu_max, n_mu)

    fname = "../input/"*input_file
    include(fname)
    info.mesh_file = file

    eienvals_real = zeros(n_mu+1,2*info.neig+1)
    eienvals_imag = zeros(n_mu+1,2*info.neig+1)
    for i = 0:n_mu
        global tbcval = [mu_min + i * (mu_max - mu_min)/n_mu]
        eienvals_real[i+1,1] = tbcval[1]
        eienvals_imag[i+1,1] = tbcval[1]

        T6, B3, nodes = readgmsh!()

        K, M, F0 = analysis(nodes,T6,B3)
        if info.beta != 0
            K0 = damping_stiffness(nodes,T6,B3) # For stiffness damping
        else
            K0 = nothing
        end

        Lambda, VR, VL, VRp = eigAB(K,M,F0,nodes,K0)
        eienvals_real[i+1,2:end] = real(diag(Lambda))
        eienvals_imag[i+1,2:end] = imag(diag(Lambda))
    end

    isdir("./" * output_folder) || mkdir("./" * output_folder)
    output_file = split(output_folder,"/")
    output_file = output_file[end]
    for i = 1:2*info.neig
        writedlm( output_folder * "/" * output_file * "_real_$(i).txt",  eienvals_real[:,[1,i+1]], ' ')
        writedlm( output_folder * "/" * output_file * "_imag_$(i).txt",  eienvals_imag[:,[1,i+1]], ' ')
    end
end