
function outgmsh(NN::Int64,nodes::Vector{Snode},T6::Vector{ST6},B3::Vector{SB3},filename::String)

  ofile = open("outgmsh/out"*filename*".msh","w")
  for i in 1:2
    mystring = "\$MeshFormat\n"*
                "2.2 0 8\n"*
                "\$EndMeshFormat\n"*
                "\$NodeData\n"*
                "1\n"*
                "\"U"*string(i)*"\"\n"*
                "1\n"*
                "0.0\n"*
                "3\n"*
                "0\n"*
                "1\n"
    write(ofile,mystring)
    write(ofile,string(NN)*"\n")
    for n in 1:NN                                  
      mystring=string(n)*" "*string(nodes[n].u[i])*"\n"
      write(ofile,mystring)
    end
    write(ofile,"\$EndNodeData\n")
  end
  
  close(ofile)
  
  ofile = open("outgmsh/opti"*filename*".geo","w")
  mystring = "General.BackgroundGradient=0;\n"*
             "Mesh.SurfaceEdges=0;\n"*
             "Mesh.SurfaceFaces=0;\n"*
             "View[0].Visible = 1;\n"*
             "View[0].VectorType = 5;\n"*
             "View[0].NbIso=21;\n"*
             "View[0].IntervalsType=3;\n"*
             "View[1].Visible=0;\n"*
             "View[1].VectorType = 5;\n"*
             "View[1].NbIso=21;\n"*
             "View[1].IntervalsType=3;\n"
  write(ofile,mystring)
  
  mystring= "Plugin(Scal2Vec).ViewX=0;\n"*
            "Plugin(Scal2Vec).ViewY=1;\n"*
            "Plugin(Scal2Vec).Run;\n"*
            "View[2].Name = \"|U|\";\n"*
            "View[2].Visible = 0;\n"*
            "View[2].IntervalsType=3;\n"*
            "View[2].VectorType = 5;\n"*
            "View[2].NbIso=21;\n"
  write(ofile,mystring)
  close(ofile)
  
  ofile = open("outgmsh/post"*filename*".msh","w")
  write(ofile,"Merge \"../input/"*info.mesh_file*"\";\n")
  write(ofile,"Merge \"./out"*filename*".msh\";\n")
  write(ofile,"Merge \"./opti"*filename*".geo\";\n")
  close(ofile)
  
return nothing
end




function outgmsheig(NN::Int64,nodes::Vector{Snode},V::Vector{Float64},filename::String)

ofile = open("outgmsh/out"*filename*".msh","w")
for i in 1:2
  mystring = "\$MeshFormat\n"*
              "2.2 0 8\n"*
              "\$EndMeshFormat\n"*
              "\$NodeData\n"*
              "1\n"*
              "\"U"*string(i)*"\"\n"*
              "1\n"*
              "0.0\n"*
              "3\n"*
              "0\n"*
              "1\n"
  write(ofile,mystring)
  write(ofile,string(NN)*"\n")
  for n in 1:NN                   
    val=0.0
    if nodes[n].dof[i]>0
      val=V[nodes[n].dof[i]]
    end  
    mystring=string(n)*" "*string(val)*"\n"
    write(ofile,mystring)
  end
  write(ofile,"\$EndNodeData\n")
end

close(ofile)

ofile = open("outgmsh/opti"*filename*".geo","w")
mystring = "General.BackgroundGradient=0;\n"*
           "Mesh.SurfaceEdges=0;\n"*
           "Mesh.SurfaceFaces=0;\n"*
           "View[0].Visible = 1;\n"*
           "View[0].VectorType = 5;\n"*
           "View[0].NbIso=21;\n"*
           "View[0].IntervalsType=3;\n"*
           "View[1].Visible=0;\n"*
           "View[1].VectorType = 5;\n"*
           "View[1].NbIso=21;\n"*
           "View[1].IntervalsType=3;\n"
write(ofile,mystring)

mystring= "Plugin(Scal2Vec).ViewX=0;\n"*
          "Plugin(Scal2Vec).ViewY=1;\n"*
          "Plugin(Scal2Vec).Run;\n"*
          "View[2].Name = \"|U|\";\n"*
          "View[2].Visible = 0;\n"*
          "View[2].IntervalsType=3;\n"*
          "View[2].VectorType = 5;\n"*
          "View[2].NbIso=21;\n"
write(ofile,mystring)          


close(ofile)

ofile = open("outgmsh/post"*filename*".msh","w")
write(ofile,"Merge \"../input/"*info.mesh_file*"\";\n")
write(ofile,"Merge \"./out"*filename*".msh\";\n")
write(ofile,"Merge \"./opti"*filename*".geo\";\n")
close(ofile)

return nothing
end



