module CreateVtu

using FileIO
export create_vtu_file

function create_vtu_file(results_folder::String)
  ff = cd(readdir, results_folder)


  ff_dir = ff[1:2:end]
  ff_file = ff[2:2:end]
  
  
  ff_split = collect(split.(ff_dir, "_"))
  time_step = String[]
  ff_split[1][2] #[1:end-4]
  for i = 1:1:length(ff_split)
      tmp = ff_split[i][2]#[1:end-4]
      push!(time_step, tmp)
  end
  
  idx_sort = sortperm(parse.(Float64,time_step))
  
  initial_string = "<?xml version=\"1.0\" encoding=\"utf-8\"?>
  <VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">
    <Collection>"
  end_string = "  </Collection>
  </VTKFile>"
  
  fname = ff_split[1][1] * ".pvd"
  io = open(fname , "w")
  println(io, initial_string)
  for i in idx_sort
      println(io, "<DataSet timestep=\"$(time_step[i])\" part=\"0\" file=\"Results_vtu\\$(ff_file[i])\"/>")
  end
  println(io, end_string)
  close(io)


end

end