


function run_case_test(case::String; meshfile=" ", t_endramp=1.0, Reynolds=1000)

   
   params = Dict(
         :N => 100,
         :D => 2, #Dimension
         :order => 1, 
         :t0 => 0.0,
         :dt => 0.25,
         :tF => 1.0,
         :case => case,
         :θ => 0.5,
         :u_in=> 1.0,
        
         :backend => with_debug,  #or with_mpi() with_debug()
         :rank_partition=>(2,2),
         :ν => 0.001,
         :petsc_options => petsc_options_default(),
         :method=>:VMS,
         :Cᵢ => [4, 36],
    
         :t_endramp=>t_endramp,
         :mesh_file => meshfile,
         :Re=> Reynolds,
         
         :c=> 1.0,
 
   )
   
   
   
   @test SegregatedVMSSolver.main(params) == true
   
   end