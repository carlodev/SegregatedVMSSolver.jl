


function run_airfoil_test()
    petsc_options = " -vel_ksp_type gmres -vel_pc_type gamg -vel_ksp_rtol 1.e-10 -vel_ksp_converged_reason \
                      -pres_ksp_type cg -pres_pc_type gamg  -pres_ksp_rtol 1.e-6 -pres_ksp_converged_reason \
                      -ksp_atol 0.0"
   
   params = Dict(
         :N => 100,
         :D => 2, #Dimension
         :order => 1, 
         :t0 => 0.0,
         :dt => 0.25,
         :tF => 1.0,
         :case => "Airfoil",
         :θ => 0.5,
         :u_in=> 1.0,
         :M=> 20, #internal iterations
         :backend => with_debug,  #or with_mpi() with_debug()
         :rank_partition=>(2,2),
         :ν => 0.001,
         :petsc_options => petsc_options,
         :method=>:VMS,
         :Cᵢ => [4, 36],
         :benchmark=>false,
         :t_endramp=> 5.0,
         :mesh_file => "DU89_2D_A1_M.msh",
         :TI => 0.001,
         :ρ=>1.0,
         :Re=> 1_000,
         :c=> 1.0,
         :restart=> false,
         :restart_file=>"BL_DU89_2D_A1_M.csv",     
   )
   
   
   
   @test SegregatedVMSSolver.main(params) == true
   
   end
   
   #mpiexecjl --project=. -n 4 julia run_SegregatedSolver.jl
   