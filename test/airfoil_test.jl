
using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays

function run_airfoil_test(backend)

   mesh_file = joinpath(@__DIR__, "..", "models", "DU89_2D_A1_M.msh")
   restart_file = joinpath(@__DIR__, "..", "restarts", "BL_DU89_2D_A1_M.csv")
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
         :backend => backend,  #or with_mpi() with_debug()
         :rank_partition=>(2,2),
         :ν => 0.001,
         :petsc_options => petsc_options_default(),
         :method=>:VMS,
         :Cᵢ => [4, 36],
         :benchmark=>false,
         :t_endramp=> 5.0,
         :mesh_file =>  mesh_file,
         :TI => 0.001,
         :ρ=>1.0,
         :Re=> 1_000,
         :c=> 1.0,
         :restart=> false,
         :restart_file=> restart_file,
         :printmodel=>false,
         :name_tags=>["airfoil","inlet"]
   )
   
   
   
   @test SegregatedVMSSolver.main(params) == true
   
   end
   
   #mpiexecjl --project=. -n 4 julia airfoil_test.jl
   