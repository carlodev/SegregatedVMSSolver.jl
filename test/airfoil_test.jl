
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
         :M=> 8, #internal iterations
         :backend => backend,  #or with_mpi() with_debug()
         :rank_partition=>(2,1),
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
         :printmodel=>true,
         :name_tags=>["airfoil"],
         :fieldexport=>[("ph","friction")],
         :time_window=>(0.25,1.0),

   )
   
   

   @test SegregatedVMSSolver.main(params) == true
   
   end

run_airfoil_test(with_mpi)

# mpiexecjl --project=. -n 4 julia test/airfoil_test.jl
