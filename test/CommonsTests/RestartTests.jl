module RestartTests
using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays
using NearestNeighbors

function main(distribute)
    restart_file = joinpath(@__DIR__, "..", "..", "restarts", "BL_DU89_2D_A1_M.csv")

    

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
    :petsc_options => petsc_options_default(),
    :sprob=>StabilizedProblem(),

    :benchmark=>false,
    :t_endramp=> 5.0,
    :mesh_file => "DU89_2D_A1_M.msh",
    :TI => 0.001,
    :ρ=>1.0,
    :Re=> 1_000,
    :c=> 1.0,
    :restart=> true,
    :restart_file=> restart_file,     
)


SegregatedVMSSolver.init_params(params)
tree = create_search_tree(params)
uh_0 = restart_uh_field(params,tree)
ph_0 = restart_ph_field(params,tree)
eval_point = Point(0.5,0.5)

@test typeof(tree)<:BruteTree
@test typeof(uh_0(eval_point))<:VectorValue
@test typeof(ph_0(eval_point))<:Float64

end

end #end module