module RestartTests
using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays

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
    :method=>:VMS,
    :Cᵢ => [4, 36],
    :benchmark=>false,
    :t_endramp=> 5.0,
    :mesh_file => "DU89_2D_A1_M.msh",
    :TI => 0.001,
    :ρ=>1.0,
    :Re=> 1_000,
    :c=> 1.0,
    :restart=> true,
    :restart_file=>"BL_DU89_2D_A1_M.csv",     
)

SegregatedVMSSolver.init_params(params)
idx_v = SegregatedVMSSolver.find_idx(VectorValue(0.0,0.0), params)

@test typeof(idx_v) == Int64
@test(typeof(SegregatedVMSSolver.uh_r(VectorValue(0.0,0.0), params, idx_v))<:VectorValue)


end #end module