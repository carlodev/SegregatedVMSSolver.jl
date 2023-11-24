module ExportUtilityTests

using SegregatedVMSSolver
using PartitionedArrays
using Gridap
using GridapDistributed

function export_utilities_test(rank_partition,distribute)
    parts  = distribute(LinearIndices((prod(rank_partition),)))
    domain = (0,1,0,1)
    mesh_partition = (4,4)
    model = CartesianDiscreteModel(parts,rank_partition,domain,mesh_partition)
    order = 2
    u((x,y)) = (x+y)^order
    f(x) = -Δ(u,x)
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe,dirichlet_tags="boundary")
    U = TrialFESpace(u,V)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*order)

    Γ = BoundaryTriangulation(model; tags="boundary") 
    n_Γ = get_normal_vector(Γ)
    params = Dict(:Γ=>Γ,:n_Γ=>n_Γ, :parts=>parts, :force_tags=>["boundary"])

    global_unique_idx = SegregatedVMSSolver.export_nodes_glob(parts, Γ)
    local_unique_idx = SegregatedVMSSolver.get_local_unique_idx(parts, Γ)
    SegregatedVMSSolver.export_n_Γ(params, local_unique_idx, global_unique_idx)
end

rank_partition = (2,2)
with_debug() do distribute
    export_utilities_test(rank_partition,distribute)
end






end
