module AddNewTagsTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays

using SegregatedVMSSolver.CreateProblem


function test_label(a::DebugArray)
    return !isempty(findall(x -> x == "centre", a.items[1].tag_to_name))

end

function test_label(a::MPIArray)
    return !isempty(findall(x -> x == "centre", a.item.tag_to_name))
end

function test_addtag(rank_partition, distribute, D)


    parts = distribute(LinearIndices((prod(rank_partition),)))
    domain = (D == 2) ? (0, 1, 0, 1) : (0, 1, 0, 1, 0, 1)
    mesh_partition = (D == 2) ? (11, 11) : (11, 11, 11)
    model = CartesianDiscreteModel(parts, rank_partition, domain, mesh_partition)
    tag_coordinate = (D == 2) ? VectorValue(0.5, 0.5) : VectorValue(0.5, 0.5, 0.5)

    SegregatedVMSSolver.CreateProblem.add_centre_tag!(model, tag_coordinate)
    labels = get_face_labeling(model)

    @testset "New Tag" begin
        @test test_label(labels.labels)
    end


end

function main(distribute)
    for D in [2]
        rank_partition = (D == 2) ? (2, 2) : (2, 2, 1)
        test_addtag(rank_partition, distribute, D)
    end

end



end #end module

