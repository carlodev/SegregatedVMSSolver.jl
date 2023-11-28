module AddNewTagsTests

using Test
using SegregatedVMSSolver
using Gridap
using GridapDistributed
using PartitionedArrays



function test_addtag(rank_partition,distribute,D)


parts  = distribute(LinearIndices((prod(rank_partition),)))
domain = (D==2) ? (0,1,0,1) : (0,1,0,1,0,1)
mesh_partition =  (D==2) ?  (40,40) : (40,40,40)
model = CartesianDiscreteModel(parts,rank_partition,domain,mesh_partition)
tag_coordinate = (D==2) ? VectorValue(0.5,0.5) : VectorValue(0.5,0.5,0.5)

SegregatedVMSSolver.add_centre_tag!(model, tag_coordinate)
labels = get_face_labeling(model)

@test !isempty(findall(x->x=="centre", labels.labels.items[1].tag_to_name))


@testset "New Tag" for lab in labels.labels.items 
    @test !isempty(findall(x->x=="centre", lab.tag_to_name))
end

end


for D in [2]
    rank_partition = (D==2) ?  (2,2) : (2,2,2)
    with_debug() do distribute
        test_addtag(rank_partition,distribute,D)
    end
end


end


