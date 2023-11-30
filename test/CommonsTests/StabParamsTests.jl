module StabParamsTests
using Test
using SegregatedVMSSolver
using Gridap
using LinearAlgebra
using PartitionedArrays


function test_field(fd, val)
    map(fd.fields) do f
        @test norm(f.cell_field.value.value - val) < 1e-10
    end
    return true
end

function test_stab_params(distribute)

    N = 100
    domains = [(0, 1, 0, 1), (0, 1, 0, 1, 0, 1)]
    rp = [(2, 2), (2, 2, 1)]

    for (i, D) in enumerate([2, 3])

        rank_partition = rp[i]

        parts = distribute(LinearIndices((prod(rank_partition),)))

        domain = domains[i]
        partition = N .* ones(D)
        model = CartesianDiscreteModel(parts, rank_partition, domain, partition)

        Ω = Triangulation(model)
        he = SegregatedVMSSolver.h_param(Ω, D)

        @testset "he param" begin @test test_field(he, 1 / N) end

        test_field(he, 1 / N)


        G, GG, gg = SegregatedVMSSolver.G_params(Ω, Dict(:D => D))

        @testset "G param" begin
            if D == 2
                @test test_field(G, TensorValue(N^2, 0, 0, N^2))
            elseif D == 3
                @test test_field(G, TensorValue(N^2, 0, 0, 0, N^2, 0, 0, 0, N^2))
            end

        end
        @testset "GG param" begin @test test_field(GG, D * N^4) end
        @testset "gg param" begin  @test test_field(gg, D .* N^2) end



    end #end for 
end #end function 

function main(distribute)
    test_stab_params(distribute)
end

end #end module 