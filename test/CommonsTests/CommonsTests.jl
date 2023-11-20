module CommonsTests

using Test

@testset "Commons" begin
  include("InitializeParamsTests.jl")
  include("LinearUtilitiesTests.jl")
  include("MatrixCreationTests.jl")
  include("StabParamstests.jl")
end


end