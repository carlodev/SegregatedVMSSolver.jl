using Test

@testset "Utils" begin
  include("ReadAirfoilResultsTests.jl")
  include("WallDistanceTests.jl")
  include("CreateVtuTests.jl")
end
