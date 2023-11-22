module Utilstests

using Test

@testset "Utils" begin
  include("ReadAirfoilResultsTests.jl")
  include("WallDistanceTests.jl")
end


end