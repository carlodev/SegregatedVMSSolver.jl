module CommonsTests

using Test

@testset "Commons" begin
  include("InitializeParamsTests.jl")
  
#   include("CommonsProcedureTests.jl")
  include("AddNewTagsTests.jl")

  include("StabParamsTests.jl")
  include("LinearUtilitiesTests.jl")
  include("StabilizedEquationsTests.jl")
  

  include("MatrixCreationTests.jl")
  
  include("RestartTests.jl")
#   include("ExportUtilitiesTests.jl")


end


end