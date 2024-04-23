module SegregatedVMSSolver


include(joinpath("Commons","SolverOptions.jl"))
include(joinpath("Commons","Interfaces.jl"))

include(joinpath("Commons","ParametersDef","ParametersDef.jl"))
include(joinpath("Commons","ExportUtility.jl"))

include(joinpath("Commons","Equations","Equations.jl"))
include(joinpath("Commons","CreateProblem","CreateProblem.jl"))

include(joinpath("Commons","VectorsOperations.jl"))

include(joinpath("Commons","MatrixCreation.jl"))


include(joinpath("Commons","SolveProblem.jl"))

include("Main.jl")


include(joinpath("utils","ReadAirfoilResults.jl"))
include(joinpath("utils","WallDistance.jl"))


include("Exports.jl")

end
