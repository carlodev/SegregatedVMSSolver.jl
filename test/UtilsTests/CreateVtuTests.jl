module CreateVtuTests

using Test
using DataFrames
using SegregatedVMSSolver.CreateVtu


fname = create_vtu_file(joinpath(@__DIR__,"Results_vtu"))
@test typeof(fname) == String

end #end module