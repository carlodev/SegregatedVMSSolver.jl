module CreateVtuTests

using Test
using DataFrames
using SegregatedVMSSolver.CreateVtu


fname = create_vtu_file("Results_vtu")
@test typeof(fname) == String

end #end module