module InitializeParamsTests

using Test
using SegregatedVMSSolver
using Parameters

c = 1.0
ρ=1.0
u_in= 1.0
Re = 10_000
ν=-1.0
case = "Airfoil"
D = 2
rank_partition=(2,2)
t0 = 0.0
dt = 0.1
tF = 10.0

params= Dict(:case=>case, :D=>D, :rank_partition=>rank_partition, 
:c=>c, :ρ=>ρ, :u_in=>u_in, :Re=>Re, :ν=>ν, :t0=>t0, :dt=>dt, :tF=>tF)
SegregatedVMSSolver.init_params(params)


@unpack printmodel, printinitial, matrix_freq_update,a_err_threshold,benchmark,M,t_endramp,mesh_file,restart,restart_file,ν = params
@test printmodel == false
@test printinitial == false
@test matrix_freq_update == 20
@test a_err_threshold == 200
@test benchmark == true
@test M == 20
@test t_endramp == 0.0
@test mesh_file == " "
@test restart == false
@test restart_file == " "
@test ν == c*ρ*u_in/Re

end