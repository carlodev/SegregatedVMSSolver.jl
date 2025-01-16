using Plots
using CSV
using DataFrames

ν=1/1600
dt = 0.01
res =  CSV.read("TGV_output.csv", DataFrame; header=[:time,:Ek,:Enstropy ])    
vol = (2*pi)^3

bk =  CSV.read("C:/Users/Carlo Brunelli/Documents/GitHub/ExoFlow-Results/BenchmarkDir/TGV/RE1600/TGVDiss.csv", DataFrame; header=[:time,:dissipation ])    

dEkdt = - diff(res.Ek)./dt ./ vol
Enst = res.Enstropy[2:end] * ν ./vol

Plots.default(linewidth=2)
plot(res.time[2:end], dEkdt, label="dEk/dt")
plot!(res.time[2:end],Enst, label="Enstrophy")
plot!(bk.time, bk.dissipation, label="Benchmark", linecolor=:black)


1/1600