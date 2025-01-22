using Gridap
using Gridap.FESpaces
using GridapDistributed
using PartitionedArrays


with_debug
function main_ex1(rank_partition,distribute)
    parts  = distribute(LinearIndices((prod(rank_partition),)))

n = 5
domain = (0,1,0,1)
partition = (n,n)

model = CartesianDiscreteModel(parts,rank_partition,domain,partition)



D = 2
order = 2
reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
V = TestFESpace(model,reffeᵤ,conformity=:H1,dirichlet_tags=["boundary"])



function ubound(x,t)
    # println(x)
    VectorValue(x[1]*t,x[2]*t)
end
ubound(t::Real) = x->ubound(x,t)
U = TransientTrialFESpace(V,[ubound])
Utn = U(0.0)
Utn1 = U(0.1)
# Utn.dirichlet_values
# Utn1.dirichlet_values

uhfe = FESpaces.interpolate( VectorValue(0.0,0.0), Utn)

# uhfe = interpolate( VectorValue(0.0,0.0), Utn)

my_f = map(uhfe.fields) do f
  interpolate( VectorValue(0.0,0.0), get_fe_space(f))

end


println(typeof(uhfe)<:FEFunction)

println(typeof(uhfe.fields))

println(typeof(Utn))

# Utn = Utn1
# Utn1 = U(0.2)

end


rank_partition = (2,2)
with_debug() do distribute
  main_ex1(rank_partition,distribute)
end