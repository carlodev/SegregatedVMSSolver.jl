
include("AnalyticalSolution.jl")
include("SpaceConditions.jl")

function run_taylorgreen(params, distribute)
  @unpack rank_partition, D, N, t0, dt, tF, ν, θ, M = params

  parts = distribute(LinearIndices((prod(rank_partition),)))


  diameter = 0.5 #0.5 [m] vortex dimension

  Vs = 1 #1[m/s]swirling speed
  Ua = 0.3 #0.3 [m/s]convective velocity in x
  Va = 0.2 #0.2 [m/s]convective velocity in y
  params[:ν] = 0.001 #0.001 m2/s 

  #Domain and mesh definition
  domain = (-diameter, diameter, -diameter, diameter)
  partition = (params[:N], params[:N])
  model = CartesianDiscreteModel(parts, rank_partition, domain, partition; isperiodic=(true, true))

  # hf_gen!(params)
  velocity, pa, ωa = analytical_solution(diameter, Vs, Ua, Va, params[:ν])
  merge!(params, Dict(:model => model, :force_tags=>nothing, :parts=>parts))
  V, Q, U, P, Y, X, model = CreateTGSpaces(model, params, pa) #We update model with the new label of the center point
  print_model(params)

  println("spaces created")

  degree = 4 * params[:order]
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  
  new_dict = Dict(:U => U,
    :P => P,
    :X => X,
    :Y => Y,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree,
    :force_params => nothing,
    :p0 => pa, :u0 => velocity)

  merge!(params, new_dict)



  #Create trials tests vectors 

  trials = [U, P]
  tests = [V, Q]

  merge!(params, Dict(:trials => trials, :tests => tests))


  solve_case(params)




end #end function



