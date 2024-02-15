using SegregatedVMSSolver.ExportUtility
include("BoundaryConditions.jl")

function run_windtunnel(params,distribute)
    @unpack rank_partition, D, N, t0, dt, tF, ν, θ, M = params

    parts  = distribute(LinearIndices((prod(rank_partition),)))

    mesh_file_path =  params[:mesh_file]
    
    model = GmshDiscreteModel(parts, mesh_file_path)
    println("model read completed")
    add_new_tag!(model, params)
    


    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_windtunnel(params) 
    merge!(params, Dict(:u0 => u0, :model => model, :parts=>parts))
    print_model(params)

    V, U, P, Q, Y, X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    

    degree = 4*params[:order]
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    
    new_dict = Dict(:U => U,
    :P => P,
    :X => X,
    :Y => Y,
    :Ω => Ω,
    :dΩ => dΩ,
    :degree => degree)
    merge!(params, new_dict)

    trials = [U, P]
    tests = [V, Q]
  
    merge!(params, Dict(:trials => trials, :tests => tests))
  
    solve_case(params)




end


