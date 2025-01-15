"""
    create_initial_conditions(simcase::SimulationCase)

It creates the initial conditions for velocity and pressure. If `restart` is `true` then the velocity and the pressure field are interpoled on the specified DataFrame.  
"""
function create_initial_conditions(simcase::VelocityBoundaryCase, params::Dict{Symbol,Any})
    @unpack Utn, Ptn = params ## initial fields

    @sunpack restart,D,t0, u0 = simcase

    Ut0 = Utn
    Pt0 = Ptn

    uh0v = VectorValue(u0[1:D]...)
    uh0 = interpolate(uh0v, Ut0)
    ph0 = interpolate(0.0, Pt0)

    if restart
        @sunpack restartfile =simcase
        restart_df = DataFrame(CSV.File(restartfile))
        tree = create_search_tree(restart_df)
        uh_0 = restart_uh_field(D,tree,restart_df)
        ph_0 = restart_ph_field(tree,restart_df)
        uh0 = interpolate(uh_0, Ut0)
        ph0 = interpolate(ph_0, Pt0)
    end
    
    
    print_initial_conditions((uh0,ph0,uh0,uh0,ph0), simcase,params)

    return uh0,ph0
end

function create_initial_∇uh(uh0,params)
    @unpack ∇U = params ## initial 

    ∇uh =  interpolate(∇(uh0),∇U)

    return ∇uh
end

#TaylorGreenCase Natural and Periodic
function create_initial_conditions(simcase::TaylorGreen,params::Dict{Symbol,Any})
    @unpack U,P = params
    @unpack a_solution = simcase.bc_type
    @sunpack t0 =  simcase

    uh0 = interpolate(a_solution[:velocity](t0), U(t0))
    ph0 = interpolate(a_solution[:pressure](t0), P(t0))
    print_initial_conditions((uh0,ph0), simcase,params)
    return uh0,ph0
end



function print_initial_conditions(fields, simcase::SimulationCase,params::Dict{Symbol,Any})
    if simcase.simulationp.exportp.printinitial
      case = typeof(simcase)
      @unpack Ω = params
      @sunpack order = simcase
      dir = "Initial_Conditions"
      mkpath(dir)
      save_path = joinpath("Initial_Conditions","InitialCondition_$(case)_.vtu")
      
      
      writesolution(simcase, Ω, order, save_path,0.0, fields)
    end
end
