module InitialConditions
using Parameters
using Gridap
using GridapDistributed
using PartitionedArrays
using CSV
using DataFrames

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.Restart
using SegregatedVMSSolver.ExportUtility

export create_initial_conditions


"""
  create_initial_conditions(simcase::SimulationCase)

It creates the initial conditions for velocity and pressure. If `restart` is `true` then the velocity and the pressure field are interpoled on the specified DataFrame.  
"""
function create_initial_conditions(simcase::VelocityBoundaryCase,params::Dict{Symbol,Any})
    @unpack U,P = params

    restart,D,t0 =  get_field(simcase,[:restart,:D,:t0])

    uh0v = VectorValue(zeros(D)...)
    uh0 = interpolate_everywhere(uh0v, U(t0))
    ph0 = interpolate_everywhere(0.0, P(t0))

    if restart
        restartfile = get_field(simcase,:restartfile)
        restart_df = DataFrame(CSV.File(restartfile))
        tree = create_search_tree(params,restart_df)
        uh_0 = restart_uh_field(D,tree,restart_df)
        ph_0 = restart_ph_field(tree,restart_df)
        uh0 = interpolate_everywhere(uh_0, U(t0))
        ph0 = interpolate_everywhere(ph_0, P(t0))
    end

    print_initial_conditions((uh0,ph0,uh0,uh0,ph0), simcase,params)

    return uh0,ph0
end


#TaylorGreenCase
function create_initial_conditions(simcase::TaylorGreen,params::Dict{Symbol,Any})
    @unpack U,P = params
    @unpack analyticalsol = simcase

    uh0 = interpolate_everywhere(analyticalsol.velocity, U(t0))
    ph0 = interpolate_everywhere(analyticalsol.pressure, P(t0))
    print_initial_conditions((uh0,ph0), simcase,params)
    return uh0,ph0
end

function print_initial_conditions(fields, simcase::SimulationCase,params::Dict{Symbol,Any})
    if simcase.simulationp.exportp.printinitial
      save_sim_dir = simcase.simulationp.exportp.save_sim_dir
      case = typeof(simcase)
      @unpack Ω = params

      save_path = joinpath(save_sim_dir,"InitialCondition_$(case)_.vtu")
      
      writesolution(simcase, Ω, save_path,0.0, fields)
    end
end



end #end module