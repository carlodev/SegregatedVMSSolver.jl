"""
    verifykey(params::Dict{Symbol,Any},keyname; val = false)

It check if the dictionary params has the entry `keyname`. If not it adds the new entry with the value val. It is used to add default values
"""
function verifykey(params::Dict{Symbol,Any},keyname; val = false)
    if !haskey(params, keyname)
        merge!(params,Dict(keyname=>val))
    end
end

"""
    init_params(params::Dict{Symbol,Any})

It initialize the parametes of the simulations. Default values are added if not specified by the user.
In case of a restarting simulation it reads the file and create a `initial_rescale_factor`<1 for CFL ramping initialization.

"""
function init_params(params::Dict{Symbol,Any})

    @unpack D,rank_partition,case,t0,dt,tF = params

    @assert length(rank_partition) == D
    
    verifykey(params,:ρ; val = 1.0)

    if case == "Airfoil" || case == "LidDriven"  || case == "Cylinder"
        @unpack c,ρ,u_in,Re,ν = params

        params[:ν] = c*ρ*u_in/Re
        println("Recompute ν = $(params[:ν])")
    end

    #Time Step vector creation
    time_step = t0+dt:dt:tF 
    merge!(params,Dict(:time_step=>time_step))


    
  

    verifykey(params,:printmodel)
    verifykey(params,:printinitial)
    verifykey(params,:matrix_freq_update; val = 20) #Matrix update every 20 time steps
    verifykey(params,:a_err_threshold; val = 200) #Norm on accelaration vector reduced by a factor of 200
    verifykey(params,:benchmark; val = true) # do not print the results
    verifykey(params,:M; val = 20)
    verifykey(params,:t_endramp; val = 0.0)
    verifykey(params,:TI; val = 0.0)
    verifykey(params,:mesh_file; val = " ")
    verifykey(params,:restart)
    verifykey(params,:restart_file; val = " ")


    if params[:restart]
        @unpack restart_file, t_endramp, t0 = params
        restart_path = joinpath(@__DIR__, "../../restarts", restart_file)
        restart_df = DataFrame(CSV.File(restart_path))
        
        initial_rescale_factor = 1.0
        
        if t_endramp>t0
            initial_rescale_factor = t0/t_endramp
        end
        params = merge!(params, Dict(:restart_df => restart_df, :initial_rescale_factor=>initial_rescale_factor))
    end

end