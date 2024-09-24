function boundary_velocities(simcase::VelocityBoundaryCase)   
    t_endramp = simcase.simulationp.timep.t_endramp
    u_in = simcase.simulationp.physicalp.u_in
    u_in_mag = simcase.simulationp.physicalp.u_in_mag
    
    D =  simcase.meshp.D
    @unpack TurbulenceInlet,Eddies,  Vboxinfo, Re_stress, R = simcase.simulationp.turbulencep



    """ No ramp
    u_free(x,t) = (D == 2) ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    u_free(t::Real) = x -> u_free(x,t)
    
    """
  
    u_in = u_in[1:D]

    uint(t) = (t < t_endramp) ? t/t_endramp : 1.0

    uin(t) = uint(t) .* u_in
    
    function u_free(x, t)
    (D == 2) ? VectorValue(uin(t)...) :  VectorValue(uin(t)...)
    end

    u_free(t::Real) = x -> u_free(x,t)

    #No generation of Eddies during ramping
    uin0(t) = (t < t_endramp) ? 0.0 : 1.0
    
    u_SEM(x,t) = (TurbulenceInlet == 0.0 ) ? u_free(x,t) : u_free(x,t) .+  uint(t) .*compute_fluctuation(x,t, (TurbulenceInlet,Eddies, u_in_mag, Vboxinfo, Re_stress,D, R))
    u_SEM(t::Real) = x -> u_SEM(x,t)

   


    function u_wall(x,t::Real)
        return VectorValue(zeros(D)...)
    end 

    u_wall(t::Real) = x -> u_wall(x,t)

    return u_free, u_SEM, u_wall
end


function create_boundary_conditions(simcase::VelocityBoundaryCase) 
    u_free,u_SEM, u_wall = boundary_velocities(simcase)
    create_boundary_conditions(simcase,u_free,u_SEM, u_wall) 
end

function create_boundary_conditions(simcase::Airfoil, u_free,u_SEM, u_wall) 
    simulationp = simcase.simulationp

    u_diri_tags,u_diri_values= create_boundary_conditions(simulationp, u_free,u_SEM, u_wall) 
    p_diri_tags = ["outlet"]

    p_diri_values = [0.0]
    return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values

end 

function create_boundary_conditions(simulationp::SimulationParameters{Inlet}, u_free,u_SEM, u_wall) 
    u_diri_tags=["inlet", "limits", "airfoil"]
    u_diri_values = [u_SEM, u_free, u_wall]

    return u_diri_tags,u_diri_values
end 

function create_boundary_conditions(simulationp::SimulationParameters{Internal}, u_free,u_SEM, u_wall) 
    u_diri_tags=["inlet", "sem", "limits", "airfoil"]
    u_diri_values = [u_free, u_SEM, u_free, u_wall]
    return u_diri_tags,u_diri_values
end 


function create_boundary_conditions(simcase::WindTunnel, u_free,u_SEM, u_wall) 
    u_diri_tags=["inlet", "limits", "airfoil"]
    u_diri_values = [u_SEM, u_wall, u_wall]
    p_diri_tags=["outlet"]
    p_diri_values = [0.0]
    return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end 


function create_boundary_conditions(simcase::Cylinder, u_free,u_SEM, u_wall) 
    u_diri_tags=["inlet", "limits", "cylinder"]
    u_diri_values = [u_SEM, u_SEM, u_wall]
    p_diri_tags=["outlet"]
    p_diri_values = [0.0]
    return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end 


function create_boundary_conditions(simcase::LidDriven, u_free,u_SEM, u_wall) 
    p_diri_values = 0.0
    u_diri_tags=["diri0", "diri1"]
    u_diri_values = [u_wall, u_free]

    p_diri_tags= ["p"]

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values 
end 


function create_boundary_conditions(simcase::TaylorGreen)
    @unpack analyticalsol = simcase
    @unpack pressure = analyticalsol

  u_diri_tags=[]
  u_diri_values = []
  p_diri_tags= ["centre"]
  p_diri_values=pressure
  return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end
