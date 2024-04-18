function boundary_velocities(simcase::VelocityBoundaryCase)
    D, u_in, t_endramp = get_field(simcase,[:D, :u_in, :t_endramp])
   

    """ No ramp
    u_free(x,t) = (D == 2) ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    u_free(t::Real) = x -> u_free(x,t)
    
    """
  

    uin(t) = (t < t_endramp) ? (u_in .*(t/t_endramp)) : u_in

    u_free(x, t) = (D == 2) ? VectorValue(uin(t), 0.0) :  VectorValue(uin(t), 0.0, 0.0)
    u_free(t::Real) = x -> u_free(x,t)

    #No generation of Eddies during ramping
    uin0(t) = (t < t_endramp) ? 0.0 : 1
    

    u_SEM(x,t) = u_free(x,t) #: u_free(x,t) #.+ uin0(t) .* generation_u_fluct!(x,t, params[:sem_cache])
    u_SEM(t::Real) = x -> u_SEM(x,t)
    


    u_wall(x,t) = (D == 2) ? VectorValue(0.0, 0.0) : VectorValue(0.0, 0.0, 0.0) 
    u_wall(t::Real) = x -> u_wall(x,t)

    return u_free,u_SEM, u_wall
end


function create_boundary_conditions(simcase::VelocityBoundaryCase) 
    u_free,u_SEM, u_wall = boundary_velocities(simcase)
    create_boundary_conditions(simcase,u_free,u_SEM, u_wall) 
end

function create_boundary_conditions(simcase::Airfoil, u_free,u_SEM, u_wall) 
    u_diri_tags=["inlet", "limits", "airfoil"]
    u_diri_values = [u_SEM, u_SEM, u_wall]
    p_diri_tags=["outlet"]
    p_diri_values = [0.0]
    return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
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
    D = get_field(simcase,[:D])

    if D == 2
        add_tag_from_tags!(labels, "diri1", [5,])
        add_tag_from_tags!(labels, "diri0", [1, 2, 4, 3, 6, 7, 8])
        add_tag_from_tags!(labels, "p", [4,])
        u_diri_tags=["diri0", "diri1"]
        u_diri_values = [u_wall, u_free]
        p_diri_tags= ["p"]
    
    elseif D == 3
        println("Model 3 dimension")

        add_tag_from_tags!(labels, "diri1", [3,4,7,8,10,12,19,20,24])
        add_tag_from_tags!(labels, "diri0", [1, 2,5,6,9,11,13,14,15,16,17,18,21,22,23,25,26])
        add_tag_from_tags!(labels, "p", [1,])
        u_diri_tags=["diri0", "diri1"]
        u_diri_values = [u_wall, u_free]
        p_diri_tags= ["p"]

    end
  

    p_diri_values = 0.0

    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values
end 


function create_boundary_conditions(simcase::TaylorGreen)

  diameter = 0.5 #0.5 [m] vortex dimension
  Vs = 1 #1[m/s]swirling speed
  Ua = 0.3 #0.3 [m/s]convective velocity in x
  Va = 0.2 #0.2 [m/s]convective velocity in y
  params[:ν] = 0.001 #0.001 m2/s 


  velocity, pa, ωa = analytical_solution(diameter, Vs, Ua, Va, params[:ν])
  u_diri_tags=[]
  u_diri_values = []
  p_diri_tags= ["centre"]
  p_diri_values=pa
  return u_diri_tags,u_diri_values,p_diri_tags,p_diri_values
end
