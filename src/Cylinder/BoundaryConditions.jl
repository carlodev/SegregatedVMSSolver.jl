function bc_cylinder(params)
    @unpack D, u_in, t_endramp, TI = params

    #u_free(x,t) = VectorValue(4*umax/0.41^2*(0.205^2 - x[2]^2),0) #parabolic inlet

    # u_free(x,t) = D == 2 ? VectorValue(u_in, 0.0) :  VectorValue(u_in, 0.0, 0.0)
    # u_free(t::Real) = x -> u_free(x,t)

    uin(t) = (t < t_endramp) ? (u_in .*(t/t_endramp)) : u_in

    u_free(x, t) = (D == 2) ? VectorValue(uin(t), 0.0) :  VectorValue(uin(t), 0.0, 0.0)
    u_free(t::Real) = x -> u_free(x,t)

    u_wall(x,t) = D == 2 ? VectorValue(0.0, 0.0) :  VectorValue(0.0, 0.0, 0.0)
    u_wall(t::Real) = x -> u_wall(x,t)

  
    u_diri_tags=["inlet", "limits", "cylinder"]
    u_diri_values = [u_free, u_wall, u_wall]
    p_diri_tags=["outlet"]
    p_diri_values = [0.0]
    force_tags = ["cylinder"]

        
    return u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u_free, force_tags


end