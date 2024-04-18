function create_model(parts, simcase::MeshFileCase)
    model = GmshDiscreteModel(parts, get_field(simcase,:meshfile))
    print_model(model,simcase)
    return model
end


function create_model(simcase::CartesianCase)
    D,N = get_field(simcase,[:D,:N])
    L = 0.5 #same length in all the 3 dimensions
    if D == 2
        domain = (-L, L, -L, L)
        partition = (N,N)
    elseif D == 3
        domain = (-L, L, -L, L, -L, L)
        partition = (N,N,N)
    end

    return domain,partition

end


function create_model(parts, simcase::LidDriven)
    domain,partition,rank_partition = create_model(simcase)

    function stretching_y_function(x)
        gamma1 = 2.5
        S = 0.5815356159649889 #for rescaling the function over the domain -0.5 -> 0.5
        -tanh.(gamma1 .* (x)) ./ tanh.(gamma1) .* S
    end

    
    function stretching(x::Point)
        m = zeros(length(x))
        m[1] = stretching_y_function(x[1])
        m[2] = stretching_y_function(x[2])
        if length(x)>2
            m[3] = stretching_y_function(x[3])
        end
        Point(m)
    end
    model =CartesianDiscreteModel(parts,rank_partition, domain, partition, map=stretching)

    return model 
end




function create_model(parts, simcase::TaylorGreen)
    domain,partition,rank_partition = create_model(simcase)

    model =CartesianDiscreteModel(parts,rank_partition, domain, partition;isperiodic=(true, true) )
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate

    return model 
end