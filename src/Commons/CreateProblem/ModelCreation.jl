function create_model(parts, simcase::VelocityBoundaryCase)
    mesh = simcase.meshp
    model = create_model(parts, mesh.meshinfo, mesh.D, mesh.rank_partition)
    print_model(model,simcase)
    return model
end



function create_model(parts, meshinfo::GmshMeshParams, D, rank_partition)
    model = GmshDiscreteModel(parts, meshinfo.filename)
    return model
end


function create_model(meshinfo::CartesianMeshParams, D::Int64)
    N = meshinfo.N
    L = meshinfo.L
    @assert length(L) == D
    if D == 2
        domain = (-L[1], L[1], -L[2], L[2])
        partition = (N[1],N[2])
    elseif D == 3
        domain = (-L[1], L[1], -L[2], L[2], -L[3], L[3])
        partition = (N[1],N[2],N[3])
    end

    return domain,partition

end


function create_model(parts, meshinfo::CartesianMeshParams, D::Int64, rank_partition)

    domain,partition = create_model(meshinfo,D)

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
    labels = get_face_labeling(model)

    if D == 2
        add_tag_from_tags!(labels, "diri1", [5,])
        add_tag_from_tags!(labels, "diri0", [1, 2, 4, 3, 6, 7, 8])
        add_tag_from_tags!(labels, "p", [4,])
    elseif D == 3
        add_tag_from_tags!(model, "diri1", [3,4,7,8,10,12,19,20,24])
        add_tag_from_tags!(model, "diri0", [1, 2,5,6,9,11,13,14,15,16,17,18,21,22,23,25,26])
        add_tag_from_tags!(model, "p", [1,])
    end

    return model 
end


function create_model(parts, simcase::TaylorGreen)
    mesh = simcase.meshp
    rank_partition = mesh.rank_partition

    domain,partition = create_model(mesh.meshinfo, mesh.D)

    model =CartesianDiscreteModel(parts,rank_partition, domain, partition;isperiodic=(true, true) )
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate
    
    print_model(model,simcase)
    return model 
end


"""
  print_model(model,simcase::SimulationCase)

It prints the model mesh
"""
function print_model(model,simcase::SimulationCase)
  @sunpack printmodel = simcase
  if printmodel
      writevtk(model,"$(typeof(simcase))")
  end
end

