"""
    create_new_tag!(model::GridapDistributed.DistributedDiscreteModel, tagname::String, is_tag::Function)

It creates the `centre` tag at the tag_coordinate (Point); if mesh extremely fine the tolrances have to be smaller (unlikely)
"""
function create_new_tag!(model::GridapDistributed.DistributedDiscreteModel, tagname::String, is_tag::Function)

    map(model.models) do model
       create_new_tag!(model,tagname::String,is_tag)
       end
println("New tag $tagname added to the model")
end

function create_new_tag!(model,tagname::String, is_tag::Function)
       labels = get_face_labeling(model)
       new_entity = num_entities(labels) + 1

       Dmax = length(labels.d_to_dface_to_entity)

    for D = 1:1:Dmax
       model_nodes = DiscreteModel(Polytope{D-1}, model)
       cell_nodes_coords = get_cell_coordinates(model_nodes)
       cell_node_centre = collect1d(lazy_map(is_tag, cell_nodes_coords))
       cell_node = findall(cell_node_centre)
       

       for centre_point in cell_node
           labels.d_to_dface_to_entity[D][centre_point] = new_entity
       end

    end

       add_tag!(labels, tagname, [new_entity])
 
end


"""
    add_new_tag!(model, tag_coordinate::Point, tagname::String)

It add a new tag named `tagname` at the specified `tag_coordinates`    
"""
function add_new_tag!(model, tag_coordinate::Point, tagname::String)
    
  
a_tol = 1e-8
    function is_tag(x::Vector{VectorValue{2,Float64}})
        isapprox(norm(getindex.(x,1)), tag_coordinate[1], atol=a_tol) && isapprox(norm( getindex.(x,2) ), tag_coordinate[2], atol=a_tol)
    end

    function is_tag(x::Vector{VectorValue{3,Float64}})
        isapprox(norm(getindex.(x,1)), tag_coordinate[1], atol=a_tol) && isapprox(norm(getindex.(x,1)), tag_coordinate[2], atol=a_tol) && isapprox(norm(getindex.(x,3)), tag_coordinate[3], atol=a_tol)
    end

    create_new_tag!(model,tagname,is_tag)

end


function add_new_tag!(model, range_coordinate::Tuple, tagname::String)
   function is_tag(x::Vector{VectorValue{2,Float64}})
        @assert length(range_coordinate) == 2 "Range new tags not the same dimension of the problem"
        range_x = getindex(range_coordinate,1)
        range_y = getindex(range_coordinate,2)
       x = getindex(x,1)

        return getindex.(x,1) .> range_x[1] && getindex.(x,1) .< range_x[2] && 
               getindex.(x,2) .> range_y[1] && getindex.(x,2) .< range_y[2]
    end

    function is_tag(x::Vector{VectorValue{3,Float64}})
        @assert length(range_coordinate) == 3 "Range new tags not the same dimension of the problem"
        range_x = getindex(range_coordinate,1)
        range_y = getindex(range_coordinate,2)
        range_z = getindex(range_coordinate,3)
        x = getindex(x,1)

        return getindex.(x,1) .> range_x[1] && getindex.(x,1) .< range_x[2] && 
               getindex.(x,2) .> range_y[1] && getindex.(x,2) .< range_y[2] && 
               getindex.(x,3) .> range_z[1] && getindex.(x,3) .< range_z[2]
    end


    create_new_tag!(model,tagname,is_tag)

end

function add_vertical_tag!(model, x_coord::Real, tagname::String)
    
    a_tol = 1/400

    function is_tag(x::Vector)
         x = getindex(x,1)
         return isapprox(getindex.(x,1), x_coord, atol=a_tol) && getindex.(x,2) .< 0.1 && getindex.(x,2) .> 0.0
     end

     create_new_tag!(model,tagname,is_tag)
  end



function add_new_tag!(model, params)
@unpack newtag = params
    if !isnothing(newtag)
        @unpack range_coordinate, tagname = newtag
        comm = MPI.COMM_WORLD

        for (xr,tg) in zip(range_coordinate,tagname)
        # add_new_tag!(model, range_coordinate, tagname)
        println("Xc = $xr, Tag = $tg")
            add_vertical_tag!(model, xr, tg)
            MPI.Barrier(comm)
        end
        
    end
end



"""
    add_centre_tag!(model, tag_coordinate::Point)

It creates the `centre` tag at the tag_coordinate (Point); if mesh extremely fine the tolrances have to be smaller (unlikely)
"""
function add_centre_tag!(model, tag_coordinate::Point)

    add_new_tag!(model, tag_coordinate, "centre")
    return model
end




"add the SEM tag to the nodes in front of the airfoil. c: chord of the airfoil, the points are on a circle at 1 chord of distance from the leading edge. 
set a_tol for selecting only one lines of points, it depends on how fine the mesh is. "
function add_SEM_tag!(model; c=1.0, a_tol = 1e-1)

    v1(v) = v[1]
    v2(v) = v[2]
 
    function is_SEM(x)
        x = x[1]
        r2 = (v1.(x)).^2 .+ (v2.(x)).^2
        x[1][1]< -c/3 && isapprox(r2, c.^2, atol=a_tol)
    end

  
    function create_SEM_tag!(model::GridapDistributed.DistributedDiscreteModel)
        map_parts(model.models) do model
            create_SEM_tag!(model)
        end

    end

    function create_SEM_tag!(model)
        labels = get_face_labeling(model)
        model_nodes = DiscreteModel(Polytope{0}, model)
        cell_nodes_coords = get_cell_coordinates(model_nodes)
        cell_node_centre = collect1d(lazy_map(is_SEM, cell_nodes_coords))
        cell_node = findall(cell_node_centre)
        new_entity = num_entities(labels) + 1
        for centre_point in cell_node
            labels.d_to_dface_to_entity[1][centre_point] = new_entity
        end
        add_tag!(labels, "sem", [new_entity])
    end


    create_SEM_tag!(model)

    return model
end

