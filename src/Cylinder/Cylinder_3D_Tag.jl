

d = 0.1 #diameter
L_back = 2
L_front = 0.2
H = 0.41

v1(v) = v[1]
v2(v) = v[2]
v3(v) = v[3]


function is_cylinder(x::Vector{VectorValue{3,Float64}})
    distance = v1.(x) .^ 2 + v2.(x) .^ 2
    isapprox(distance[1], (d / 2)^2)
end

function is_inlet(x::Vector{VectorValue{3,Float64}})
    isapprox(v1.(x)[1], -L_front)
end



function is_outlet(x::Vector{VectorValue{3,Float64}})
    isapprox(norm(v1.(x)), L_back)
end


function is_wall(x::Vector{VectorValue{3,Float64}})
    isapprox(norm(v2.(x)), H)
end


function create_3d_tag!(model)
 
    labels = get_face_labeling(model)
    model_nodes = DiscreteModel(Polytope{0}, model)
    cell_nodes_coords = get_cell_coordinates(model_nodes)


    cylinder = collect1d(lazy_map(is_cylinder, cell_nodes_coords))
    cylinder_node = findall(cylinder)
    
    inlet = collect1d(lazy_map(is_inlet, cell_nodes_coords))
    inlet_node = findall(inlet)

    outlet = collect1d(lazy_map(is_outlet, cell_nodes_coords))
    outlet_node = findall(outlet)

    wall = collect1d(lazy_map(is_wall, cell_nodes_coords))
    wall_node = findall(wall)
        cylinder_entity = num_entities(labels) + 1
        inlet_entity = cylinder_entity + 1
        outlet_entity = inlet_entity + 1
        wall_entity = outlet_entity + 1

        for i in cylinder_node
            labels.d_to_dface_to_entity[1][i] = cylinder_entity
        end
        add_tag!(labels, "Cylinder", [cylinder_entity])

        for i in inlet_node
            labels.d_to_dface_to_entity[1][i] = inlet_entity
        end
        add_tag!(labels, "Inlet", [inlet_entity])

        for i in outlet_node
            labels.d_to_dface_to_entity[1][i] = outlet_entity
        end
        add_tag!(labels, "Outlet", [outlet_entity])

        for i in wall_node
            labels.d_to_dface_to_entity[1][i] = wall_entity
        end
        add_tag!(labels, "Walls", [wall_entity])

           
    
    return model

end





function create_3d_tag!(model::GridapDistributed.DistributedDiscreteModel)
 
    
    map_parts(model.models) do model
        create_3d_tag!(model)
    end
return model
end
