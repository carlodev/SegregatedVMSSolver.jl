include("BoundaryConditions.jl")



function run_liddriven(params,distribute)
    @unpack rank_partition, N, D, order = params

    parts  = distribute(LinearIndices((prod(rank_partition),)))

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

        L = 0.5
        if D == 2
            domain = (-L, L, -L, L)
            partition = (N,N)
        elseif D == 3
            domain = (-L, L, -L, L, -L, L)
            partition = (N,N,N)
        end

    
        model = CartesianDiscreteModel(parts,rank_partition, domain, partition, map=stretching)
        
  


    u_diri_tags, u_diri_values, p_diri_tags, p_diri_values, u0 = bc_liddriven(model, params) #u_top in bc_liddriven taken as initial velocity
    merge!(params, Dict(:u0 => u0, :model => model, :force_tags=>nothing, :parts=>parts))
    print_model(params)
    
    
     V, U, P, Q, Y, X = creation_fe_spaces(params, u_diri_tags, u_diri_values, p_diri_tags, p_diri_values)
    println("spaces created")

    degree = 4*order
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    new_dict = Dict(:U => U,
                    :P => P,
                    :X => X,
                    :Y => Y,
                    :Ω => Ω,
                    :dΩ => dΩ,
                    :degree => degree,
                    :force_params => nothing)
    merge!(params, new_dict)

    trials = [U, P]
    tests = [V, Q]
  
    merge!(params, Dict(:trials => trials, :tests => tests))
    
    solve_case(params)


   

    

  

end