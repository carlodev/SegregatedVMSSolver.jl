#It is a unique case, with pressure that is time dependent and no boundary conditions on the velocity, it is periodic in all dimensions
function CreateTGSpaces(model, params, pa)
    
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D], Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(model, reffeᵤ, conformity=:H1)
    U = TransientTrialFESpace(V)

    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags="centre")
    P = TransientTrialFESpace(Q, pa)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])
    
    return V, Q, U, P, Y, X, model
end

#It is a unique case, with pressure that is time dependent and no boundary conditions on the velocity, it is periodic in all dimensions
function CreateTGSpaces_va(model, params, va)
    
    model = add_centre_tag!(model, Point(0.0, 0.0)) #(0.0, 0.0) is the centre coordinate
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{params[:D], Float64}, params[:order])
    reffeₚ = ReferenceFE(lagrangian, Float64, params[:order])


    V = TestFESpace(model, reffeᵤ, conformity=:H1,dirichlet_tags="centre")
    U = TransientTrialFESpace(V, va)

    Q = TestFESpace(model, reffeₚ, conformity=:H1)
    P = TransientTrialFESpace(Q)

    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])
    
    return V, Q, U, P, Y, X, model
end