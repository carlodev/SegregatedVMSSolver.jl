function creation_fe_spaces(simcase::TaylorGreen{Periodic}, model, boundary_conditions)

    u_diri_tags,u_diri_values,p_diri_tags,p_diri_values = boundary_conditions

    @sunpack D,order = simcase
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D, Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order)


    V = TestFESpace(model, reffeᵤ, conformity=:H1)
    U = TransientTrialFESpace(V)
    if D == 2
        Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
        P = TransientTrialFESpace(Q, p_diri_values)
    elseif D ==3
        Q = TestFESpace(model, reffeₚ, conformity=:H1)
        P = TransientTrialFESpace(Q)        
    end

    
  
    return V, U, P, Q

end


function creation_fe_spaces(simcase::TaylorGreen{Natural}, model, boundary_conditions)

    u_diri_tags,u_diri_values,p_diri_tags,p_diri_values = boundary_conditions

    @sunpack D,order = simcase
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D, Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order)


    V = TestFESpace(model, reffeᵤ, conformity=:H1)
    U = TransientTrialFESpace(V)

    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
    P = TransientTrialFESpace(Q, p_diri_values)

    return V, U, P, Q

end


function creation_fe_spaces(simcase, model, boundary_conditions)

    u_diri_tags,u_diri_values,p_diri_tags,p_diri_values = boundary_conditions
    @sunpack D,order = simcase

    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    reffeₚ = ReferenceFE(lagrangian, Float64, order)


    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=u_diri_tags)
    U = TransientTrialFESpace(V, u_diri_values)

    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=p_diri_tags)
    P = TrialFESpace(Q, p_diri_values)

    return V, U, P, Q

end

