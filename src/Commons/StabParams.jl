
function h_param(Ω::GridapDistributed.DistributedTriangulation, D::Int64)
    h = map(Ω.trians) do trian
        h_param(trian, D)
    end
    h = CellData.CellField(h, Ω)
    h
end


function h_param(Ω::Triangulation, D::Int64)
    h = lazy_map(h -> h^(1 / D), get_cell_measure(Ω))

    h
end



function G_params(Ω::GridapDistributed.DistributedTriangulation, params)
    G = map(Ω.trians) do trian
      compute_G(trian,params)
    end
    GG = map(Ω.trians) do trian
      compute_GG(trian,params)
    end
    gg = map(Ω.trians) do trian
      compute_gg(trian,params)
    end

    G = CellData.CellField(G, Ω)
    GG = CellData.CellField(GG, Ω)
    gg = CellData.CellField(gg, Ω)
    return G, GG, gg
  end

  function G_params(Ω::Triangulation, params)
    G = compute_G(Ω,params)
    GG = compute_GG(Ω,params)
    gg = compute_gg(Ω,params)
    return G, GG, gg
end
  
  function compute_d(trian::Gridap.Geometry.BodyFittedTriangulation, params) #trian == Ω
    D = params[:D]
  
    ξₖ = get_cell_map(trian)
    Jt = lazy_map(Broadcasting(∇), ξₖ)
    inv_Jt = lazy_map(Operation(inv), Jt)
  
    if D == 2
      eval_point = Point(0.5, 0.5)
    else
      eval_point = Point(0.5, 0.5, 0.5)
    end
  
    d = lazy_map(evaluate, inv_Jt, Fill(eval_point, num_cells(trian)))
    return d
  end

  function compute_G(trian::Gridap.Geometry.BodyFittedTriangulation, params) #trian == Ω
    d = compute_d(trian, params) #trian == Ω
    dtrans = lazy_map(Broadcasting(transpose), d)
    G = lazy_map(Broadcasting(⋅), d, dtrans)
    return G
  end

  function compute_GG(trian::Gridap.Geometry.BodyFittedTriangulation, params) #trian == Ω
    G = compute_G(trian, params) #trian == Ω
    GG = lazy_map(Broadcasting(⊙), G, G)
    return GG
  end

  function compute_gg(trian::Gridap.Geometry.BodyFittedTriangulation, params) #trian == Ω
    @unpack D = params

    d = compute_d(trian, params) #trian == Ω

    function gg_operation(d)
  
      if D == 2
  
        return (d[1] + d[3])^2 + (d[2] + d[4])^2
  
  
      elseif D == 3
  
        return (d[1] + d[4] + d[7])^2 + (d[2] + d[5] + d[8])^2 + (d[3] + d[6] + d[9])^2
  
      end
  
    end
  
    gg = lazy_map(Broadcasting(gg_operation), d)
    return gg
  end