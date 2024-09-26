
"""
    analytical_solution(diameter::Int64, Vs::Float64, Ua::Float64, Va::Float64, ν::Float64) 

It provides the anlytical solution for the Taylor Green Vortex case. Solution for the 2D periodic case
"""
function analytical_solution(diameter, Vs, Ua, Va, ν) 
  
  Tx(x, t) = pi / diameter * (x[1] - Ua * t)
  Ty(x, t) = pi / diameter * (x[2] - Va * t)
  Et(t) = exp(-(2 * ν * t * pi^2) / (diameter^2))
  ua(x, t) = Ua - Vs * cos(Tx(x, t)) * sin(Ty(x, t)) * Et(t)
  va(x, t) = Va + Vs * sin(Tx(x, t)) * cos(Ty(x, t)) * Et(t)
  velocity(x, t) = VectorValue(ua(x, t), va(x, t))
  pa(x, t) = -(Vs^2 / 4) * (cos(2 * Tx(x, t)) + cos(2 * Ty(x, t))) * Et(t)^2
  ωa(x, t) = 2 * Vs * pi / diameter * cos(Tx(x, t)) * cos(Ty(x, t)) * Et(t)^2

  ω₀ = 2 * Vs * pi / diameter

  ua(t::Real) = x -> ua(x, t)
  va(t::Real) = x -> va(x, t)
  velocity(t::Real) = x -> velocity(x, t)
  pa(t::Real) = x -> pa(x, t)
  ωa(t::Real) = x -> ωa(x, t)

  velocity, pa, ωa
end


function TGV_initial(u0::Float64, D::Int64, L::Vector)

  Npi(x::VectorValue,i::Int64) = x[i] * (pi)/(2*L[i])

  ux(x,t) = (D==2) ? u0*cos(Npi(x,1))*sin(Npi(x,2)) : u0*cos(Npi(x,1))*sin(Npi(x,2))*cos(Npi(x,3)) 
  uy(x,t) = (D==2) ? -u0*sin(Npi(x,1))*cos(Npi(x,2)) : -u0*sin(Npi(x,1))*cos(Npi(x,2))*cos(Npi(x,3)) 
  uz(x,t) = 0.0
  velocity(x,t) = (D==2) ? VectorValue(ux(x,t), uy(x,t)) :  VectorValue(ux(x,t), uy(x,t), uz(x,t))
  p0(x,t) =  0.0  # (D==2) ? 1/16 * (cos(2*Npi(x,1)) + cos(2*Npi(x,2)))*(2) : 1/16 * (cos(2*Npi(x,1)) + cos(2*Npi(x,2)))*(cos(2*Npi(x,3))+2)

  ux(t::Real) = x -> ux(x, t)
  vy(t::Real) = x -> uy(x, t)
  uz(t::Real) = x -> uy(x, t)

  velocity(t::Real) = x -> velocity(x, t)
  p0(t::Real) = x -> p0(x, t)

  return velocity, p0
end