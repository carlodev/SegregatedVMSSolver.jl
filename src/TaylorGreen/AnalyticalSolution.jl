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

