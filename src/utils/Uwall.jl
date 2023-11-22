G = 1.12
H = 0.3504
N = 90
Ps = 0.0
for i = 0:1:(N-1)
  Ps = Ps + G^i
end
h0 = H/Ps

# h = 3.83e-6
τw = μ .* abs.(Friction)
uτ = sqrt.(τw)

yplus = h0 .* uτ ./ μ

plot(1:1:177, yplus.p)
