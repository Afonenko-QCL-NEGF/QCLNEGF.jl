using QCLNEGF
using LinearAlgebra
using Unitful

p = reference_parameters(V_period = 56.0u"mV", Tᴸ = 200u"K")
n = baseline_numerics()
s = ScaleSystem()
g = build_grids(p, n, s)
profiles = build_profiles(p, n, s, g)

println("Lₚ = ", sum(layer.d for layer in p.layers) |> x -> uconvert(u"nm", x))
println("Eₚ = ", uconvert(u"meV", p.F_bias * sum(layer.d for layer in p.layers) * u"eV/V"))
println("Σ wˣ Nᴰ = ", dot(g.wˣ, profiles.Nᴰ), " (dimensionless sheet density)")
println("doped cell indices (Julia) = ", findall(!iszero, profiles.Nᴰ))
