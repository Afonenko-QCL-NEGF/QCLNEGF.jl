using QCLNEGF
using LinearAlgebra

p = reference_parameters()
n = tutorial_numerics()
s = ScaleSystem()
g = build_grids(p, n, s)
profiles = build_profiles(p, n, s, g)
H = build_bdd_hamiltonian(profiles, g, s)
basis = build_localized_basis(p, n, s, g, profiles)

println("BDD shape: ", size(H))
println("basis shape: ", size(basis.Φ))
println("r_orth = ", basis.r_orth)
println("κ₂(S_cell) = ", basis.κ_overlap)
println("||T₋-T₊†|| = ", norm(basis.T₋ - basis.T₊'))
