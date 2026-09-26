using QCLNEGF

ε = collect(range(-1.0, 1.0; length = 401))
h = reshape(ComplexF64[0.15], 1, 1, 1)
Σᴿ = fill(-0.025im, length(ε), 1, 1, 1)
Gᴿ, κD = retarded_green(ε, h, Σᴿ)
A = spectral_function(Gᴿ)

println("max A = ", maximum(real.(A[:, 1, 1, 1])))
println("max κ₂(D) = ", maximum(κD))
