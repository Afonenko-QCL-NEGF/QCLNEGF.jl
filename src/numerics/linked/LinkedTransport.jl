"""Finite-chain Q-linked transport with packed local collision sources."""
module LinkedTransport

using LinearAlgebra
using SparseArrays: SparseMatrixCSC, sparse
using ..QCLDomain: SelfEnergyFamily
import ..QCLNumerics as N

export LinkedProblem, LocalCollisionKernels, linked_fresh_map,
    project_linked, project_linked_economical, COLLISION_CHANNELS

# Numerical stages carried over from the validated linked-map-9de1c83 release.
include("tail.jl")
include("causal.jl")
include("collision.jl")
include("quadrature.jl")
include("selected_rows.jl")
include("projection.jl")
include("map.jl")

end
