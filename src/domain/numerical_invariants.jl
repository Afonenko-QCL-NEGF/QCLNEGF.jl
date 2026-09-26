"""Versioned structural gates for representation sanity, not run tolerances."""
const _EFFECTIVE_HAMILTONIAN_NORM_FLOOR = 1.0e-14
const _EFFECTIVE_HAMILTONIAN_HERMITICITY_LIMIT = 1.0e-12
const _OPTICAL_POSITION_HERMITICITY_LIMIT = 1.0e-11
const _NUMBER_FUNCTIONAL_REALNESS_LIMIT = 1.0e-10

# Exact early-return guard for the representable tails of the Fermi function.
# It is an implementation overflow invariant, not a physical truncation: at
# this exponent the discarded occupation is below 4.3e-18.
const _FERMI_EXPONENT_SATURATION = 40.0
