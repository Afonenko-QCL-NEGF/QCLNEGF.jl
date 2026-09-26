"""Versioned backward-error budget for a PSD test in scaled solver units.

For an n×n Float64 Hermitian eigensolve, 64 n eps bounds the declared
construction/evaluation roundoff allowance at max(block norm, construction
scale). This is a numerical noise allowance, not a claim that all negative
values are roundoff. `construction_scale=1` is one solver unit; callers testing
a rescaled operator must transform it with the matrix. No clipping occurs.
"""
const PSD_METRIC_VERSION = 3
const _PSD_BACKWARD_ERROR_FACTOR = 64.0

function psd_error_budget(
    block_norm::Real,
    dimension::Integer;
    construction_scale::Real = 1.0,
    construction_error::Real = 0.0,
)
    dimension >= 1 || throw(ArgumentError("PSD dimension must be positive"))
    all(x -> isfinite(x) && x >= 0, (block_norm, construction_scale, construction_error)) ||
        throw(ArgumentError("PSD scales and error must be finite and nonnegative"))
    return Float64(
        construction_error +
        _PSD_BACKWARD_ERROR_FACTOR *
        dimension *
        eps(Float64) *
        max(block_norm, construction_scale),
    )
end

"""PSD gate witness; ratio is the remaining defect after a roundoff budget.

Acceptance r_PSD < τ is equivalent to absolute_defect < backward_error +
τ block_norm. Physical dimensions are A/G: 1/E₀, broadening: E₀. Coordinates
are indices in the owning solution's declared basis/discretization.
"""
struct _PositivityWitness
    matrix_kind::Symbol
    energy_index::Int
    momentum_index::Int
    minimum_eigenvalue::Float64
    block_norm::Float64
    floor::Float64
    ratio::Float64
    absolute_defect::Float64
    relative_defect::Float64
    backward_error::Float64
    hermiticity_defect::Float64
    construction_scale::Float64
end

# Legacy witness readers still accept the seven original fields.
_PositivityWitness(kind, e, m, minimum_value, block_norm, floor, ratio) =
    _PositivityWitness(
        kind,
        e,
        m,
        minimum_value,
        block_norm,
        floor,
        ratio,
        max(0.0, -minimum_value),
        max(0.0, -minimum_value)/max(block_norm, floatmin(Float64)),
        floor,
        0.0,
        1.0,
    )
_empty_positivity_witness() =
    _PositivityWitness(:none, 0, 0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)

function _positivity_witness(
    kind,
    e,
    m,
    values;
    floor = 1e-14,
    construction_scale = 1.0,
    construction_error = 0.0,
    hermiticity_defect = 0.0,
)
    all(isfinite, values) || return _PositivityWitness(kind, e, m, NaN, Inf, floor, Inf)
    minimum_eigenvalue = minimum(values)
    block_norm = maximum(abs, values)
    absolute = max(0.0, -minimum_eigenvalue)
    budget =
        psd_error_budget(block_norm, length(values); construction_scale, construction_error)
    denominator = max(block_norm, floatmin(Float64))
    ratio = max(0.0, absolute-budget)/denominator
    return _PositivityWitness(
        kind,
        e,
        m,
        minimum_eigenvalue,
        block_norm,
        budget,
        ratio,
        absolute,
        absolute/denominator,
        budget,
        Float64(hermiticity_defect),
        Float64(construction_scale),
    )
end

"""Worst PSD blocks on the complete grid, including raw Dyson/Keldysh products.

Raw occupation and hole matrices are checked independently of any charge
constraint update. This prevents normalization from hiding inadmissible input.
"""
function _green_positivity_witnesses(green::GreenState, total::SelfEnergyFamily)
    return _scan_positivity_witnesses(
        green,
        6,
        (e, m)->begin
            GR=_matrix_block(green.Gᴿ, e, m)
            lesser=_matrix_block(total.Σˡ, e, m)
            greater=_matrix_block(total.Σᵍ, e, m)
            _positivity_blocks(green, e, m, lesser, greater, GR)
        end,
    )
end

function _positivity_blocks(green, e, m, lesser, greater, GR)
    return (
        (:spectral, _matrix_block(green.A, e, m)),
        (:occupied, -im .* _matrix_block(green.Gˡ, e, m)),
        (:unoccupied, im .* _matrix_block(green.Gᵍ, e, m)),
        (:broadening, im .* (greater-lesser)),
        (:raw_occupied, GR*(-im .* lesser)*GR'),
        (:raw_unoccupied, GR*(im .* greater)*GR'),
    )
end

# A caller supplying only broadening has no measured raw Keldysh products.
# Preserve this four-kind interface instead of inventing raw matrices.
function _green_positivity_witnesses(green::GreenState, broadening_block::Function)
    return _scan_positivity_witnesses(
        green,
        4,
        (e, m)->(
            (:spectral, _matrix_block(green.A, e, m)),
            (:occupied, -im .* _matrix_block(green.Gˡ, e, m)),
            (:unoccupied, im .* _matrix_block(green.Gᵍ, e, m)),
            (:broadening, broadening_block(e, m)),
        ),
    )
end

function _scan_positivity_witnesses(green, count, blocks_at)
    witnesses=fill(_empty_positivity_witness(), count)
    for e in axes(green.A, 1), m in axes(green.A, 2)
        for (index, (kind, block)) in enumerate(blocks_at(e, m))
            hermitian=(block+block')/2
            values=if all(isfinite, hermitian)
                try
                    eigvals(Hermitian(hermitian))
                catch error
                    error isa LinearAlgebra.LAPACKException || rethrow()
                    [NaN]
                end
            else
                [NaN]
            end
            witness=_positivity_witness(
                kind,
                e,
                m,
                values;
                hermiticity_defect = all(isfinite, block) ? norm(block-block')/2 : Inf,
            )
            witness.ratio>witnesses[index].ratio && (witnesses[index]=witness)
        end
    end
    return witnesses
end

# Only one local self-energy block is assembled: no duplicate full-state array.
function _scba_positivity_blocks(scba::SCBAResult, mechanism_order, e, m)
    lesser=_matrix_block(scba.embedding.Σˡ, e, m)
    greater=_matrix_block(scba.embedding.Σᵍ, e, m)
    for name in mechanism_order
        family=scba.scattering[name]
        lesser .+= view(family.Σˡ, e, m, :, :)
        greater .+= view(family.Σᵍ, e, m, :, :)
    end
    return _positivity_blocks(
        scba.green,
        e,
        m,
        lesser,
        greater,
        _matrix_block(scba.green.Gᴿ, e, m),
    )
end

function _scba_positivity_witnesses(scba::SCBAResult, mechanism_order)
    return _scan_positivity_witnesses(
        scba.green,
        6,
        (e, m)->_scba_positivity_blocks(scba, mechanism_order, e, m),
    )
end
