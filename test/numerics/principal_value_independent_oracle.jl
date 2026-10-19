module PrincipalValueIndependentOracle
include("../support/common.jl")

const ATOL = 2e-13
const RTOL = 2e-13
const OPERATORS = (direct_hilbert_transform, QCLNEGF.fft_hilbert_transform)


# O08-06 ratio propagation: retain the scalar functional budget, then enclose
# its quotient. Expected errors/endpoints never depend on operator output.
strict_error_gate(x, a) = isfinite(x) && isfinite(a) && 0 < x < 1 &&
    0 < a < 1 && isapprox(x, a; atol=ATOL, rtol=RTOL)

function positive_ratio_interval(a, b)
    B = nextfloat(ATOL)
    alo, ahi = prevfloat(a-B), nextfloat(a+B)
    blo, bhi = prevfloat(b-B), nextfloat(b+B)
    finite_endpoints = all(isfinite, (alo, ahi, blo, bhi))
    positive_endpoints = 0 < alo <= ahi < 1 && 0 < blo <= bhi < 1
    @test finite_endpoints
    @test positive_endpoints
    finite_endpoints && positive_endpoints || error(
        "DEV08 STOP: unresolved ratio interval; finite positive endpoints required",
    )
    lo, hi = prevfloat(alo/bhi), nextfloat(ahi/blo)
    finite_positive_ratio = isfinite(lo) && isfinite(hi) && 0 < lo <= hi
    @test finite_positive_ratio
    finite_positive_ratio || error("DEV08 STOP: invalid outward quotient interval")
    return (; lo, hi, alo, ahi, blo, bhi)
end

@testset "Independent literal omitted-diagonal direct and FFT PV map" begin
    energy = [-1.0, 0.0, 1.0]
    weights = [1/4, 1/2, 3/4]
    gamma = [1.0, 2.0, 4.0]
    # Exactly two independently expanded scalar terms per row; no map helper.
    expected = [-5/(4π), -11/(8π), 9/(16π)]
    @test expected ≈ [(-1-3/2)/(2π), (1/4-3)/(2π), (1/8+1)/(2π)] atol=ATOL rtol=RTOL
    for operator in OPERATORS
        actual = vec(operator(reshape(gamma, 3, 1, 1, 1), energy, weights))
        @test all(isfinite, actual)
        @test actual ≈ expected atol=ATOL rtol=RTOL
        @test operator(zeros(ComplexF64, 3, 1, 1, 1), energy, weights) == zeros(ComplexF64, 3, 1, 1, 1)
    end
    # Deliberate operator variants keep the same literal RHS. Diagonal cell
    # addition is finite here; a divide-by-zero variant also fails finiteness.
    wrong_sign = -expected
    inserted_diagonal = expected .+ weights .* gamma ./ (2π)
    cyclic = zeros(3)
    for i in 1:3, j in 1:3
        i == j && continue
        # Incorrect circular nearest-index distance, rather than i-j.
        distance = mod(i-j+1, 3)-1
        cyclic[i] += weights[j]*gamma[j]/(2π*distance)
    end
    for mutated in (wrong_sign, inserted_diagonal, cyclic)
        @test !isapprox(mutated, expected; atol=ATOL, rtol=RTOL)
        @test !isapprox(mutated[1], expected[1]; atol=ATOL, rtol=RTOL)
    end
end

@testset "Center-only continuum PV separates omitted cell and finite support" begin
    # Gamma(x)=4+x on [-3,3], zero outside. Not a physical SCBA spectrum.
    # PV at E=0 is -min(W,3)/pi; the omitted central cell adds h/(2pi).
    # h=2/N, boundaries +/-3 are cell edges for every odd N below.
    full_continuum = -3/π
    for operator in OPERATORS
        for N in (17, 33, 65)
            h = 2/N
            centres = Dict{Int,Float64}()
            for W in (1, 3, 5)
                count = W*N
                energy = [(i-(count+1)/2)*h for i in 1:count]
                weights = fill(h, count)
                gamma = [abs(x)<3 ? 4+x : 0.0 for x in energy]
                actual = operator(reshape(gamma, count, 1, 1, 1), energy, weights)
                centre = real(actual[(count+1)÷2, 1, 1, 1])
                @test all(isfinite, actual)
                finite_continuum = -min(W, 3)/π
                omitted_cell = h/(2π)
                omitted_tail = (3-min(W, 3))/π
                @test centre ≈ finite_continuum+omitted_cell atol=ATOL rtol=RTOL
                @test centre-finite_continuum ≈ omitted_cell atol=ATOL rtol=RTOL
                @test centre-full_continuum ≈ omitted_tail+omitted_cell atol=ATOL rtol=RTOL
                @test centre-finite_continuum > 0
                centres[W] = centre
            end
            @test centres[1]-centres[3] ≈ 2/π atol=ATOL rtol=RTOL
            @test centres[3] ≈ centres[5] atol=ATOL rtol=RTOL
        end
        # Exact center error has first-order omitted-cell behavior. Compare
        # measured errors against the analytic h ratio, not an assumed order.
        errors = Float64[]
        for N in (17, 33, 65)
            h = 2/N
            energy = [(i-(3N+1)/2)*h for i in 1:3N]
            gamma = 4 .+ energy
            actual = operator(reshape(gamma, 3N, 1, 1, 1), energy, fill(h, 3N))
            push!(errors, real(actual[(3N+1)÷2, 1, 1, 1])-full_continuum)
        end
        nodes = (17, 33, 65)
        expected_errors = [(2/N)/(2π) for N in nodes]
        # These are separately recomputed operator calls. Earlier window
        # checks do not establish this loop's strict per-error premises.
        for i in eachindex(errors)
            @test isfinite(errors[i]) && 0 < errors[i] < 1
            @test 0 < expected_errors[i] < 1
            @test isapprox(errors[i], expected_errors[i]; atol=ATOL, rtol=RTOL)
            strict_error_gate(errors[i], expected_errors[i]) || error(
                "DEV08 STOP: strict per-error functional precondition failed",
            )
        end
        @test errors[1] > errors[2] > errors[3] > 0
        for (i, j) in ((1, 2), (2, 3))
            a, b = expected_errors[i], expected_errors[j]
            box = positive_ratio_interval(a, b)
            @test box.bhi < box.alo
            rho = nodes[j]/nodes[i]
            @test box.lo <= rho <= box.hi
            measured_ratio = errors[i]/errors[j]
            @test isfinite(measured_ratio)
            @test box.lo <= measured_ratio <= box.hi

            # Scalar controls add no operator evaluations. Ratio evidence
            # complements, and never replaces, signed amplitude evidence.
            @test strict_error_gate(a, a) && strict_error_gate(b, b)
            @test box.lo <= a/b <= box.hi
            x, y = a+ATOL/2, b-ATOL/2
            @test strict_error_gate(x, a) && strict_error_gate(y, b)
            @test box.lo <= x/y <= box.hi
            constant_x, constant_y = b, b
            @test !(box.lo <= constant_x/constant_y <= box.hi)
            second_x, second_y = rho^2*b, b
            @test !(box.lo <= second_x/second_y <= box.hi)
            scaled_x, scaled_y = 2a, 2b
            @test box.lo <= scaled_x/scaled_y <= box.hi
            @test !strict_error_gate(scaled_x, a) && !strict_error_gate(scaled_y, b)
            negative_x, negative_y = -a, -b
            @test box.lo <= negative_x/negative_y <= box.hi
            @test !(0 < negative_x < 1) && !(0 < negative_y < 1)
            @test !isapprox(negative_x, a; atol=ATOL, rtol=RTOL)
            @test !isapprox(negative_y, b; atol=ATOL, rtol=RTOL)
        end
    end
end
end
