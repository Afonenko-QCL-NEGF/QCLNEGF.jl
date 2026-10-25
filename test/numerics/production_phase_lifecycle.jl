module ProductionPhaseLifecycleTests

using Test
using Logging

@test VERSION == v"1.13.0"
include(joinpath(@__DIR__, "..", "..", "src", "numerics", "optimized", "production_phase_lifecycle.jl"))

mutable struct Sentinel <: Exception
    tag::Symbol
end

mutable struct Observer
    events::Vector{Tuple{Symbol,Symbol}}
    stack::Vector{Symbol}
    attempts::Vector{Tuple{Symbol,Symbol}}
    failure::Any
end
Observer(failure = nothing) = Observer(Tuple{Symbol,Symbol}[], Symbol[], Tuple{Symbol,Symbol}[], failure)
function inject(o, name, operation, moment)
    o.failure === nothing && return
    target, op, when, error = o.failure
    (name, operation, moment) == (target, op, when) && throw(error)
end
function begin_phase(o, name)
    inject(o, name, :begin, :before)
    push!(o.stack, name)
    push!(o.events, (name, :begin))
    inject(o, name, :begin, :after)
    return name
end
function end_phase(o, name, status)
    push!(o.attempts, (name, status))
    inject(o, name, :end, :before)
    isempty(o.stack) && error("unmatched end")
    last(o.stack) == name || error("wrong phase")
    pop!(o.stack)
    push!(o.events, (name, status))
    inject(o, name, :end, :after)
    return nothing
end
Base.@noinline throwing_body(error) = throw(error)
function caught(f)
    try
        return (; value = f(), error = nothing, backtrace = nothing)
    catch error
        return (; value = nothing, error, backtrace = catch_backtrace())
    end
end
has_body_frame(bt) = any(frame -> frame.func === :throwing_body, stacktrace(bt))
function pair(o; scattering_error = nothing, embedding_error = nothing)
    calls = Symbol[]
    scattering_value = (Ref(:scattering), Ref(:roundoff))
    embedding_value = (Ref(:embedding), Ref(:plus), Ref(:minus))
    result = caught() do
        _production_candidate_phases(
            token -> begin
                push!(calls, :scattering)
                scattering_error === nothing || throwing_body(scattering_error)
                scattering_value
            end,
            (token, value) -> begin
                @test value === scattering_value
                push!(calls, :embedding)
                embedding_error === nothing || throwing_body(embedding_error)
                embedding_value
            end,
            () -> begin_phase(o, :scattering),
            () -> begin_phase(o, :embedding),
            (token, status) -> end_phase(o, token, status),
        )
    end
    return (; result, calls, scattering_value, embedding_value)
end

struct ThrowingLogger <: AbstractLogger
    error::Sentinel
end
Logging.min_enabled_level(::ThrowingLogger) = Logging.Debug
Logging.shouldlog(::ThrowingLogger, args...) = true
Logging.catch_exceptions(::ThrowingLogger) = false
Logging.handle_message(logger::ThrowingLogger, args...; kwargs...) = throw(logger.error)

@testset "production phase lifecycle" begin
    @testset "healthy production pair preserves values and order" begin
        o = Observer()
        p = pair(o)
        @test p.result.error === nothing
        value = p.result.value
        @test value.scattering === p.scattering_value
        @test value.embedding === p.embedding_value
        @test p.calls == [:scattering, :embedding]
        @test o.events == [(:scattering, :begin), (:scattering, :completed),
                           (:embedding, :begin), (:embedding, :completed)]
        @test o.attempts == [(:scattering, :completed), (:embedding, :completed)]
        @test isempty(o.stack)
        @test value.started_ns > 0
        @test isfinite(value.scattering_seconds) && value.scattering_seconds >= 0
        @test isfinite(value.embedding_seconds) && value.embedding_seconds >= 0
    end
    @testset "SCI4 embedding admission preserves primary without replay" begin
        primary = Sentinel(:admission)
        o = Observer((:embedding, :begin, :before, primary))
        p = pair(o)
        @test p.result.error === primary
        @test p.calls == [:scattering]
        @test o.events == [(:scattering, :begin), (:scattering, :completed)]
        @test o.attempts == [(:scattering, :completed)]
        @test isempty(o.stack)
    end
    @testset "begin observer failure before and after external effect" begin
        for moment in (:before, :after)
            primary = Sentinel(:begin)
            o = Observer((:scattering, :begin, moment, primary))
            p = pair(o)
            @test p.result.error === primary
            @test isempty(p.calls)
            @test isempty(o.attempts)
            @test o.stack == (moment === :after ? [:scattering] : Symbol[])
        end
    end
    @testset "completed end observer failure is attempted once" begin
        for name in (:scattering, :embedding), moment in (:before, :after)
            primary = Sentinel(:end)
            o = Observer((name, :end, moment, primary))
            p = pair(o)
            @test p.result.error === primary
            @test p.calls == (name === :scattering ? [:scattering] : [:scattering, :embedding])
            expected = name === :scattering ? [(:scattering, :completed)] :
                       [(:scattering, :completed), (:embedding, :completed)]
            @test o.attempts == expected
            @test o.stack == (moment === :after ? Symbol[] : [name])
        end
    end
    @testset "body primary and original throw frame survive failed cleanup" begin
        primary = Sentinel(:body)
        secondary = Sentinel(:cleanup)
        o = Observer((:scattering, :end, :before, secondary))
        p = with_logger(NullLogger()) do
            pair(o; scattering_error = primary)
        end
        @test p.result.error === primary
        @test has_body_frame(p.result.backtrace)
        @test p.calls == [:scattering]
        @test o.attempts == [(:scattering, :failed)]
    end
    @testset "only body DomainError becomes a captured scientific failure" begin
        for name in (:scattering, :embedding), secondary in (nothing, Sentinel(:cleanup), DomainError(:cleanup))
            primary = DomainError(Ref(:solver), "fixture")
            o = secondary === nothing ? Observer() : Observer((name, :end, :after, secondary))
            p = with_logger(NullLogger()) do
                name === :scattering ? pair(o; scattering_error = primary) : pair(o; embedding_error = primary)
            end
            @test p.result.error === nothing
            failure = p.result.value
            @test failure isa _ProductionPhaseBodyFailure
            @test failure.error === primary
            @test has_body_frame(failure.backtrace)
            @test p.calls == (name === :scattering ? [:scattering] : [:scattering, :embedding])
            @test o.attempts == (name === :scattering ? [(:scattering, :failed)] :
                                [(:scattering, :completed), (:embedding, :failed)])
            @test isempty(o.stack)
        end
    end
    @testset "callback DomainError never enters scientific fallback" begin
        for name in (:scattering, :embedding), operation in (:begin, :end), moment in (:before, :after)
            primary = DomainError(Ref(:callback), "fixture")
            o = Observer((name, operation, moment, primary))
            p = pair(o)
            @test p.result.error === primary
            @test !(p.result.value isa _ProductionPhaseBodyFailure)
            @test count(row -> first(row) === name, o.attempts) == (operation === :end ? 1 : 0)
            operation === :begin && @test !(name in p.calls)
        end
    end
    @testset "throwing diagnostics cannot replace primary" begin
        for cleanup_fails in (false, true)
            primary = Sentinel(:body)
            o = cleanup_fails ? Observer((:scattering, :end, :before, Sentinel(:cleanup))) : Observer()
            result = with_logger(ThrowingLogger(Sentinel(:logger))) do
                caught() do
                    try
                        _production_phase_scope(
                            token -> throwing_body(primary),
                            () -> begin_phase(o, :scattering),
                            (token, status) -> end_phase(o, token, status),
                        )
                    catch error
                        bt = catch_backtrace()
                        # Mirrors the protected existing Dyson/candidate warning, even
                        # when failed-end cleanup itself completed successfully (F1).
                        _production_cleanup_preserving_primary(error, bt) do
                            @warn "existing failure warning" exception=(error, bt)
                        end
                        rethrow()
                    end
                end
            end
            @test result.error === primary
            @test has_body_frame(result.backtrace)
            @test o.attempts == [(:scattering, :failed)]
        end
    end
end

end # module
