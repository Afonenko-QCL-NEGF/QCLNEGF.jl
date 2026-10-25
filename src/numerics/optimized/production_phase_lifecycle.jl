# Private phase ownership. A returned begin token receives at most one end
# attempt; callbacks that throw after an external effect are never replayed.
struct _ProductionPhaseBodyFailure{B}
    error::DomainError
    backtrace::B
end

function _production_cleanup_preserving_primary(
    cleanup::F,
    primary,
    primary_backtrace,
) where {F}
    try
        cleanup()
    catch secondary
        secondary_backtrace = catch_backtrace()
        try
            @warn "production cleanup failed while preserving primary exception" exception=(
                secondary,
                secondary_backtrace,
            ) primary_exception=(primary, primary_backtrace)
        catch
            # A diagnostic observer must not replace an already known primary.
        end
    end
    return nothing
end

function _production_phase_scope(
    body::F,
    begin_phase::B,
    end_phase::E;
    capture_domain_error::Bool = false,
) where {F,B,E}
    token = begin_phase()
    started_ns = time_ns()
    value = try
        body(token)
    catch error
        bt = catch_backtrace()
        _production_cleanup_preserving_primary(error, bt) do
            end_phase(token, :failed)
        end
        if capture_domain_error && error isa DomainError
            return _ProductionPhaseBodyFailure(error, bt)
        end
        rethrow()
    end
    # End callback errors are operational primary failures, outside the body
    # catch, even when their exception type is DomainError.
    end_phase(token, :completed)
    return (; value, started_ns, seconds = (time_ns()-started_ns)*1e-9)
end

function _production_candidate_phases(
    scattering_body::S,
    embedding_body::E,
    begin_scattering::BS,
    begin_embedding::BE,
    end_phase::END,
) where {S,E,BS,BE,END}
    scattering = _production_phase_scope(
        scattering_body,
        begin_scattering,
        end_phase;
        capture_domain_error = true,
    )
    scattering isa _ProductionPhaseBodyFailure && return scattering
    embedding = _production_phase_scope(
        token -> embedding_body(token, scattering.value),
        begin_embedding,
        end_phase;
        capture_domain_error = true,
    )
    embedding isa _ProductionPhaseBodyFailure && return embedding
    return (;
        scattering = scattering.value,
        embedding = embedding.value,
        started_ns = scattering.started_ns,
        scattering_seconds = scattering.seconds,
        embedding_seconds = embedding.seconds,
    )
end
