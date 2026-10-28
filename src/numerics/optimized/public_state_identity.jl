# Canonical v1 encoding is in-memory only. Tag bytes and network byte order are normative.
# s/y: UInt64 byte length + UTF-8; i/f: exact 64 bits; b: 0/1; n: no payload;
# c/u: two recursively tagged values; q/a/d/t/r: lengths then ordered values.
_exact_u64(io,value)=write(io,hton(UInt64(value)))
function _exact_canonical!(io,value)
    if value===nothing
        write(io,UInt8('n'))
    elseif value isa Bool
        write(io,UInt8('b'),UInt8(value))
    elseif value isa Integer
        write(io,UInt8('i')); _exact_u64(io,reinterpret(UInt64,Int64(value)))
    elseif value isa Float64
        write(io,UInt8('f')); _exact_u64(io,reinterpret(UInt64,value))
    elseif value isa Float32
        write(io,UInt8('F')); write(io,hton(reinterpret(UInt32,value)))
    elseif value isa Complex
        write(io,UInt8('c')); _exact_canonical!(io,real(value)); _exact_canonical!(io,imag(value))
    elseif value isa Union{String,Symbol,VersionNumber}
        write(io,UInt8(value isa Symbol ? 'y' : 's'))
        bytes=codeunits(string(value)); _exact_u64(io,length(bytes)); write(io,bytes)
    elseif value isa Unitful.AbstractQuantity
        write(io,UInt8('u')); _exact_canonical!(io,string(unit(value)))
        _exact_canonical!(io,ustrip(value))
    elseif value isa AbstractArray
        write(io,UInt8('a')); _exact_canonical!(io,string(eltype(value)))
        _exact_u64(io,ndims(value))
        for dimension in size(value); _exact_u64(io,dimension); end
        for item in value; _exact_canonical!(io,item); end
    elseif value isa Tuple
        write(io,UInt8('q')); _exact_u64(io,length(value))
        for item in value; _exact_canonical!(io,item); end
    elseif value isa NamedTuple
        write(io,UInt8('r')); _exact_canonical!(io,keys(value))
        _exact_canonical!(io,Tuple(value))
    elseif value isa AbstractDict
        write(io,UInt8('d')); _exact_u64(io,length(value))
        ordered=sort!(collect(keys(value));by=x->bytes2hex(_exact_canonical_bytes(x)))
        for key in ordered; _exact_canonical!(io,key); _exact_canonical!(io,value[key]); end
    elseif value isa AbstractSet
        write(io,UInt8('t')); _exact_u64(io,length(value))
        for item in sort!(collect(value);by=x->bytes2hex(_exact_canonical_bytes(x)))
            _exact_canonical!(io,item)
        end
    elseif value isa Function || value isa Type || !isstructtype(typeof(value))
        throw(ArgumentError("unsupported exact identity value $(typeof(value))"))
    else
        write(io,UInt8('r')); _exact_canonical!(io,string(typeof(value)))
        names=fieldnames(typeof(value)); _exact_canonical!(io,names)
        for name in names; _exact_canonical!(io,getfield(value,name)); end
    end
    return io
end
function _exact_canonical_bytes(value)
    io=IOBuffer(); _exact_canonical!(io,value); return take!(io)
end
_exact_digest(value)=Tuple(sha256(_exact_canonical_bytes(value)))

function _validate_execution_runtime(execution::ExactExecutionIdentity)
    execution.julia==string(VERSION) || throw(ArgumentError("execution.julia differs from the loaded runtime"))
    execution.scalar_type=="Float64/ComplexF64" || throw(ArgumentError("execution.scalar_type differs"))
    execution.index_type==string(Int) || throw(ArgumentError("execution.index_type differs"))
    execution.reduction=="static_chunks_v1" || throw(ArgumentError("execution.reduction is unsupported"))
    execution.blas_threads==BLAS.get_num_threads() || throw(ArgumentError("execution.blas_threads differs"))
    execution.worker_count<=Base.Threads.nthreads(:default) || throw(ArgumentError(
        "execution.worker_count exceeds the actual Julia default thread pool; exact width would be clamped"))
    return nothing
end

function exact_restart_identity(problem::NEGFProblem,options::SolverOptions,
    production_options::ProductionOptions,execution::ExactExecutionIdentity;
    domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})
    _check_options(options); _check_production_options(production_options)
    _validate_execution_runtime(execution)
    domain_adaptation.mode !== :none && kernel_options===nothing &&
        throw(ArgumentError("identity.kernel_options is required for adaptive rebuild"))
    widths=(production_options.worker_count,production_options.parallel_backend,
        production_options.energy_chunk,production_options.hilbert_columns,production_options.residual_chunk)
    expected=(execution.worker_count,execution.parallel_backend,execution.energy_chunk,
        execution.hilbert_columns,execution.residual_chunk)
    widths==expected || throw(ArgumentError("identity.execution differs from numerical widths/chunks"))
    policy=(options,production_options.algorithms,production_options.physics_markers,
        production_options.verify_fft_roundoff,domain_adaptation,kernel_options)
    return ExactRestartIdentity(EXACT_STATE_VERSION,_exact_digest(problem),_exact_digest(policy),
        _exact_digest(execution))
end
