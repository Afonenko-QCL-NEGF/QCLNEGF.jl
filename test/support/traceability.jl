function rendered_public_doc(value)
    rendered = String[]
    # Public facade aliases retain documentation at the defining module.
    # Walk only package-owned modules, never unrelated imported dependencies.
    modules = Module[QCLNEGF]
    for mod in modules
        for name in names(mod; all = true, imported = false)
            isdefined(mod, name) || continue
            child = getfield(mod, name)
            child isa Module &&
                child !== mod &&
                parentmodule(child) === mod &&
                !(child in modules) &&
                push!(modules, child)
        end
    end
    for mod in modules, (binding, multidoc) in Base.Docs.meta(mod)
        binding isa Base.Docs.Binding || continue
        isdefined(binding.mod, binding.var) || continue
        getfield(binding.mod, binding.var) === value || continue
        push!(rendered, repr(multidoc))
    end
    isempty(rendered) && error("no documentation binding found for $(value)")
    return join(rendered, "\n")
end


function valid_public_doc_reference(reference, document_source)
    occursin("(@id $reference)", document_source) && return true
    startswith(reference, "QCLNEGF.") || return false
    value = QCLNEGF
    for part in split(chopprefix(reference, "QCLNEGF."), '.')
        value isa Module || return false
        name = Symbol(part)
        isdefined(value, name) || return false
        value = getfield(value, name)
    end
    return !isempty(rendered_public_doc(value))
end
