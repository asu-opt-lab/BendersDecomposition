function exception_message(e::Exception)
    if hasproperty(e, :msg)
        return getproperty(e, :msg)
    elseif e isa CompositeException
        return join(exception_message.(e.exceptions), "; ")
    elseif e isa TaskFailedException
        stack = Base.current_exceptions(e.task)
        return isempty(stack) ? sprint(showerror, e) : join((exception_message(item.exception) for item in stack), "; ")
    else
        return sprint(showerror, e)
    end
end

function exception_contains(e::Exception, ::Type{T}) where {T <: Exception}
    e isa T && return true
    e isa CompositeException && return any(child -> exception_contains(child, T), e.exceptions)
    e isa TaskFailedException && return any(item -> exception_contains(item.exception, T), Base.current_exceptions(e.task))
    return false
end

function last_lower_bound(log; default = Inf)
    return isempty(log.iterations) ? default : log.iterations[end].LB
end
