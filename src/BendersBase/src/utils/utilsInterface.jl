export copy_variables!

function copy_variables!(model::Model, x::NamedTuple)
    for (key, value) in pairs(x)
        if value isa JuMP.VariableRef
            # needs test
            expr = :( @variable(model, $(key), base_name = $(string(key))) )
            eval(expr)
        elseif value isa JuMP.Containers.Array || value isa JuMP.Containers.DenseAxisArray
            nd = length(axes(value))
            axes_tuple = ntuple(d -> axes(value, d), nd)
            ranges = [axes_tuple...]

            expr = :( @variable($(model), $(key)[$(ranges...)], base_name = $(string(key))) )
            eval(expr)
        elseif value isa JuMP.Containers.SparseAxisArray    
            expr = :( @variable($(model), $(Symbol(key,"temp"))[$(eachindex(value))], base_name = $(string(key))) )
            eval(expr)
            
            varmap = Dict(index => model[Symbol(key,"temp")][index] for index in eachindex(value))
            var = Containers.SparseAxisArray(varmap)

            model[key] = var
            unregister(model, Symbol(key,"temp"))
        else
            @error "We are not currently accepting a custom JuMP variable container: $(typeof(value))"
        end
    end

    x_copy = (; (key => model[key] for (key, _) in pairs(x))...)

    return x_copy
end