export Mip

mutable struct Mip <: AbstractMip
    model::Model

    function Mip(data::Data)

        @debug "Building Mip Problem"
    
        model = Model()

        @variable(model, x[1:data.dim_x], Bin)
        
        new(model)
    end
end
