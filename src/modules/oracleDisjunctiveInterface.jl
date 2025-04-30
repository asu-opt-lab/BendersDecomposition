

function select_disjunctive_inequality(x_value::Vector{Float64}, ::LargestFractional; zero_tol = 1e-9)
    
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : maximum(frac_indices)

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Largest fractional simple split index: $index, phi: $phi, phi_0: $phi_0"
    
    return phi, phi_0
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::MostFractional; zero_tol = 1e-9)

    gap_x = @. abs(x_value - 0.5)

    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : argmin(gap_x)
    
    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Most fractional simple split index: $index, phi: $phi, phi_0: $phi_0"

    return phi, phi_0
end
function select_disjunctive_inequality(x_value::Vector{Float64}, ::RandomFractional; zero_tol = 1e-9)
    
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : rand(frac_indices)
    
    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Random simple split index: $index, phi: $phi, phi_0: $phi_0"
    
    return phi, phi_0
end

function add_disjunctive_cuts!(oracle::AbstractDisjunctiveOracle, ::NoDisjunctiveCuts)
    # do nothing
end
function add_disjunctive_cuts!(oracle::AbstractDisjunctiveOracle, ::AllDisjunctiveCuts)
    # do nothing; added at the time of generation
end
function add_disjunctive_cuts!(oracle::AbstractDisjunctiveOracle, ::DisjunctiveCutsSmallerIndices)
    
    @assert typeof(oracle.oracle_param.split_index_selection_rule) <: SimpleSplit

    dcglp = oracle.dcglp
    # remove all disjunctive cuts from DCGLP
    if haskey(dcglp, :con_disjunctive)
        delete.(dcglp, dcglp[:con_disjunctive]) 
        unregister(dcglp, :con_disjunctive)
    end
    
    # get variable index used for current split
    index = get_split_index(oracle)

    disjunctiveCuts = index > 1 ? reduce(vcat, [oracle.disjunctiveCutsByIndex[i] for i = 1:index-1]) : Vector{Hyperplane}()
    cuts = Vector{AffExpr}()
    for k = 1:2 # add to both kappa and nu systems
        append!(cuts, hyperplanes_to_expression(dcglp, disjunctiveCuts, dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
    end

    if typeof(oracle.oracle_param.norm) <: LpNorm
        @constraint(dcglp, con_disjunctive, 0 .>= cuts)
    elseif typeof(oracle.oracle_param.norm) <: StandardNorm 
        @constraint(dcglp, con_disjunctive, dcglp[:tau] .>= cuts)
    else
        throw(UndefError("update `add_disjunctive_cuts!` for $(typeof(oracle.oracle_param.norm))"))
    end
end

