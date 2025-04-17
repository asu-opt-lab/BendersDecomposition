function add_normalization_constraint(dcglp::Model, norm::LpNorm)
    # CPLEX only accepts p=1,2,Inf
    # if conic solver, we can use the following line
    # @constraint(dcglp, concone, var_vec in MOI.NormCone(norm.p, data.dim_x + data.dim_t + 1))
    var_vec = [dcglp[:tau]; dcglp[:sx]; dcglp[:st]]
    
    if norm.p == 1.0
        @constraint(dcglp, concone, var_vec in MOI.NormOneCone(length(var_vec)))
    elseif norm.p == 2.0
        @constraint(dcglp, concone, var_vec in MOI.SecondOrderCone(length(var_vec)))
    elseif norm.p == Inf
        @constraint(dcglp, concone, var_vec in MOI.NormInfinityCone(length(var_vec)))
    else
        throw(UndefError("Unsupported LpNorm: p=$(norm.p)"))
    end
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::LargestFractional; zero_tol = 1e-2)
    
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : argmax(frac_indices)
    
    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Largest fractional simple split index: $index"
    
    return phi, phi_0
end
function select_disjunctive_inequality(x_value::Vector{Float64}, ::MostFractional; zero_tol = 1e-2)

    gap_x = @. abs(x_value - 0.5)

    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : argmin(gap_x)
    
    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Most fractional simple split index: $index"

    return phi, phi_0
end
function select_disjunctive_inequality(x_value::Vector{Float64}, ::RandomFractional; zero_tol = 1e-2)
    
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(1:length(x_value))) : rand(frac_indices)
    
    phi = spzeros(length(x_value))
    phi[index] = 1.0
    phi_0 = 0.0

    @debug "Random simple split index: $index"
    
    return phi, phi_0
end

function add_disjunctive_cuts!(oracle::DisjunctiveOracle, ::NoDisjunctiveCuts)
    # do nothing
end
function add_disjunctive_cuts!(oracle::DisjunctiveOracle, ::AllDisjunctiveCuts)
    # do nothing; added at the time of generation
end
function add_disjunctive_cuts!(oracle::DisjunctiveOracle, ::DisjunctiveCutsSmallerIndices)
    
    @assert typeof(oracle.oracle_param.split_index_selection_rule) <: SimpleSplit

    dcglp = oracle.dcglp
    # remove all disjunctive cuts from DCGLP
    if haskey(dcglp, :con_disjunctive)
        delete(dcglp, dcglp[:con_disjunctive]) 
        unregister(dcglp, :con_disjunctive)
    end
    
    # get variable index used for current split
    index = get_split_index(oracle)

    disjunctiveCuts = index > 1 ? reduce(vcat, [oracle.disjunctiveCutsByIndex[i] for i = 1:index-1]) : Vector{Hyperplane}()
    cuts = Vector{AffExpr}()
    for k = 1:2 # add to both kappa and nu systems
        append!(cuts, hyperplanes_to_expression(dcglp, disjunctiveCuts, dcglp[:omega_x][k,:], dcglp[:omega_t][k,:], dcglp[:omega_0][k]))
    end

    @constraint(dcglp, con_disjunctive, 0 .>= cuts)
end

