function select_disjunctive_inequality(x_value::Vector{Float64}, split_selection_rule::SplitIndexSelectionRule; zero_tol = 1.0e-2)
    throw(UndefError("update select_disjunctive_inequality for $(typeof(split_selection_rule))"))
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::LargestFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : maximum(frac_indices)

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::MostFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : frac_indices[argmin(abs.(x_value[frac_indices] .- 0.5))]

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::RandomFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : rand(frac_indices)

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function get_split_index(oracle::SplitOracle)
    oracle.param.split_index_selection_rule isa SimpleSplit ||
        throw(AlgorithmException("get_split_index is only valid for simple split rules."))
    isempty(oracle.splits) &&
        throw(AlgorithmException("get_split_index requires at least one selected split."))
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::SplitOracle)
    dcglp = oracle.dcglp
    phi, phi_0 = oracle.splits[end]

    delete_registered_constraints!(dcglp, :con_split_kappa)
    delete_registered_constraints!(dcglp, :con_split_nu)

    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])
end

function retrieve_zero_one(x_value::Vector{Float64}, zero_tol)
    zero_indices = findall(x -> isapprox(x, 0.0; atol = zero_tol), x_value)
    one_indices = findall(x -> isapprox(x, 1.0; atol = zero_tol), x_value)
    return zero_indices, one_indices
end

function add_lifting_constraints!(dcglp::Model, zero_indices::Vector{Int}, one_indices::Vector{Int})
    delete_registered_constraints!(dcglp, :con_zeta)
    delete_registered_constraints!(dcglp, :con_xi)

    !isempty(zero_indices) &&
        @constraint(dcglp, con_zeta[i in 1:2, j in eachindex(zero_indices)], 0 >= dcglp[:omega_x][i, zero_indices[j]])
    !isempty(one_indices) &&
        @constraint(dcglp, con_xi[i in 1:2, j in eachindex(one_indices)], 0 >= dcglp[:omega_0][i] - dcglp[:omega_x][i, one_indices[j]])
end

function choose_split_and_update_lifting!(
    oracle::SplitOracle,
    x_value::Vector{Float64},
)
    push!(
        oracle.splits,
        select_disjunctive_inequality(
            x_value,
            oracle.param.split_index_selection_rule;
            zero_tol = oracle.param.zero_tol,
        ),
    )
    replace_disjunctive_inequality!(oracle)

    zero_indices, one_indices = oracle.param.lift ? retrieve_zero_one(x_value, oracle.param.zero_tol) : (Int[], Int[])
    add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices)

    return zero_indices, one_indices
end