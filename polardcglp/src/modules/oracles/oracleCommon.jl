function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return
    registered = model[sym]
    if registered isa AbstractArray
        delete.(model, registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end

function format_sparse_terms(coeffs::SparseVector{Float64, Int}, var_name::String; zero_tol::Float64 = 1e-10)
    terms = String[]
    for p in sortperm(coeffs.nzind)
        idx = coeffs.nzind[p]
        val = coeffs.nzval[p]
        abs(val) <= zero_tol && continue
        sign = val >= 0 ? "+" : "-"
        push!(terms, @sprintf("%s %.6g %s[%d]", sign, abs(val), var_name, idx))
    end
    return terms
end

function format_hyperplane(h::BendersX.Hyperplane; zero_tol::Float64 = 1e-10)
    pieces = String[]
    if abs(h.a_0) > zero_tol || (nnz(h.a_x) == 0 && nnz(h.a_t) == 0)
        push!(pieces, @sprintf("%.6g", h.a_0))
    else
        push!(pieces, "0")
    end
    append!(pieces, format_sparse_terms(h.a_x, "x"; zero_tol = zero_tol))
    append!(pieces, format_sparse_terms(h.a_t, "t"; zero_tol = zero_tol))
    return join(pieces, " ") * " <= 0"
end

hyperplane_violation(h::BendersX.Hyperplane, x_value::Vector{Float64}, t_value::Vector{Float64}) =
    h.a_0 + dot(h.a_x, x_value) + dot(h.a_t, t_value)

function store_dcglp_disjunctive_cut!(
    oracle,
    cut::BendersX.Hyperplane,
    master_hyperplanes::Vector{BendersX.Hyperplane},
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    append_current_disjunctive_cut!(oracle, cut)
    include_disjunctive_cuts_to_hyperplanes && push!(master_hyperplanes, cut)
    return nothing
end

function BendersX.root_node_processing!(master::BendersX.AbstractMaster, preprocessing::BendersX.DisjunctiveRootNodePreprocessing)
    root_param = deepcopy(preprocessing.params)
    undo = JuMP.relax_integrality(master.model)

    start_time = time()
    try
        env_root_typical = preprocessing.seq_type(master, preprocessing.typical_oracle; param = root_param)
        BendersX.solve!(env_root_typical)

        root_param.time_limit -= time() - start_time

        env_root_disjunctive = preprocessing.seq_type(master, preprocessing.disjunctive_oracle; param = root_param)
        BendersX.solve!(env_root_disjunctive)
    finally
        undo()
    end

    return time() - start_time
end
