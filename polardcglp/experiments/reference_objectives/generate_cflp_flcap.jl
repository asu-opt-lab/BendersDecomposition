using BendersX
using CSV
using DataFrames
using JuMP
using CPLEX
using MathOptInterface
using Printf

const MOI = MathOptInterface

include(normpath(joinpath(@__DIR__, "..", "..", "src", "read_flcap.jl")))

const INSTANCE_NAMES = [
    "cap61", "cap62", "cap63", "cap64",
    "cap71", "cap72", "cap73", "cap74",
    "cap91", "cap92", "cap93", "cap94",
    "cap101", "cap102", "cap103", "cap104",
    "cap121", "cap122", "cap123", "cap124",
    "cap131", "cap132", "cap133", "cap134",
    "capa1", "capa2", "capa3", "capa4",
    "capb1", "capb2", "capb3", "capb4",
    "capc1", "capc2", "capc3", "capc4",
]

function build_full_cflp_mip(data::CFLPData; time_limit::Float64 = 3600.0)
    model = Model(optimizer_with_attributes(
        CPLEX.Optimizer,
        "CPXPARAM_Threads" => 7,
        "CPX_PARAM_EPGAP" => 0.0,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPOPT" => 1e-9,
        "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    ))
    set_time_limit_sec(model, time_limit)

    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min,
        sum(data.fixed_costs[i] * x[i] for i in 1:I) +
        sum(cost_demands .* y)
    )
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J], y[i, j] <= x[i])
    @constraint(model, capacity[i in 1:I],
        sum(data.demands[j] * y[i, j] for j in 1:J) <= data.capacities[i] * x[i]
    )

    return model
end

function load_existing_records(csv_path::AbstractString)
    isfile(csv_path) || return Set{String}()
    df = DataFrame(CSV.File(csv_path))
    return Set(String.(df.instance_name))
end

function ensure_csv_header(csv_path::AbstractString)
    if !isfile(csv_path)
        open(csv_path, "w") do io
            println(io, "instance_name,objective_value")
        end
    end
end

function append_record(csv_path::AbstractString, name::AbstractString, obj_value::Float64)
    open(csv_path, "a") do io
        @printf(io, "%s,%.10g\n", name, obj_value)
    end
end

function solve_and_record(name::AbstractString, csv_path::AbstractString;
                          time_limit::Float64 = 3600.0)
    @info "Solving CFLP FLCAP MIP" instance = name time_limit = time_limit
    data = read_flcap_data(name)
    model = build_full_cflp_mip(data; time_limit = time_limit)
    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL
        obj = objective_value(model)
        append_record(csv_path, name, obj)
        @info "Recorded optimal objective" instance = name objective = obj
        return :ok
    else
        @warn "Skipping non-optimal instance" instance = name status = status
        return :skipped
    end
end

function main(; time_limit::Float64 = 3600.0,
              csv_path::AbstractString = joinpath(@__DIR__, "cflp_flcap.csv"))
    ensure_csv_header(csv_path)
    existing = load_existing_records(csv_path)

    n_total = length(INSTANCE_NAMES)
    n_done = 0
    n_skipped = 0
    n_solved = 0

    for name in INSTANCE_NAMES
        if name in existing
            @info "Skipping already-recorded instance" instance = name
            n_done += 1
            continue
        end
        result = solve_and_record(name, csv_path; time_limit = time_limit)
        if result == :ok
            n_solved += 1
        else
            n_skipped += 1
        end
    end

    @info "Reference generation finished" total = n_total already_done = n_done newly_solved = n_solved skipped_non_optimal = n_skipped csv_path = csv_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
