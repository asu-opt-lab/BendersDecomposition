using BendersX
using CSV
using DataFrames
using JuMP

const OUTPUT_DIR = @__DIR__

function solve_reference_objective(data, instance_name::AbstractString)
    mip_model = Model()
    customize_mip_model!(mip_model, data)
    optimize!(mip_model)
    termination_status(mip_model) == OPTIMAL || error("Reference MILP for $(instance_name) did not terminate optimally: $(termination_status(mip_model))")
    return objective_value(mip_model)
end

function write_reference_csv(filename::AbstractString, rows)
    output_path = joinpath(OUTPUT_DIR, filename)
    df = DataFrame(rows)
    temp_path, temp_io = mktemp(dirname(output_path))
    close(temp_io)
    try
        CSV.write(temp_path, df)
        mv(temp_path, output_path; force = true)
    finally
        isfile(temp_path) && rm(temp_path; force = true)
    end
    @info "Wrote $(nrow(df)) reference objectives to $(output_path)"
end

function build_uflp_rows()
    rows = NamedTuple[]
    for i in setdiff(1:71, [67])
        instance_name = "p$i"
        data = read_uflp_benchmark_data(instance_name)
        push!(rows, (instance_name = instance_name, objective_value = solve_reference_objective(data, instance_name)))
    end
    return rows
end

function build_cflp_rows()
    rows = NamedTuple[]
    for i in setdiff(1:71, [67])
        instance_name = "p$i"
        data = read_cflp_benchmark_data(instance_name)
        push!(rows, (instance_name = instance_name, objective_value = solve_reference_objective(data, instance_name)))
    end
    return rows
end

function build_scflp_rows()
    rows = NamedTuple[]
    for i in 1:5
        instance_name = "f25-c50-s64-r10-$i"
        data = read_stochastic_capacited_facility_location_problem(instance_name)
        push!(rows, (instance_name = instance_name, objective_value = solve_reference_objective(data, instance_name)))
    end
    return rows
end

function build_snip_rows()
    rows = NamedTuple[]
    instance_name = "instance=0;snipno=0;budget=30.0"
    data = read_snip_data(0, 0, 30.0)
    push!(rows, (instance_name = instance_name, objective_value = solve_reference_objective(data, instance_name)))
    return rows
end

function main()
    mkpath(OUTPUT_DIR)
    uflp_rows = build_uflp_rows()
    cflp_rows = build_cflp_rows()
    scflp_rows = build_scflp_rows()
    snip_rows = build_snip_rows()

    write_reference_csv("uflp.csv", uflp_rows)
    write_reference_csv("cflp.csv", cflp_rows)
    write_reference_csv("scflp.csv", scflp_rows)
    write_reference_csv("snip.csv", snip_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
