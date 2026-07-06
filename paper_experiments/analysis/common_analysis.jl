using CSV
using DataFrames
using Statistics

const ANALYSIS_ROOT = normpath(joinpath(@__DIR__, ".."))

function arg_value(name::AbstractString, default::AbstractString)
    prefix = string(name, "=")
    for (i, arg) in pairs(ARGS)
        startswith(arg, prefix) && return abspath(arg[length(prefix)+1:end])
        if arg == name && i < length(ARGS)
            return abspath(ARGS[i+1])
        end
    end
    return abspath(default)
end

const RAW_DIR = arg_value("--input_dir", joinpath(ANALYSIS_ROOT, "results", "raw"))
const PROCESSED_DIR = arg_value("--output_dir", joinpath(ANALYSIS_ROOT, "results", "processed"))

ensure_dir(path::AbstractString) = (isdir(path) || mkpath(path); path)

function raw_csv(name::AbstractString)
    path = joinpath(RAW_DIR, name)
    if isfile(path)
        return DataFrame(CSV.File(path))
    end

    paths = String[]
    for (root, _, files) in walkdir(RAW_DIR)
        for file in files
            file == name && push!(paths, joinpath(root, file))
        end
    end
    isempty(paths) && error("Missing raw result file: $(path), and no nested $(name) files found under $(RAW_DIR)")
    return vcat((DataFrame(CSV.File(p)) for p in sort(paths))...; cols = :union)
end

finite_values(xs) = [Float64(x) for x in skipmissing(xs) if isfinite(Float64(x))]
median_finite(xs) = isempty(finite_values(xs)) ? NaN : median(finite_values(xs))
maximum_finite(xs) = isempty(finite_values(xs)) ? NaN : maximum(finite_values(xs))
sum_bool(xs) = sum(Bool.(xs))

function write_processed(name::AbstractString, df::DataFrame)
    ensure_dir(PROCESSED_DIR)
    path = joinpath(PROCESSED_DIR, name)
    CSV.write(path, df)
    @info "wrote $(path)" rows = nrow(df)
    return path
end
