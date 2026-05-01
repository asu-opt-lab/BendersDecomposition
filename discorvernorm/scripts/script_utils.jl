function parse_script_args(args::Vector{String})
    options = Dict{String, String}()
    positionals = String[]
    for arg in args
        if startswith(arg, "--")
            keyval = split(arg[3:end], "="; limit = 2)
            if length(keyval) == 1
                options[keyval[1]] = "true"
            else
                options[keyval[1]] = keyval[2]
            end
        else
            push!(positionals, arg)
        end
    end
    return options, positionals
end

get_string_option(options::Dict{String, String}, key::String, default::String) = get(options, key, default)
get_int_option(options::Dict{String, String}, key::String, default::Int) = parse(Int, get(options, key, string(default)))
get_float_option(options::Dict{String, String}, key::String, default::Float64) = parse(Float64, get(options, key, string(default)))

function get_bool_option(options::Dict{String, String}, key::String, default::Bool)
    raw = lowercase(get(options, key, string(default)))
    if raw in ("true", "1", "yes", "y")
        return true
    elseif raw in ("false", "0", "no", "n")
        return false
    end
    throw(ArgumentError("Unable to parse boolean option --$(key)=$(raw)"))
end
