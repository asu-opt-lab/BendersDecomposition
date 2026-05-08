import BendersX: CFLPData

"""
    read_flcap_data(filename; filepath) -> CFLPData

Read an FLCAP-format CFLP instance.

File layout:
- line 1: `n_facilities n_customers`
- next `n_facilities` lines: `capacity fixed_cost`
- next line: `n_customers` demand values
- next `n_facilities` lines: each line lists total assignment costs from that
  facility to all customers

The original Beasley/OR-Library files store the cost of assigning all demand of
customer `j` to facility `i`. Internally `CFLPData.costs[i, j]` is interpreted
as a unit cost, so we divide each stored total cost by `demands[j]` here.
"""
function read_flcap_data(filename::AbstractString;
                         filepath = joinpath(POLARDCGLP_ROOT, "data", "FLCAP"))
    fullpath = joinpath(filepath, filename)
    f = open(fullpath)

    vals = split(readline(f))
    n_facilities = parse(Int, vals[1])
    n_customers = parse(Int, vals[2])

    capacities = zeros(Float64, n_facilities)
    fixed_costs = zeros(Float64, n_facilities)
    for i in 1:n_facilities
        vals = split(readline(f))
        capacities[i] = parse(Float64, vals[1])
        fixed_costs[i] = parse(Float64, vals[2])
    end

    demands = zeros(Float64, n_customers)
    vals = split(readline(f))
    for j in 1:n_customers
        demands[j] = parse(Float64, vals[j])
    end

    costs = zeros(Float64, n_facilities, n_customers)
    for i in 1:n_facilities
        vals = split(readline(f))
        for j in 1:n_customers
            costs[i, j] = parse(Float64, vals[j]) / demands[j]
        end
    end
    close(f)

    return CFLPData(n_facilities, n_customers, capacities, demands, fixed_costs, costs)
end
