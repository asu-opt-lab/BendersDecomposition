using Random
using JSON

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP


const FLCAP_DIR = normpath(joinpath(@__DIR__, "..", "data", "FLCAP"))
const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "data", "SCFLP"))
const SCENARIO_SIZES = [256, 512, 1024]

mkpath(OUTPUT_DIR)

flcap_files = sort(filter(f -> startswith(f, "cap"), readdir(FLCAP_DIR)))

@assert length(flcap_files) == 36 "Expected 36 FLCAP files, found $(length(flcap_files))"

counter = 0
for capname in flcap_files
    for S in SCENARIO_SIZES
        payload = generate_scflp_from_flcap(capname, S; flcap_dir=FLCAP_DIR)
        out_path = joinpath(OUTPUT_DIR, string(capname, "-s", S, ".json"))
        write_scflp_json(payload, out_path)
        global counter += 1
        @info "Wrote $(counter)/$(length(flcap_files) * length(SCENARIO_SIZES)): $(basename(out_path))"
    end
end

@info "Total files written: $counter"
@assert counter == length(flcap_files) * length(SCENARIO_SIZES) "Expected $(length(flcap_files) * length(SCENARIO_SIZES)) files, wrote $counter"
