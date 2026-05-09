using Random
using JSON

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP


const FLCAP_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "FLCAP"))
const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "SCFLP_bodur"))
const SCENARIO_SIZES = [250, 500, 1500]

mkpath(OUTPUT_DIR)

flcap_files = sort(filter(f -> startswith(f, "cap"), readdir(FLCAP_DIR)))

@assert length(flcap_files) == 36 "Expected 36 FLCAP files, found $(length(flcap_files))"

counter = 0
for capname in flcap_files
    for S in SCENARIO_SIZES
        payload = generate_scflp_bodur(capname, S; flcap_dir=FLCAP_DIR)
        out_path = joinpath(OUTPUT_DIR, string(capname, "-s", S, ".json"))
        write_scflp_json(payload, out_path)
        global counter += 1
        @info "Wrote $(counter)/$(length(flcap_files) * length(SCENARIO_SIZES)): $(basename(out_path))"
    end
end

@info "Total files written: $counter"
@assert counter == length(flcap_files) * length(SCENARIO_SIZES) "Expected $(length(flcap_files) * length(SCENARIO_SIZES)) files, wrote $counter"
