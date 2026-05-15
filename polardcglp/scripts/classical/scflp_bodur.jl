using BendersX
using JuMP
using CPLEX
using Printf
using Random

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))

function build_scflp_bodur_oracle(data::SCFLPBodurData, master::Master, oracle_name::String; optimizer)
    normalized = lowercase(strip(oracle_name))
    if normalized in ("classical", "cls")
        return SeparableOracle(
            data,
            master,
            ClassicalOracle(),
            data.n_scenarios;
            customize = customize_sub_model!,
            optimizer = optimizer,
        )
    elseif normalized in ("knapsack", "cfl_knapsack", "cflp_knapsack")
        return SeparableOracle(
            data,
            master,
            CFLKnapsackOracle(),
            data.n_scenarios;
            customize = customize_sub_model!,
            optimizer = optimizer,
        )
    end

    throw(ArgumentError("Unsupported SCFLP Bodur oracle `$(oracle_name)`. Use `classical` or `knapsack`."))
end

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "cap102-s250")
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 14400.0)
threads = get_int_option(options, "threads", 7)
oracle_name = get_string_option(options, "oracle", "knapsack")
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "SCFLP_bodur"))

@info "classical SCFLP Bodur script" instance = instance dataset_dir = SCFLP_DIR seed = seed time_limit = time_limit threads = threads oracle = oracle_name build_only = build_only

data = read_scflp_bodur(instance; filepath = SCFLP_DIR)

master_optimizer_local = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => threads,
    "CPX_PARAM_EPINT" => 1e-9,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    MOI.Silent() => true,
)

sub_optimizer_local = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => threads,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    MOI.Silent() => true,
)

master = Master(data; customize = customize_master_model!, optimizer = master_optimizer_local)

typical_oracle = build_scflp_bodur_oracle(data, master, oracle_name; optimizer = sub_optimizer_local)

root_preprocessing = RootNodePreprocessing(
    typical_oracle,
    BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = min(500.0, time_limit),
        gap_tolerance = 1e-9,
        stabilizing_x = ones(data.n_facilities),
        α = 0.9,
        λ = 0.1,
        verbose = true,
    ),
)
# root_preprocessing = NoRootNodePreprocessing()

lazy_callback = LazyCallback(typical_oracle)
user_callback = NoUserCallback()

benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
)

env = BendersBnB(
    master,
    root_preprocessing,
    lazy_callback,
    user_callback;
    param = benders_param,
)

if build_only
    @info "classical SCFLP Bodur script build completed without solve." instance = instance oracle = oracle_name
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "classical SCFLP Bodur script finished" instance = instance oracle = oracle_name termination_status = env.termination_status objective_value = obj_value
end
