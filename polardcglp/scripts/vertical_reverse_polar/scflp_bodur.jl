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
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 1000.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 250)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 250)
threads = get_int_option(options, "threads", 7)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
strengthened = get_bool_option(options, "strengthened", true)
lift = get_bool_option(options, "lift", false)
oracle_name = get_string_option(options, "oracle", "knapsack")
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

const SCFLP_DIR = normpath(joinpath(@__DIR__, "..", "..", "data", "SCFLP_bodur"))

@info "vertical_reverse_polar SCFLP Bodur script" instance = instance dataset_dir = SCFLP_DIR seed = seed time_limit = time_limit frequency = frequency threads = threads reuse_dcglp = reuse_dcglp strengthened = strengthened lift = lift oracle = oracle_name build_only = build_only

data = read_scflp_bodur(instance; filepath = SCFLP_DIR)

benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
)

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

dcglp_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => threads,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    "CPX_PARAM_EPOPT" => 1e-9,
    MOI.Silent() => true,
)

dcglp_param = DcglpParam(
    dcglp_optimizer;
    time_limit = dcglp_time_limit,
    gap_tolerance = 1e-3,
    halt_limit = dcglp_halt_limit,
    iter_limit = dcglp_iter_limit,
    verbose = true,
)

oracle_param = VerticalReversePolarDCGLPParam(
    dcglp_param;
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master = 2,
    fraction_of_benders_cuts_to_master = 0.05,
    reuse_dcglp = reuse_dcglp,
    strengthened = strengthened,
    lift = lift,
    zero_tol = 1e-9,
)

master = Master(data; customize = customize_master_model!, optimizer = master_optimizer_local)

typical_oracles = [
    build_scflp_bodur_oracle(data, master, oracle_name; optimizer = sub_optimizer_local),
    build_scflp_bodur_oracle(data, master, oracle_name; optimizer = sub_optimizer_local),
]
disjunctive_oracle = VerticalReversePolarDCGLPOracle(master, typical_oracles, oracle_param)

lazy_oracle = build_scflp_bodur_oracle(data, master, oracle_name; optimizer = sub_optimizer_local)
# root_preprocessing = RootNodePreprocessing(
#     lazy_oracle,
#     BendersSeqInOut,
#     BendersSeqInOutParam(
#         time_limit = min(300.0, time_limit),
#         gap_tolerance = 1e-6,
#         stabilizing_x = ones(data.n_facilities),
#         α = 0.9,
#         λ = 0.1,
#         verbose = true,
#     ),
# )

root_preprocessing = NoRootNodePreprocessing()
lazy_callback = LazyCallback(lazy_oracle)
user_callback = UserCallback(disjunctive_oracle; params = UserCallbackParam(frequency = frequency))

env = BendersBnB(
    master,
    root_preprocessing,
    lazy_callback,
    user_callback;
    param = benders_param,
)

if build_only
    @info "vertical_reverse_polar SCFLP Bodur script build completed without solve." instance = instance oracle = oracle_name
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "vertical_reverse_polar SCFLP Bodur script finished" instance = instance oracle = oracle_name termination_status = env.termination_status objective_value = obj_value
end
