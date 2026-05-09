using BendersX
using JuMP
using CPLEX
using Printf
using Random

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

include(normpath(joinpath(@__DIR__, "..", "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "..", "script_utils.jl")))

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "T100x100_5_3")
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 200.0)
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 1000.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 250)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 250)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

@info "SimplexNormDCGLP hard CFL script" instance = instance seed = seed time_limit = time_limit frequency = frequency reuse_dcglp = reuse_dcglp build_only = build_only

data = read_cfl_file(instance)

benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
)

dcglp_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
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

oracle_param = SimplexNormDCGLPParam(
    dcglp_param;
    split_index_selection_rule = MostFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master = true,
    fraction_of_benders_cuts_to_master = 0.5,
    reuse_dcglp = reuse_dcglp,

    zero_tol = 1e-9,
)

master = Master(data; customize = customize_master_model!, optimizer = mip_optimizer)

typical_oracles = [
    CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
    CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer),
]
disjunctive_oracle = SimplexNormDCGLPOracle(master, typical_oracles, oracle_param)

lazy_oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!, optimizer = optimizer)
root_preprocessing = RootNodePreprocessing(
    lazy_oracle,
    BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-6,
        stabilizing_x = ones(data.n_facilities),
        α = 0.9,
        λ = 0.1,
        verbose = true,
    ),
)
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
    @info "SimplexNormDCGLP hard CFL script build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "SimplexNormDCGLP hard CFL script finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
