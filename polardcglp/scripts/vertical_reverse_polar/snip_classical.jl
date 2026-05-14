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

instance_no = get_int_option(options, "instance_no", 0)
snip_no = get_int_option(options, "snip_no", 3)
budget = get_float_option(options, "budget", 30.0)
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 14400.0)
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 100.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 250)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 500)
threads = get_int_option(options, "threads", 7)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
strengthened = get_bool_option(options, "strengthened", true)
lift = get_bool_option(options, "lift", true)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

sub_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 1,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    MOI.Silent() => true,
)

@info "vertical_reverse_polar SNIP script" instance_no = instance_no snip_no = snip_no budget = budget seed = seed time_limit = time_limit frequency = frequency threads = threads reuse_dcglp = reuse_dcglp strengthened = strengthened lift = lift build_only = build_only

data = read_snip_data(instance_no, snip_no, budget)

benders_param = BendersBnBParam(
    time_limit = time_limit,
    gap_tolerance = 1e-6,
    verbose = true,
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
    verbose = false,
)

oracle_param = VerticalReversePolarDCGLPParam(
    dcglp_param;
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master = 1,
    fraction_of_benders_cuts_to_master = 0.05,
    reuse_dcglp = reuse_dcglp,
    strengthened = strengthened,
    lift = lift,
    zero_tol = 1e-9,
)

master = Master(data; customize = customize_master_model!, optimizer = master_optimizer)

typical_oracles = [
    SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!, optimizer = sub_optimizer),
    SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!, optimizer = sub_optimizer),
]
disjunctive_oracle = VerticalReversePolarDCGLPOracle(master, typical_oracles, oracle_param)

lazy_oracle = SeparableOracle(data, master, ClassicalOracle(), data.num_scenarios; customize = customize_sub_model!, optimizer = sub_optimizer)
root_preprocessing = RootNodePreprocessing(
    lazy_oracle,
    BendersSeq,
    BendersSeqParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-6,
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
    @info "vertical_reverse_polar SNIP script build completed without solve." instance_no = instance_no snip_no = snip_no budget = budget
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "vertical_reverse_polar SNIP script finished" instance_no = instance_no snip_no = snip_no budget = budget termination_status = env.termination_status objective_value = obj_value
end
