using BendersX
using JuMP
using CPLEX
using Printf
using Random

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

include(normpath(joinpath(@__DIR__, "solver_defaults.jl")))
include(normpath(joinpath(@__DIR__, "script_utils.jl")))

options, _ = parse_script_args(ARGS)

instance = get_string_option(options, "instance", "ga250a-2")
seed = get_int_option(options, "seed", 1)
time_limit = get_float_option(options, "time_limit", 14400.0)
dcglp_time_limit = get_float_option(options, "dcglp_time_limit", 100.0)
dcglp_iter_limit = get_int_option(options, "dcglp_iter_limit", 250)
dcglp_halt_limit = get_int_option(options, "dcglp_halt_limit", 3)
frequency = get_int_option(options, "frequency", 500)
threads = get_int_option(options, "threads", 7)
reuse_dcglp = get_bool_option(options, "reuse_dcglp", false)
build_only = get_bool_option(options, "build_only", false)

Random.seed!(seed)

@info "PolarDCGLP hard UFL script" instance = instance seed = seed time_limit = time_limit frequency = frequency threads = threads reuse_dcglp = reuse_dcglp build_only = build_only

data = read_Simple_data(instance)

function customize_master_knapsack!(model::Model, data::UFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer,
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPINT" => 1e-9,
        "CPX_PARAM_EPRHS" => 1e-9,
        "CPX_PARAM_EPGAP" => 1e-6,
        MOI.Silent() => true,
    )
    set_optimizer(model, optimizer)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t[1:data.n_customers] >= -1e6)
    @constraint(model, sum(x) >= 2)
    @objective(model, Min, data.fixed_costs' * x + sum(t))
    return (x = x,), t
end

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
    verbose = true,
)

oracle_param = PolarDCGLPParam(
    dcglp_param;
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master = 2,
    fraction_of_benders_cuts_to_master = 0.05,
    reuse_dcglp = reuse_dcglp,

    zero_tol = 1e-9,
)

master = Master(data; customize = customize_master_knapsack!, optimizer = mip_optimizer)
set_optimizer_attribute(master.model, "CPX_PARAM_BRDIR", 1)

typical_oracles = [UFLKnapsackOracle(data), UFLKnapsackOracle(data)]
for oracle in typical_oracles
    set_parameter!(oracle, "add_only_violated_cuts", true)
end
disjunctive_oracle = PolarDCGLPOracle(master, typical_oracles, oracle_param)

lazy_oracle = UFLKnapsackOracle(data)
set_parameter!(lazy_oracle, "add_only_violated_cuts", true)
root_preprocessing = RootNodePreprocessing(
    lazy_oracle,
    BendersSeq,
    BendersSeqParam(
        time_limit = min(100.0, time_limit),
        gap_tolerance = 1e-9,
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
    @info "PolarDCGLP hard UFL script build completed without solve." instance = instance
else
    solve!(env)
    obj_value = try
        env.obj_value
    catch
        NaN
    end
    @info "PolarDCGLP hard UFL script finished" instance = instance termination_status = env.termination_status objective_value = obj_value
end
