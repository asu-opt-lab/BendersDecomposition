# tight tolerance for dcglp; adjust tol for generate_cut with omega_0
# dcglp is not solved -> cannot proceed the algorithm, and should be terminated. Added an optimal argument to dcglp `throw_typical_cuts_for_errors=false`
# remove other_constraints from classical oracle.
# add a method for generating optimal vertex for SpecializedBendersSeq

using BendersX
using CSV
using DataFrames
using Test
using JuMP
using CPLEX
using Logging
include(joinpath(@__DIR__, "solver_defaults.jl"))
global_logger(ConsoleLogger(stderr, Logging.Debug))

# loop parameters
specialized_benders_param = SpecializedBendersSeqParam(;
                                                        time_limit = 1000.0,
                                                        lp_gap_tolerance = 1e-9,
                                                        integrality_tolerance = 1e-9,
                                                        verbose = true
                                                        )
dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_LPMETHOD" => 3, MOI.Silent() => true)
dcglp_param = DcglpParam(dcglp_optimizer;
                        time_limit = 1000.0,
                        gap_tolerance = 1e-9,
                        halt_limit = Int(1e9),
                        iter_limit = Int(1e9),
                        verbose = true
                        )

@testset verbose = true "UFLP Specialized Sequential Benders Tests" begin
    ufl_reference_path = normpath(joinpath(@__DIR__, "reference_objectives", "uflp.csv"))
    ufl_reference_df = DataFrame(CSV.File(ufl_reference_path))
    @assert nrow(ufl_reference_df) == length(unique(ufl_reference_df.instance_name)) "Duplicate UFLP reference objectives found in $(ufl_reference_path)"
    ufl_reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(ufl_reference_df))
    # instances = setdiff(1:71, [67])
    instances = [49] # only instance that uses SpecializedBendersSeq
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            data = read_uflp_benchmark_data("p$(i)")

            instance_name = "p$i"
            @assert haskey(ufl_reference_objectives, instance_name) "Missing UFLP reference objective for $(instance_name) in $(ufl_reference_path)"
            mip_opt_val = ufl_reference_objectives[instance_name]

            @testset "Classical oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [Inf] # too slow

                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]
                # for strengthened in [false], add_benders_cuts_to_master in [true], reuse_dcglp in [true], p in [1.0]

                # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] # fail
                # for strengthened in [false], add_benders_cuts_to_master in [false], reuse_dcglp in [true], p in [1.0] #fail
                    @info "solving p$i - begin oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = SplitOracleParam(LpDistanceNormalization(LpNorm(p)); dcglp_param = dcglp_param,
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
                                                                strengthened=strengthened,
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master,
                                                                fraction_of_benders_cuts_to_master = 0.5,
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                        set_optimizer_attribute(master.model, "CPX_PARAM_LPMETHOD", 1)
                        typical_oracles = [ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer); ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)] # for kappa & nu
                        disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles), oracle_param)
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)

                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end

            @testset "fat knapsack oracle" begin
                function update_master_model!(model::Model, data::UFLPData)
                    optimizer = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)
                    set_optimizer(model, optimizer)
                    @variable(model, x[1:data.n_facilities], Bin)
                    @variable(model, t[1:data.n_customers] >= -1e6)
                    @constraint(model, sum(x) >= 2)
                    @objective(model, Min, data.fixed_costs'* x + sum(t))
                    return (x = x, ), t
                end

                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                    @info "solving p$i - fat Knapsack oracle - Specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = SplitOracleParam(LpDistanceNormalization(LpNorm(p)); dcglp_param = dcglp_param,
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
                                                                strengthened=strengthened,
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master,
                                                                fraction_of_benders_cuts_to_master = 0.5,
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                        set_optimizer_attribute(master.model, "CPX_PARAM_LPMETHOD", 1)
                        typical_oracles = [UFLKnapsackOracle(data); UFLKnapsackOracle(data)] # for kappa & nu
                        disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles), oracle_param)
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)

                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end
        end
    end
end

@testset verbose = true "CFLP Specialized Sequential Benders Tests" begin
    cfl_reference_path = normpath(joinpath(@__DIR__, "reference_objectives", "cflp.csv"))
    cfl_reference_df = DataFrame(CSV.File(cfl_reference_path))
    @assert nrow(cfl_reference_df) == length(unique(cfl_reference_df.instance_name)) "Duplicate CFLP reference objectives found in $(cfl_reference_path)"
    cfl_reference_objectives = Dict(String(row.instance_name) => Float64(row.objective_value) for row in eachrow(cfl_reference_df))
    # instances = setdiff(1:71, [67])
    instances = [25 32 34 36 49 51]
    # numerical issue: 29 30 31 35
    # success: 25 32 34 36 49 51
    # lp: 1:24, 26:28, 33, 37-48, 50, 52-66, 68-71
    for i in instances
        @testset "Instance: p$i" begin
            # Load problem data if necessary
            data = read_cflp_benchmark_data("p$i")
            # data = read_GK_data("f100-c100-r3-1")

            instance_name = "p$i"
            @assert haskey(cfl_reference_objectives, instance_name) "Missing CFLP reference objective for $(instance_name) in $(cfl_reference_path)"
            mip_opt_val = cfl_reference_objectives[instance_name]

            @testset "Classical oracle" begin
                    # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                    @info "solving CFLP p$i - disjunctive oracle/classical - specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                    @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                        oracle_param = SplitOracleParam(LpDistanceNormalization(LpNorm(p)); dcglp_param = dcglp_param,
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
                                                                strengthened=strengthened,
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master,
                                                                fraction_of_benders_cuts_to_master = 0.5,
                                                                reuse_dcglp=reuse_dcglp)

                        master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                        typical_oracles = [ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer); ClassicalOracle(data, master; model = update_sub_model!, optimizer = optimizer)] # for kappa & nu
                        disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles), oracle_param)
                        env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)

                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end

            @testset "Knapsack oracle" begin
                # for strengthened in [true; false], add_benders_cuts_to_master in [true; false], reuse_dcglp in [true; false], p in [1.0; Inf]
                for strengthened in [true], add_benders_cuts_to_master in [true], reuse_dcglp in [false], p in [1.0]
                    @info "solving CFLP p$i - disjunctive oracle/knapsack- specialized seq - strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p"
                        @testset "strgthnd $strengthened; benders2master $add_benders_cuts_to_master reuse $reuse_dcglp p $p" begin
                            oracle_param = SplitOracleParam(LpDistanceNormalization(LpNorm(p)); dcglp_param = dcglp_param,
                                                                split_index_selection_rule = LargestFractional(),
                                                                disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
                                                                strengthened=strengthened,
                                                                add_benders_cuts_to_master=add_benders_cuts_to_master,
                                                                fraction_of_benders_cuts_to_master = 0.5,
                                                                reuse_dcglp=reuse_dcglp)

                                                                master = Master(data; model = update_master_model!, optimizer = mip_optimizer)
                            typical_oracles = [CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer); CFLKnapsackOracle(data, master; model = update_sub_model!, optimizer = optimizer)]
                            disjunctive_oracle = SplitOracle(master, Tuple(typical_oracles), oracle_param)
                            env = SpecializedBendersSeq(master, disjunctive_oracle; param = specialized_benders_param)
                        log = solve!(env)
                        @test env.termination_status == Optimal() ? isapprox(mip_opt_val, env.obj_value, atol=1e-5) : false
                    end
                end
            end
        end
    end
end
