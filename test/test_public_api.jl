using Test
using BendersX

@testset "Public API Contract" begin
    exported = [
        :AbstractData,
        :Master,
        :BendersSeq, :BendersSeqInOut, :SpecializedBendersSeq, :BendersBnB,
        :ClassicalOracle, :UnifiedOracle, :ParetoOracle, :SeparableOracle,
        :SplitOracle,
        :solve!,
        :BasicOracleParam, :ClassicalOracleParam, :UnifiedOracleParam, :ParetoOracleParam,
        :SeparableOracleParam, :DcglpParam,
        :SplitOracleParam, :LpDistanceNormalization, :ReversePolarNormalization,
        :BendersSeqParam, :BendersSeqInOutParam, :BendersBnBParam, :SpecializedBendersSeqParam,
        :LazyCallback, :UserCallback, :NoUserCallback, :UserCallbackParam,
        :RootNodePreprocessing, :NoRootNodePreprocessing, :DisjunctiveRootNodePreprocessing,
        :TerminationStatus, :NotSolved, :TimeLimit, :Optimal, :InfeasibleOrNumericalIssue,
        :GBCBoundType, :UpperBound, :LowerBound, :Fixed,
        :RandomFractional, :MostFractional, :LargestFractional,
        :NoDisjunctiveCuts, :AllDisjunctiveCuts, :DisjunctiveCutsSmallerIndices,
        :set_parameter!, :set_core_point!,
        :update_master_model!, :update_sub_model!,
        :CFLPData, :UFLPData, :SCFLPData, :SNIPData,
        :CFLKnapsackOracle, :CFLKnapsackOracleParam, :UFLKnapsackOracle, :UFLKnapsackOracleParam,
        :read_GK_data, :read_cfl_file, :read_cflp_benchmark_data,
        :read_uflp_benchmark_data, :read_Simple_data,
        :read_stochastic_capacited_facility_location_problem, :read_snip_data,
    ]

    public_only = [
        :AbstractBendersEnv, :AbstractBendersSeq, :AbstractBendersBnB,
        :AbstractMaster, :AbstractOracle, :AbstractOracleParam, :AbstractTypicalOracle,
        :AbstractDisjunctiveOracle, :AbstractSplitOracle,
        :AbstractRootNodePreprocessing,
        :AbstractDisjunctiveNormalization,
        :SplitIndexSelectionRule, :DisjunctiveCutsAppendRule,
        :generate_cuts, :add_normalization_constraint!, :update_dcglp_upper_bound_and_gap!,
        :Hyperplane, :aggregate, :evaluate_violation,
        :select_top_fraction, :hyperplanes_to_expression, :add_constraints,
        :copy_variables!, :var_from_tuple, :transfer_scaled_linear_rows_and_bounds_with_types!,
        :infeasibility_report,
        :TimeLimitException, :UnexpectedModelStatusException, :UndefError,
        :AlgorithmException, :UnsupportedModelException,
    ]

    private_symbols = [
        :lazy_callback, :user_callback, :root_node_processing!,
        :get_sec_remaining, :record_iteration!, :update_upper_bound_and_gap!,
        :is_terminated, :check_lb_improvement!, :print_iteration_info, :to_dataframe,
        :calculate_KP_value,
        :AbstractLoopState, :AbstractLoopLog, :AbstractLoopParam,
        :AbstractBendersSeqState, :AbstractBendersSeqLog, :AbstractBendersSeqParam,
        :AbstractBendersBnBState, :AbstractBendersBnBLog, :AbstractBendersBnBParam,
        :BendersSeqState, :BendersSeqLog, :BendersBnBState, :BendersBnBLog,
        Symbol("Abstract" * "Dcglp" * "Oracle"),
        Symbol("Dcglp" * "Oracle"),
        Symbol("Dcglp" * "OracleParam"),
    ]

    visible_names = Set(names(BendersX))

    for sym in exported
        @test Base.isexported(BendersX, sym)
        @test Base.ispublic(BendersX, sym)
        @test sym in visible_names
    end

    for sym in public_only
        @test !Base.isexported(BendersX, sym)
        @test Base.ispublic(BendersX, sym)
        @test sym in visible_names
    end

    for sym in private_symbols
        @test !Base.ispublic(BendersX, sym)
        @test !(sym in visible_names)
    end
end
