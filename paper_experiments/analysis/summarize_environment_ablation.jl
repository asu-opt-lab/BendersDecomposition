include("common_analysis.jl")

df = raw_csv("03_environment_ablation_summary.csv")
benders = filter(:config_type => ==("benders"), df)

summary = combine(
    groupby(benders, [:problem, :env, :oracle]),
    :solved => sum_bool => :solved,
    :solved => length => :runs,
    :total_time => median_finite => :median_time,
    :iterations => median_finite => :median_iterations,
    :oracle_share => median_finite => :median_oracle_share,
    :final_gap => median_finite => :median_final_gap,
)

write_processed("table_environment_ablation.csv", summary)

