include("common_analysis.jl")

df = raw_csv("01_correctness_summary.csv")
benders = filter(:config_type => ==("benders"), df)

summary = combine(
    groupby(benders, [:problem, :oracle]),
    :solved => sum_bool => :solved,
    :solved => length => :runs,
    :obj_error => maximum_finite => :max_obj_error,
    :total_time => median_finite => :median_time,
)

write_processed("table_correctness.csv", summary)

