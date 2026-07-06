include("common_analysis.jl")

df = raw_csv("04_parallel_scflp_summary.csv")
benders = filter(:config_type => ==("benders"), df)

summary = combine(
    groupby(benders, [:env, :oracle, :n_scenarios, :julia_threads]),
    :solved => sum_bool => :solved,
    :solved => length => :runs,
    :total_time => median_finite => :median_time,
    :iterations => median_finite => :median_iterations,
    :oracle_share => median_finite => :median_oracle_share,
)

summary.speedup = fill(NaN, nrow(summary))
summary.efficiency = fill(NaN, nrow(summary))
for sub in groupby(summary, [:env, :oracle, :n_scenarios])
    t1_rows = filter(:julia_threads => ==(1), sub)
    isempty(t1_rows) && continue
    t1 = t1_rows.median_time[1]
    for idx in parentindices(sub)[1]
        summary.speedup[idx] = t1 / summary.median_time[idx]
        summary.efficiency[idx] = summary.speedup[idx] / summary.julia_threads[idx]
    end
end

write_processed("table_parallel_scflp.csv", summary)
