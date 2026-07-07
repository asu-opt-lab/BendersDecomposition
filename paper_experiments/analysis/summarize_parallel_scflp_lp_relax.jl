include("common_analysis.jl")

df = raw_csv("05_parallel_scflp_lp_relax_summary.csv")
benders = filter(:config_type => ==("benders_lp_relax"), df)

summary = combine(
    groupby(benders, [:env, :oracle, :n_customers, :n_scenarios, :julia_threads]),
    :solved => sum_bool => :solved,
    :solved => length => :runs,
    :total_time => median_finite => :median_time,
    :master_time => median_finite => :median_master_time,
    :oracle_time => median_finite => :median_oracle_time,
    :iterations => median_finite => :median_iterations,
    :oracle_share => median_finite => :median_oracle_share,
)

summary.oracle_speedup = fill(NaN, nrow(summary))
summary.oracle_efficiency = fill(NaN, nrow(summary))
for sub in groupby(summary, [:env, :oracle, :n_customers, :n_scenarios])
    t1_rows = filter(:julia_threads => ==(1), sub)
    isempty(t1_rows) && continue
    t1 = t1_rows.median_oracle_time[1]
    for idx in parentindices(sub)[1]
        summary.oracle_speedup[idx] = t1 / summary.median_oracle_time[idx]
        summary.oracle_efficiency[idx] = summary.oracle_speedup[idx] / summary.julia_threads[idx]
    end
end

write_processed("table_parallel_scflp_lp_relax.csv", summary)
