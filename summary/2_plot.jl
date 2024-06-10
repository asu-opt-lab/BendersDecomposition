using CSV,DataFrames,Plots

files = readdir("results/Gurobi")

# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"

# "L1GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "LINFGAMMANORM"

# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "MOST_FRAC_INDEX"

# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"

# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "NO_SPLIT_BENDERS_STRATEGY"

function pplot(StrengthenCutStrategy, SplitBendersStrategy)
    iters = []
    LBs = []
    # for k in [100,200,300,500,700,1000]
    for k in [100]
        for i in 1:10
            filename = "results/Gurobi/result_f$(k)-c$(k)-r5.0-p$(i)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv"
            if isfile(filename)
                # @info "here"
                data = CSV.read("$(filename)", DataFrame)
                iter = data.iter[end]
                LB = data.LB[end]
                push!(iters, iter)
                push!(LBs, LB)
            else
                push!(iters, 0)
                push!(LBs, -Inf)
            end
        end
    end
    return iters, LBs
end


x = range(1,10, length=10)
StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"
SplitBendersStrategy = "NO_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
plot = scatter(x, LBs, label="LB_Pure_No", xlabel="instance", ylabel="LB", title="LB", dpi=300, legend=:bottomleft, color=:red, markersize=10)

StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"
SplitBendersStrategy = "ALL_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
scatter!(plot, x, LBs, label="LB_Pure_All", xlabel="instance", ylabel="LB", title="LB", markersize=8)

StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"
SplitBendersStrategy = "TIGHT_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
scatter!(plot, x, LBs, label="LB_Pure_Tight", xlabel="instance", ylabel="LB", title="LB", markersize=6)

StrengthenCutStrategy = "SPLIT_STRENGTHEN_CUT_STRATEGY"
SplitBendersStrategy = "NO_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
scatter!(plot, x, LBs, label="LB_Strengthen_No", xlabel="instance", ylabel="LB", title="LB", markersize=5)

StrengthenCutStrategy = "SPLIT_STRENGTHEN_CUT_STRATEGY"
SplitBendersStrategy = "ALL_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
scatter!(plot, x, LBs, label="LB_Strengthen_All", xlabel="instance", ylabel="LB", title="LB", markersize=4)

StrengthenCutStrategy = "SPLIT_STRENGTHEN_CUT_STRATEGY"
SplitBendersStrategy = "TIGHT_SPLIT_BENDERS_STRATEGY"
iters, LBs = pplot(StrengthenCutStrategy, SplitBendersStrategy)
scatter!(plot, x, LBs, label="LB_Strengthen_Tight", xlabel="instance", ylabel="LB", title="LB", markersize=3)

savefig(plot, "summary/plot_strength_benders_100.png")


