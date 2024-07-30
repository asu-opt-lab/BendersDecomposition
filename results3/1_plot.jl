using CSV,DataFrames,Plots

function pplot(files)
    iters = []
    LBs = []
    for k in [200]
        # for i in [2,5,6,7,8,9,10]
        for i in 1:10
            filename = "$(files)/result_f$(k)-c$(k)-r5.0-p$(i).csv"
            data = CSV.read("$(filename)", DataFrame)
            iter = data.iter[end]
            LB = data.LB[end] 
            push!(iters, iter)
            push!(LBs, LB)
        end
    end
    return iters, LBs
end


x = range(1,10, length=10)
files = "results3/Ordinary_4hr"
iters, LBs = pplot(files)
figure = scatter(x, LBs, label="Standard Benders Cut", xlabel="instance", ylabel="Lower Bound", title="Results (TimeLimit=1000s)", dpi=300, legend=:right)

# files = "results3/Ordinary_4hr"
# iters, LBs = pplot(files)
# scatter!(figure, x, LBs, label="Standard Benders Cut")

# files = "results3/Split_all_L1_cplex"
# iters, LBs = pplot(files)
# scatter!(plot, x, LBs, label="Split_all")

# files = "results3/Split_all_L1_iter"
# iters, LBs = pplot(files)
# scatter!(plot, x, LBs, label="Split_all_iter100")

files = "results3/Split_all_L1_iter2_1hr"
iters, LBs = pplot(files)
scatter!(figure, x, LBs, label="Disjuctive Benders Cut")

# files = "results3/Advanced_4hr"
# iters, LBs = pplot(files)
# scatter!(figure, x, LBs, label="Standard Benders Cut")

savefig(figure, "results3/plot_200.png")
