using CSV,DataFrames,Plots

function pplot(files)
    iters = []
    LBs = []
    for k in [700]
        for i in [1,2,4,5,6,7,8,9,10]
        # for i in 1:10
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
files = "results2/Advanced"
iters, LBs = pplot(files)
plot = scatter(x, LBs, label="Advance", xlabel="instance", ylabel="LB", title="LB", dpi=300, legend=:bottomleft)

# files = "results2/Ordinary"
# iters, LBs = pplot(files)
# scatter!(plot, x, LBs, label="Ordinary")

# files = "results2/Split"
# iters, LBs = pplot(files)
# scatter!(plot, x, LBs, label="Split")

files = "results2/Split_all"
iters, LBs = pplot(files)
scatter!(plot, x, LBs, label="Split_all")

files = "results2/Split_all_mip"
iters, LBs = pplot(files)
scatter!(plot, x, LBs, label="Split_all_mip")

savefig(plot, "results2/plot_700_.png")
