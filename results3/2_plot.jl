using CSV,DataFrames,Plots

function pplot(files)
    k = 700
    i = 1
    filename = "$(files)/result_f$(k)-c$(k)-r5.0-p$(i).csv"
    data = CSV.read("$(filename)", DataFrame)
    iter = data.iter
    LB = data.LB

    return iter, LB
end


files = "results3/Advanced_lp_2hr"
iters1, LBs1 = pplot(files)
figure = plot(iters1, LBs1, label="Enhanced Benders Cut", xlabel="# of iteration", ylabel="Lower Bound", title="Solving Process", marker=:rect,markersize=2, dpi=300, legend=:right)

# files = "results3/Advanced_4hr"
# iters2, LBs2 = pplot(files)
# plot!(figure, iters2, LBs2, label="Enhanced Benders Cut",marker=:rect,markersize=2)

files = "results3/Split_all_L1_iter2_2hr"
iters3, LBs3 = pplot(files)
plot!(figure, iters3, LBs3, label="Disjuctive Benders Cut",marker=:rect,markersize=2)

savefig(figure, "results3/plot_700_1.png")
