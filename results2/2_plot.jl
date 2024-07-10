using CSV,DataFrames,Plots

function pplot()
    iters = []
    LBs = []
    df = DataFrame(Split_time = Float64[], Advanced_iter=Int[], Advanced_time=Float64[])
    for k in [700]
        for i in [1,2,4,5,6,7,8,9,10]
        # for i in 1:10
            iter_ = 0
            all_time = 0
            filename = "results2/Split_all/result_f$(k)-c$(k)-r5.0-p$(i).csv"
            data = CSV.read("$(filename)", DataFrame)
            std_LB = data.LB[2]
            std_time = data.sub_time[2]
            filename = "results2/Advanced/result_f$(k)-c$(k)-r5.0-p$(i).csv"
            data = CSV.read("$(filename)", DataFrame)
            LBs = data.LB
            for j in eachindex(LBs)
                if LBs[j] >= std_LB
                    iter_ = j
                    all_time = sum(data.sub_time[1:j]) + sum(data.sub_time[1:j])
                    break
                end
            end
            new_row = (std_time, iter_, all_time)
            push!(df, new_row)
        end
    end
    return df
end

df = pplot()
CSV.write("results2/data_700.csv", df)