using CSV,DataFrames,Plots

function pplot()
    iters = []
    LBs = []
    df = DataFrame(Split_time = Float64[], Split_LB = Float64[], Split_iter = Int[], Advanced_time=Float64[], Advanced_LB = Float64[], Advanced_iter=Int[], diff = Float64[])
    for k in [1000]
        # for i in [1,2,3,5,6,7,8,9]
        for i in 1:10
            iter_ = 0
            all_time = 0
            filename = "results3/Split_all_L1_iter2_2hr/result_f$(k)-c$(k)-r5.0-p$(i).csv"
            data = CSV.read("$(filename)", DataFrame)
            std_LB = data.LB[end]
            std_time = sum(data.master_time[1:end-1] + data.sub_time[1:end-1])
            std_iter = data.iter[end]
            filename = "results3/Advanced_lp_2hr/result_f$(k)-c$(k)-r5.0-p$(i).csv"
            data = CSV.read("$(filename)", DataFrame)
            # LBs = data.LB
            ad_LB = data.LB[end]
            ad_time = sum(data.master_time[1:end-1] + data.sub_time[1:end-1])
            ad_iter = data.iter[end]
            # for j in eachindex(LBs)
            #     if LBs[j] >= std_LB
            #         iter_ = j
            #         all_time = sum(data.sub_time[1:j]) + sum(data.sub_time[1:j])
            #         break
            #     end
            # end
            diff = std_LB - ad_LB
            new_row = (std_time, std_LB, std_iter, ad_time, ad_LB, ad_iter, diff)
            push!(df, new_row)
        end
    end
    return df
end

df = pplot()
CSV.write("results3/data_2hr_1000.csv", df)