using CSV,DataFrames

files = readdir("results4/Split_all_L1_iter50_2hr_")

df = DataFrame(name=String[], iter=Int[], LB = Float64[], UB=Float64[], gap=Float64[], time=Float64[])
for file in files
    filename = basename(file)
    @info filename
    contains = occursin(".csv", filename)
    if contains
        data = CSV.read("results4/Split_all_L1_iter50_2hr_/$(filename)", DataFrame)
        iter = data.iter[end-1]
        LB = data.LB[end-1]
        UB = data.UB[end-1]
        gap = data.gap[end-1]
        master_time = sum(data.master_time[1:end-1])
        sub_time = sum(data.sub_time[1:end-1])
        time = master_time + sub_time
        new_row = (filename, iter, LB, UB, gap, time)
        push!(df, new_row)
    end
end

CSV.write("summary/results_new.csv", df)


# files = readdir("results4/Split_all_L1_iter50_2hr")

# df = DataFrame(name=String[], iter=Int[], LB = Float64[], UB=Float64[], gap=Float64[], time=Float64[])
# for file in files
#     filename = basename(file)
#     @info filename
#     contains = occursin(".csv", filename)
#     if contains
#         data = CSV.read("results4/Split_all_L1_iter50_2hr/$(filename)", DataFrame)
#         iter = data.iter[end-1]
#         LB = data.LB[end-1]
#         UB = data.UB[end-1]
#         gap = data.gap[end-1]
#         master_time = sum(data.master_time[1:end-1])
#         sub_time = sum(data.sub_time[1:end-1])
#         time = master_time + sub_time
#         new_row = (filename, iter, LB, UB, gap, time)
#         push!(df, new_row)
#     end
# end

# CSV.write("summary/results.csv", df)