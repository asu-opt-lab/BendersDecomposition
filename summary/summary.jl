using CSV,DataFrames

files = readdir("results/gurobi")

cut_strategy = settings["cut_strategy"]
SplitCGLPNormType = settings["SplitCGLPNormType"]
SplitSetSelectionPolicy = settings["SplitSetSelectionPolicy"]
StrengthenCutStrategy = settings["StrengthenCutStrategy"]
SplitBendersStrategy = settings["SplitBendersStrategy"]

df = DataFrame(name=String[], iter=Int[], LB = Float64[], master_time=Float64[], sub_time=Float64[], time=Float64[])
for file in files
    filename = basename(file)
    contains = occursin("_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", filename)
    if contains
        data = CSV.read("results/gurobi/$(filename)", DataFrame)
        if nrow(data) != 0
            iter = data.iter[end]
            LB = data.LB[end]
            master_time = sum(data.master_time)
            sub_time = sum(data.sub_time)
            time = master_time + sub_time
            new_row = (filename, iter, LB, master_time, sub_time, time)
            push!(df, new_row)
        end
    end
end


CSV.write("summary/gurobi_summary_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)