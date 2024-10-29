export SequentialBenders, solve!

struct SequentialBenders <: AbstractBendersAlgorithm
    data::Any
    master::AbstractMasterProblem
    sub::AbstractSubProblem
    cut_strategy::CutGenerationStrategy
    params::BendersParams
    iteration_data::DataFrame  
end


function SequentialBenders(data, master, sub, cut_strategy, params)
    iteration_data = DataFrame(
        iter = Int[], 
        LB = Float64[], 
        UB = Float64[], 
        gap = Float64[], 
        master_time = Float64[], 
        sub_time = Float64[], 
        elapsed_time = Float64[]
    )
    return SequentialBenders(data, master, sub, cut_strategy, params, iteration_data)
end

function solve!(algo::SequentialBenders)
    UB, LB = Inf, -Inf
    iter = 0
    start_time = time()
    master_time, sub_time = 0.0, 0.0
    
    while true
        iter += 1
        
        # Solve master problem and record time
        master_start = time()
        solve_master!(algo.master)
        master_time += time() - master_start
        
        # Solve sub problem and record time
        sub_start = time()
        solve_sub!(algo.sub, algo.master)
        sub_time += time() - sub_start
        
        cuts, sub_obj_value = generate_cuts(algo, algo.cut_strategy)
        @constraint(algo.master.model, cuts .<= 0)

        # Update bounds
        LB = algo.master.obj_value
        UB_temp = sum(algo.data.fixed_costs .* algo.master.x_value) + sub_obj_value
        UB = min(UB, UB_temp)
        
        # Calculate gap and elapsed time
        gap = (UB - LB) / UB * 100
        elapsed_time = time() - start_time
        
        # Store iteration data
        push!(algo.iteration_data, (
            iter = iter,
            LB = LB,
            UB = UB,
            gap = gap,
            master_time = master_time,
            sub_time = sub_time,
            elapsed_time = elapsed_time
        ))
        
        # Print iteration information

        # @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.2f%% | Time: (M: %6.2f, S: %6.2f) | Elapsed: %6.2f\n",
        #            iter, LB, UB, gap, master_time, sub_time, elapsed_time)

        # Check termination criteria
        if gap < algo.params.gap_tolerance || elapsed_time > algo.params.time_limit
            break
        end
    end

    # Print final time breakdown
    # @printf("Total time: %.2f s (Master: %.2f s, Sub: %.2f s)\n", 
    #         time() - start_time, master_time, sub_time)

    return algo.iteration_data
end

# solve master problem
function solve_master!(master::AbstractMasterProblem)
    optimize!(master.model)
    master.obj_value = objective_value(master.model)
    master.x_value = value.(master.var[:x])
end

# solve sub problem
function solve_sub!(sub::AbstractSubProblem, master::AbstractMasterProblem)

    set_normalized_rhs.(sub.fixed_x_constraints, master.x_value)
    optimize!(sub.model)
end

function solve_sub!(sub::KnapsackUFLPSubProblem, master::AbstractMasterProblem)
    # no need to solve sub problem (knapsack)
end



