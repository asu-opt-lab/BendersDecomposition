# Load order matters: oracleDisjunctiveSplit.jl defines SplitOracle before the
# normalization-specific files add methods for it.
include("common.jl")
include("split_disjunction.jl")
include("cut_pool.jl")
include("lifting_strengthening.jl")
include("lp_distance_normalization.jl")
include("reverse_polar_normalization.jl")
include("logging.jl")
include("solve_dcglp.jl")
