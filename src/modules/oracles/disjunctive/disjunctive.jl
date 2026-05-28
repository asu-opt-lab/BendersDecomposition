# Load order matters: strategy types are needed by split_oracle.jl; split_oracle.jl
# defines the const aliases used by the strategy-specific files.
include("common.jl")
include("strategy.jl")
include("split_oracle.jl")
include("split_disjunction.jl")
include("cut_pool.jl")
include("lifting_strengthening.jl")
include("distance_norm.jl")
include("simplex_norm.jl")
include("vertical_reverse_polar.jl")
include("directional_polar.jl")
include("logging.jl")
include("solve_dcglp.jl")
