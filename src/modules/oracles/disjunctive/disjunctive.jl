# Load order matters: normalization types are needed by split_oracle.jl; split_oracle.jl
# defines the generic oracle type used by the normalization-specific files.
include("common.jl")
include("normalization.jl")
include("split_oracle.jl")
include("split_disjunction.jl")
include("cut_pool.jl")
include("lifting_strengthening.jl")
include("lp_distance_normalization.jl")
include("vertical_reverse_polar_normalization.jl")
include("directional_reverse_polar_normalization.jl")
include("logging.jl")
include("solve_dcglp.jl")
