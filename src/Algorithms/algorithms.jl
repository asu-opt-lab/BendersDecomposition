export run_Benders


abstract type AbstractCutStrategy end
struct OrdinaryCutStrategy <: AbstractCutStrategy end
const ORDINARY_CUTSTRATEGY = OrdinaryCutStrategy()
struct KNCutStrategy <: AbstractCutStrategy end
const KN_CUTSTRATEGY = KNCutStrategy()
struct SplitCutStrategy <: AbstractCutStrategy end
const SPLIT_CUTSTRATEGY = SplitCutStrategy()

# abstract type AbstractNormType end
# struct StandardNorm <: AbstractNormType end
# const STANDARDNORM = StandardNorm()
# abstract type GammaNorm <: AbstractNormType end
# struct L1GammaNorm <: GammaNorm end
# const L1GAMMANORM = L1GammaNorm()
# struct L2GammaNorm <: GammaNorm end
# const L2GAMMANORM = L2GammaNorm()
# struct LInfGammaNorm <: GammaNorm end
# const LINFGAMMANORM = LInfGammaNorm()

abstract type AbstractSplitSetSelectionPolicy end
struct MostFracIndex <: AbstractSplitSetSelectionPolicy end
const MOST_FRAC_INDEX = MostFracIndex()
struct RandomIndex <: AbstractSplitSetSelectionPolicy end
const RANDOM_INDEX = RandomIndex()

abstract type AbstractSplitBendersPolicy end
struct NoSplitBendersStrategy <: AbstractSplitBendersPolicy end
const NO_SPLIT_BENDERS_STRATEGY = NoSplitBendersStrategy()
struct AllSplitBendersStrategy <: AbstractSplitBendersPolicy end
const ALL_SPLIT_BENDERS_STRATEGY = AllSplitBendersStrategy()
struct TightSplitBendersStrategy <: AbstractSplitBendersPolicy end
const TIGHT_SPLIT_BENDERS_STRATEGY = TightSplitBendersStrategy()

abstract type AbstractSplitStengtheningPolicy end
struct SplitPureCutStrategy <: AbstractSplitStengtheningPolicy end
const SPLIT_PURE_CUT_STRATEGY = SplitPureCutStrategy()
struct SplitStrengthenCutStrategy <: AbstractSplitStengtheningPolicy end
const SPLIT_STRENGTHEN_CUT_STRATEGY = SplitStrengthenCutStrategy()

struct AlgorithmParams
    cut_strategy::Union{AbstractCutStrategy, Nothing}
    SplitCGLPNormType::Union{AbstractNormType, Nothing}
    SplitSetSelectionPolicy::Union{AbstractSplitSetSelectionPolicy, Nothing}
    StrengthenCutStrategy::Union{AbstractSplitStengtheningPolicy, Nothing}
    SplitBendersStrategy::Union{AbstractSplitBendersPolicy, Nothing} 
end

function run_Benders(

)
    @error "wrong type set" 
end

function generate_cut(

)
    @error "wrong type set"
end

include("Loop/loop.jl")
include("Cut/cut.jl")
