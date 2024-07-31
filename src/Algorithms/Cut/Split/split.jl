# include("LGamma.jl")
# include("LStandard.jl")
# include("LGamma_ordinary.jl")
include("LGamma_ordinary_.jl")
# include("LGamma_advance.jl")
# include("LGamma_kn.jl")
include("LStandard_ordinary.jl")

function update_UB!(UB,_sx,g₁,g₂,t̂,::L1GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], Inf)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::L2GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 2)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::LInfGammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 1)) end
