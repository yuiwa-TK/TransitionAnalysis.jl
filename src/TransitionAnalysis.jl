module TransitionAnalysis

module Consts
    const γ=1.4::Float64    # gas specific ratio
    const pr=0.72::Float64  # Prandtl number
end
# Write your package code here.
include("BasicFunction/TAbase.jl")
include("EmpericalRelation/EmpericalRelation.jl")

end
