module TransitionAnalysis

module Consts
    const γ=1.4::Float64    # gas specific ratio
    const pr=0.72::Float64  # Prandtl number
end
# Write your package code here.
# export 
# Modules
# nothing
# Functions
# nothing


# Fundamental functions
include("BasicFunction/TAbase.jl")

# Functions on Bondary Layer physics
include("BoundaryLayer/BoundaryLayer.jl")

# EmpericalRelation (correlations)
include("EmpericalRelation/EmpericalRelation.jl")

end
