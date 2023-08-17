module TransitionAnalysis

# Write your package code here.
# using Reexport

include("BasicFunction/TAbase.jl")
import .Base

@info constatns=Base.Constant(
    1.4, #gas specific ratio
    0.72 #Prandtl number
    )

include("EmpericalRelation/EmpericalRelation.jl")
import .EmpericalRelation

end
