module Base

struct Constant
    γ :: Float64 # gas specific ratio
    pr:: Float64 # Prandtl number
end
include("Sutherland.jl")
end