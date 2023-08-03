"""
    cf_incompressible(Rex)
returns an emperical relation of incompressible boundary layer's skin friction coefficient

## Reference
White, "Viscous Fluid Flow"(2006), pp436, eq(6-78)
"""
function cf_incompressible(Rex)
    return 0.455/(ln(0.06*Rex))/(ln(0.06*Rex))
end

"""
    cf_compressible_WhiteCristoph(Rex, Twall, Mach; verbose=0)
computes the theoretical friction coef. of the compressible boundary layer on a flat plate, proposed White&Christoph.

## Reference
White, "Viscous Fluid Flow"(2006), pp558, eq(7-133)
"""
function cf_compressible_WhiteCristoph(Rex, Twall, Mach; verbose=0)
    c1=0.455
    γ = 1.4
    γ1 = 0.4
    Te = 1.0
    μe = 1.0
    if verbose==1
        @info c1,γ,Te,μe
    end

    # eq.7-107 
    a = sqrt(0.5*γ1*Ma^2*Te/Tw)
    b = sqrt(Taw/Tw)
    
    # eq.7-119
    A = 2(a*a-b)/sqrt(b^2+4a^2)
    B = b / sqrt(b*b + 4*a*a)

    # eq.7-130
    S = sqrt(Taw/Te-1)/ (asin(A)+asin(B))
    
    # eq.7-133
    Fc = S*S
    μw = TAbase.Sutherland(Tw)
    F_Rex = (μe/μw)*sqrt(Te/Tw)/(S*S)    
    Cf = cf_incompressible(Rex*F_Rex)/Fc

    return Cf
end


