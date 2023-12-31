using ..TransitionAnalysis

"""
    cf_lami__incompressible(Rex)
returns an emperical relation of *incompressible* *turbulent* boundary layer's skin friction coefficient

## Reference
White, "Viscous Fluid Flow"(2006), pp436, Fig 6-20
"""
function cf_lami__incompressible(Rex)
    return 0.664/sqrt(Rex)
end

"""
    cf_lami__incompressible(Rex)
returns an emperical relation of *compressible* *turbulent* boundary layer's skin friction coefficient

## Reference
White, "Viscous Fluid Flow"(2006), pp517, Eq 7-39
"""
function cf_lami__compressible(Rex,Twall, Te=1.0)
    cw_root = sqrt((Twall/Te)^(-1/3))
    return 0.664*cw_root/sqrt(Rex)
end

"""
    cf_turb__incompressible(Rex)
returns an emperical relation of incompressible *turbulent* boundary layer's skin friction coefficient

## Reference
White, "Viscous Fluid Flow"(2006), pp436, eq(6-78)
"""
function cf_turb__incompressible(Rex)
    a = (log(0.06*Rex))*(log(0.06*Rex))
    return 0.455/a
end

"""
    cf_turb__compressible_WhiteCristoph(Rex, Twall, Mach; verbose=0)
computes the theoretical friction coef. of the compressible boundary layer on a flat plate, proposed White&Christoph.

## Reference
White, "Viscous Fluid Flow"(2006), pp558, eq(7-133)
"""
function cf_turb__compressible_WhiteCristoph(Rex, Twall, Mach; verbose=0)
    c1=0.455
    γ  = TransitionAnalysis.Consts.γ
    γ1 = γ-1
    Te = 1.0
    μe = 1.0
    if verbose==1
        @info c1,γ,Te,μe
    end
    Taw = TransitionAnalysis.recovery_temp(Mach,"turbulent")
    # eq.7-107 
    a = sqrt( 0.5*γ1*Mach*Mach*Te/Twall )
    b = Taw/Twall-1
    
    # eq.7-119
    A = (2*a*a-b)/sqrt(b^2+4a^2)
    B = b / sqrt(b*b + 4*a*a)

    # eq.7-130
    S = sqrt(Taw/Te-1)/ (asin(A)+asin(B))
    
    # eq.7-133
    Fc = S*S
    μw = TransitionAnalysis.Sutherland(Twall)
    F_Rex = (μe/μw)*sqrt(Te/Twall)/S
    Cf = cf_turb__incompressible(Rex*F_Rex)/Fc

    # @info "Eq 7-132"
    # Cf = 0.455/ (S*S)/ (log(0.06/S*Rex*μw/μe*sqrt(Te/Twall)))^2
    return Cf
end
"""
    Cf_conv_comp2incomp(Cf_comp, Twall, T∞)
convert the friction coefficient of adiabatic compressible turbulent boundary layer into incompressible turbulent boundary laeyr.
"""
function Cf_conv_comp2incomp(Cf_comp; Twall, T∞=1.0)
    r = Twall/T∞
    α = (r-1)/ sqrt(r*(r-1))
    p = asin(α)*asin(α)
    return Cf_comp*(r-1)/p
end