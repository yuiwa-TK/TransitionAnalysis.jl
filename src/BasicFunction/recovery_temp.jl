using ..TransitionAnalysis
function recovery_temp_laminar(mach)
    γ  = TransitionAnalysis.Consts.γ
    Pr = TransitionAnalysis.Consts.pr
    r  = Pr^(1/3)

    recovery_temp_laminar = 1.0 + 0.5*r*(γ-1)*mach*mach

    return recovery_temp_laminar
end

function recovery_temp_turbulent(mach)
    γ  = TransitionAnalysis.Consts.γ
    Pr = TransitionAnalysis.Consts.pr
    r  = Pr^(1/2)

    recovery_temp_turb = 1.0 + 0.5*r*(γ-1)*mach*mach

    return recovery_temp_turb
end

function recovery_temp(mach, mode::AbstractString)
    if mode=="laminar"
        return recovery_temp_laminar(mach)
    elseif mode=="turbulent"
        return recovery_temp_turbulent(mach)
    end
end