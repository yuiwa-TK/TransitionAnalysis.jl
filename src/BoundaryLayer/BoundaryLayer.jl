module BoundaryLayer
using ..TransitionAnalysis
deriv1d=TransitionAnalysis.deriv1d

"""
    boundarylayer_thickness(u, y; uinf, α)
"""
function boundarylayer_thickness(u,y;uinf,α)
    utarget=uinf
    lmax = length(y)
    for l in 2:lmax
        # @show l,abs(u[l]- utarget)/utarget
        if (u[l])/uinf>α
            lp = l
            lm = l-1
            up = u[lp]
            um = u[lm]

            du1 = utarget-um
            du2 = up - utarget
            duu = up-um

            yp = y[lp]
            ym = y[lm]

            delta = (yp*du1 + ym*du2)/duu

            return delta
        end
    end
    println("delta $α was not found." )
    println("u[lmax]/uinf", u[lmax]/uinf)
    return nothing
end
"""
    delta99(u,y,uinf) returns 99% thickness of the boundary layer

- This is the alias function of boundarylayer_thickness(u,y;uinf,α=0.99)
"""
delta99(u,y,uinf)=boundarylayer_thickness(u,y;uinf,α=0.99)

"""

"""
function displacement_thickness(ρu, y; Mach,ρinf=1.0)
    lmax=length(y)
    sum=0.0
    minv=1.0f0/Mach/ρinf
    for l in 2:lmax
        fm=1.0f0 - ρu[l-1]*minv
        fp=1.0f0 - ρu[l]*minv
        sum= sum+ 0.5*(fp+fm)*(y[l]-y[l-1])
    end
    return sum
end
function displacement_thickness(ρ,u, y; Mach)
    return displacement_thickness(ρ.*u,y;Mach=Mach)
end

function displacement_thickness_incompressible(u,y;Mach)
    lmax=length(y)
    sum=0.0
    minv=1.0f0/Mach
    for l in 2:lmax
        fm=1.0f0 - u[l-1]*minv
        fp=1.0f0 - u[l]*minv
        sum= sum+ 0.5*(fp+fm)*(y[l]-y[l-1])
    end
    return sum
end

function momentum_thickness(ρ,u, y; Mach,ρinf=1.0)
    v = (1.0.-u./(Mach)).*ρ.*u/(ρinf*Mach)
    lmax=length(y)
    sum=0.0
    for j in 2:lmax
        dm=v[j-1]
        dp=v[j]
        sum+= 0.5*(dp+dm)*(y[j]-y[j-1])
    end
    return sum
end

"""
    inflection_point(u,y, uinf; verbose)
    returns the first inflection point.

    u: velocity
    y: distance from the wall
    uinf: free stream speed
    verbose: if >0, the u-y, du/dy-y, d2u/dy2-y plots are shown. (defalut=0)
"""
function inflection_point(u, y, uinf; verbose=0)
    @assert length(u)==length(y)

    d99 = boundarylayer_thickness(u,y; uinf=uinf,α=0.99)
    idend=argmin(abs.(y.-d99*1.5))
    du = deriv1d(u,y, type="2ndcentral")
    d2 = deriv1d(du,y, type="2ndcentral");

    if verbose>0
        p=plot(u,y, label="u")
        p=plot!(p,du./maximum(abs.(du)),y,label="du/max|du|")
        p=plot!(p,d2./maximum(abs.(d2)),y,label="d2u/max|d2u|")
        p=ylims!(0, y[idend])
    end
    idspan=1
    N=idend
    inflection_point=zeros(5)
    # Even if there exists d2max s.t. d2>0 but d2max is small sufficiently, such point is not regarded as an inflection point.
    d2max=maximum(d2)./maximum(abs.(d2))
    if d2max < 0.01
        if verbose>0
            println("no point")
            display(p)
        end
        return [0.0]
    end

    for j in idspan+1:N-idspan
            dp=d2[j+idspan]
            dd=d2[j]
            dm=d2[j-idspan]
            dmh=0.5*(dm+dd)
            dph=0.5*(dp+dd)
            if dph*dmh<0
                # @show j, y[j]
                push!(inflection_point, y[j])
                if verbose>0
                    println("j:",j, " y:",y[j], " dp*dm:", dp*dm, " d2u:",d2[j])
                    p=scatter!(p, [d2[j]./maximum(abs.(d2))], [y[j]], color="black",  label=false)
                    p=scatter!(p, [u[j]], [y[j]], color="black",  label=false, pt=:diamond)
                end
            end
    end
    if verbose>0
        display(p)
    end
    return Float32.(inflection_point)
end

"""
    Fjortoft_unstable_inflection_point(u,y, uinf; verbose)
    returns the unstable inflection point based on Fjotoft's theorem. [U"(U-Us)<0]

    u: velocity
    y: distance from the wall
    uinf: free stream speed
    verbose: if >0, the u-y, du/dy-y, d2u/dy2-y plots are shown. (defalut=0)
"""
function Fjortoft_unstable_inflection_point(u, y, uinf; verbose=0)
    @assert length(u)==length(y)

    d99 = boundarylayer_thickness(u,y; uinf=uinf,α=0.99)
    idend=argmin(abs.(y.-d99*1.5))
    du = deriv1d(u,y, type="2ndcentral")
    d2 = deriv1d(du,y, type="2ndcentral");

    if verbose>0
        p=plot(u,y, label="u")
        p=plot!(p,du./maximum(abs.(du)),y,label="du/max|du|")
        p=plot!(p,d2./maximum(abs.(d2)),y,label="d2u/max|d2u|")
        p=ylims!(0, y[idend])
    end
    inflection_point=zeros(20)

    # Even if there exists d2max s.t. d2>0 but d2max is small sufficiently, such point is not regarded as an inflection point.
    d2max=maximum(d2)./maximum(abs.(d2))
    if d2max < 0.01
        if verbose>0
            println("no point")
            display(p)
        end
        return inflection_point
    end


    N=idend
    Nip=0
    idspan=1
    bandwidth=idspan*2
    for j in bandwidth+1:idspan:N-bandwidth
            dp=sum( d2[j:j+bandwidth])/(bandwidth+1)
            dd=sum( d2[j-idspan:j+idspan])/(bandwidth+1)
            dm=sum( d2[j-bandwidth:j])/(bandwidth+1)
            dmh=0.5*(dm+dd)
            dph=0.5*(dp+dd)
            dy=0.5(y[j+1]-y[j-1])
            tangentval= ( 8.0*( d2[j+1]-d2[j-1]) - (d2[j+2]-d2[j-2]))/(12*dy)./maximum(abs.(d2))
        
            if dph*dmh<0 && abs(tangentval) > tan(π/2/10)
                Us=u[j]
                cond=false
                for jj in eachindex(y)
                    cond= d2[j]*(u[jj]-Us)<0.0
                    if cond && abs(tangentval) > tan(π/2/6)#15deg
                        Nip+=1
                        inflection_point[Nip]=y[j]
                        break
                    end
                end
                if verbose>0
                    println("j:",j, " y:",y[j], " dp*dm:", dp*dm, " d2u:",d2[j])
                    p=scatter!(p, [d2[j]./maximum(abs.(d2))], [y[j]], color="black",  label=false)
                    p=scatter!(p, [u[j]], [y[j]], color="black",  label=false, pt=:diamond)
                end
            end
    end
    if verbose>0
        display(p)
    end
    return inflection_point
end





end