function deriv1d(u,y; type::AbstractString="defalut")
    if type == "2ndcentral"
        return deriv1d_2ndcentral_general(u,y)

    elseif type === "defalut"
        println("2nd central diff is performed.")
        return deriv1d_2ndcentral_general(u,y)
    end
end

function deriv1d_2ndcentral(vec)
    kmax=length(vec)
    dvec=zeros(kmax)
    for k=2:kmax-1
        dvec[k] = 0.5*(vec[k+1]-vec[k-1])
    end
    dvec[1] = (-3vec[1]+4vec[2]-vec[3])*0.5 
    dvec[end] = (3vec[end]-4vec[end-1]+vec[end-2])*0.5
    return dvec
end

function deriv1d_2ndcentral_general(vec, y)
    dydη = deriv1d_2ndcentral(y)
    lmax = length(vec)
    dvecdy = zeros(lmax)
    for l in 2:lmax-1
        dvecdy[l] = 0.5*(vec[l+1]-vec[l-1])/dydη[l]
    end
    dvecdy[1] = (-3vec[1]+4vec[2]-vec[3])*0.5/dydη[1]
    dvecdy[end] = (3vec[end]-4vec[end-1]+vec[end-2])*0.5/dydη[end]
    return dvecdy
end