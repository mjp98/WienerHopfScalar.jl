struct SDuctKernel{T,S} <: WienerHopfKernel
    γ::GammaKernel{T}
    h::S
end
SDuctKernel(k::Number, h::Number) = SDuctKernel(GammaKernel(k), h)

GammaKernel(K::SDuctKernel) = K.γ
GammaKernel(K::SDuctKernel,z) = K.γ(z)
wavenumber(K::SDuctKernel) = wavenumber(GammaKernel(K))
halfheight(K::SDuctKernel, args...) = K.h
halfheights(K::SDuctKernel) = halfheight(K, true), halfheight(K, false)
totalheight(K::SDuctKernel) = sum(halfheights(K)...)

function evaluate(K::SDuctKernel,z)
    γ = GammaKernel(K,z)
    h = halfheight(K)
    return -γ*tanh(h*γ)/2
end

root(K::SDuctKernel,u::Bool,n::Int) = γsinhγ_root(u,n,wavenumber(K),halfheight(K))

function root_residue(K::SDuctKernel,u::Bool,n::Int)
    h = halfheight(K)
    k = wavenumber(K)
    return -2*overγsinhγ_residue(u,n,k,h)*((-1)^n)
end

pole(K::SDuctKernel,u::Bool,n::Int) = coshγ_root(u,n,wavenumber(K),halfheight(K))

function pole_residue(K::SDuctKernel,u::Bool,n::Int)
    x = pole(K,u,n)
    h = halfheight(K)
    k = wavenumber(K)
    γ = K.γ(x)
    return -γ*sinh(γ*h)*overcoshγ_residue(u,n,k,h)/2
end

# ------------------------------------------------------------
# sinhγ
# ------------------------------------------------------------

sinhγ_root(u,n,k,h) = forceabove(sqrt(complex(k^2 - (n*π/h)^2)),Line(),u)
oversinhγ_residue(u,n,k,h) = (im*((-1)^(n+1))*n*π)/(sinhγ_root(u,n,k,h)*h^2)

# ------------------------------------------------------------
# γsinhγ
# ------------------------------------------------------------

γsinhγ_root(u,n,k,h) = sinhγ_root(u,n,k,h)

function overγsinhγ_residue(u,n,k,h)
    n == 0 && return inv(2*h*γsinhγ_root(u,n,k,h))
    return inv(h*γsinhγ_root(u,n,k,h))*((-1)^n)
end

# ------------------------------------------------------------
# coshγ
# ------------------------------------------------------------

coshγ_root(u,n,k,h) = forceabove(sqrt(complex(k^2 - ((n+0.5)*π/h)^2)),Line(),u)
overcoshγ_residue(u,n,k,h) = ((-1)^n)*(n+0.5)*π/(coshγ_root(u,n,k,h)*h^2)

# struct ADuctKernel{T,S} <: WienerHopfKernel
#     γ::GammaKernel{T}
#     h₊::S
#     h₋::S
# end
# ADuctKernel(k::Number,h₊::Number,h₋) = ADuctKernel(GammaKernel(k),h₊,h₋)

# GammaKernel(K::ADuctKernel) = K.γ
# wavenumber(K::ADuctKernel) = wavenumber(GammaKernel(K))
# halfheight(K::ADuctKernel,u) = u ? K.h₊ : K.h₋
# halfheights(K::ADuctKernel) =  halfheight(K,true),  halfheight(K,false)
# totalheight(K::ADuctKernel) = sum(halfheights(K)...)

# function (K::ADuctKernel)(z)
#     γ = GammaKernel(z)
#     h₊,h₋ = halfheight(K)
#    return -γ*inv(coth(h₋*γ)+coth(h₊*γ))
# end

# function checkRatioNotRational(K::ADuctKernel,n)
#     # Note we still need to treat the case h₋/h₊ rational separately, since if h₋/h₊=p/q coprime then the qth pole is counted twice in the current form.
#     # n odd corresponds to poles arising from sinh(γh₋) and n even corresponds to poles from sinh(γh₊)
#     h₊,h₋ = halfheight(K)
#     m=min(denominator(convert(Rational,h₋/h₊)),numerator(convert(Rational,h₋/h₊)))
#     @assert m/2>=n "Please ensure h₊/h₋ is not rational."
#     # A very basic fix that makes sure we dont count poles multiple times for rational h₋/h₊. Should be improved
# end

# function pole(K::ADuctKernel,u::Bool,n::Int)
#     k = wavenumber(K)
#     n +=1
#     checkRatioNotRational(K,n)
#     s = sqrt(k^2 - (n*π/totalheight(K))^2)
#     return forceabove(s,Line(),u)
# end

# function pole_residue(K::ADuctKernel,u::Bool,n::Int)
#     @assert n >= 0
#     n += 1
#     γ = -n*iπ/totalheight(K)
#     h₊,h₋ = halfheights(K)
#     return -γ*sinh(h₋*γ)*sinh(h₊*γ)*oversinhγ_residue(u,n,k,totalheight(K))
# end

# function root(K::ADuctKernel,u::Bool,n::Int)
#     checkRatioNotRational(K,n)
#     h₊,h₋ = halfheight(K)
#     s = isodd(n) ? sqrt(k^2 - ((n+1)*π/(2h₋))^2) : sqrt(k^2 - (n*π/(2h₊))^2)
#     return forceabove(s,Line(),u)
# end

# function root_residue(K::ADuctKernel,u::Bool,n::Int)
#     @assert n >= 0
#     x = root(K,u,n)
#     h₊,h₋ = halfheight(K)
#     n == 0 && return -(h₋+h₊)/(2*h₋*h₊*x)
#     isodd(n) && return -inv(h₋*x)*(-1)^((n+1)/2)*sinh((h₋+h₊)/h₋*iπ*(n+1)/2)/sinh(h₊/h₋*iπ*(n+1)/2)
#     return -inv(h₊*x)*(-1)^(n/2)*sinh((h₋+h₊)/h₊*iπ*n/2)/sinh(h₋/h₊*iπ*n/2)
# end


# # function fracsinh_der(z,k,h₋,h₊)
# # γ = γK(z,k)
# #     (z/γ)*(h₋*csch(h₋*γ)^2+h₊*csch(h₊*γ)^2)/(coth(h₋*γ)+coth(h₊*γ))^2
# # end
