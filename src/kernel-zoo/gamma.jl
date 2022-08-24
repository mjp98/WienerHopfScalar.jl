"""
    GammaKernel(k)

Wiener--Hopf gamma
"""
struct GammaKernel{T} <: WienerHopfKernel
    k::T

    function GammaKernel(k::T) where {T}
        return new{complex(float(T))}(complex(float(k)))
    end
end
wavenumber(γ::GammaKernel) = γ.k

function evaluate(K::GammaKernel, z, u)
    if u
        return sqrt(-im * (z + wavenumber(K)))
    else
        return sqrt(im * (z - wavenumber(K)))
    end
end
evaluate(K::GammaKernel, z) = K(z, true) * K(z, false)

factorise(K::GammaKernel, args...) = K

function literal_pow(::typeof(^), K::GammaKernel{T}, ::Val{n}) where {T,n}
    n < 0 && return ScalarReciprocal(literal_pow(^, K, -n))
    n == 0 && return ScalarConstant(one(T))
    n == 1 && return K
    n == 2 && return Polynomial((-wavenumber(K)^2, zero(T), one(T)))
    return z -> K(z)^n
end

getindex(K::GammaKernel,u::Bool) = u ? Sqrt{true}(-wavenumber(K),-0.5,-1/sqrt(im)) : Sqrt{true}(wavenumber(K),0.5,sqrt(im))

function GammaKernel(K::WienerHopfKernel, z)
    γ = GammaKernel(K)
    return γ(z)
end

defaultscale(K::GammaKernel) = wavenumber(K)
