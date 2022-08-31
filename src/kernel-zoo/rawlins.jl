struct RawlinsKernel{T} <: WienerHopfKernel
    γ::GammaKernel{T}
    μ::T
    δ::T
end

function RawlinsKernel(k::T,μ::S,δ::V) where {T<:Number,S,V}
    k,μ,δ = promote(complex(float(k)),μ,δ)
    return RawlinsKernel(GammaKernel(k),μ,δ)
end
GammaKernel(K::RawlinsKernel) = K.γ
GammaKernel(K::RawlinsKernel,z) = K.γ(z)
compliance(K::RawlinsKernel) = K.μ
pressurefree_wavenumber(K::RawlinsKernel) = K.δ
wavenumber(K::RawlinsKernel) = wavenumber(GammaKernel(K))

function evaluate(K::RawlinsKernel,α)
    return 1 + compliance(K)*(α-pressurefree_wavenumber(K))/GammaKernel(K,α)
end

function roots_poly(K::RawlinsKernel)
    @unpack μ,γ,δ = K
    μ² = μ^2
    u² = @SVector [δ^2,-2,1]
    γ² = γ^2
    return Polynomial((u²*μ²)...) - γ²
end

poles(::RawlinsKernel{T},args...) where T = T[]

leftlimit(K::RawlinsKernel) = 1-compliance(K)
rightlimit(K::RawlinsKernel) = 1+compliance(K)
