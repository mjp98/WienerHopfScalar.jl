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
wavenumber(K::RawlinsKernel) = wavenumber(GammaKernel(K))

evaluate(K::RawlinsKernel,α) = 1 + K.μ*(α-K.δ)/K.γ(α)

function roots_poly(K::RawlinsKernel)
    @unpack μ,γ,δ = K
    μ² = μ^2
    u² = @SVector [δ^2,-2,1]
    γ² = γ^2
    return Polynomial((u²*μ²)...) - γ²
end

poles(::RawlinsKernel{T},args...) where T = T[]

leftlimit(K::RawlinsKernel) = 1-K.μ
rightlimit(K::RawlinsKernel) = 1+K.μ
