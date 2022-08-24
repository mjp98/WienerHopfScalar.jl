struct PoroelasticK{T,S} <: WienerHopfKernel
    Ω::S
    R::S# = 0.026
    αₕ::S# = 0.0014
    ε::S# = 0.002
    Kᵣ::S# = 2R
    μ::S# = αₕ*Kᵣ/(π*R^2)
    ψ₀::S# = π/2
    ν::S # = 0.3
    k::T ## = complexfloat(Ω*sin(ψ₀))
    p²::SVector{2,T} # = SVector(((k^2 - Ω^2 - 1), (k^2 - Ω^2 + 1)))

    function PoroelasticK(;
        Ω=0.014,    # 0.1242 < Ω < 0.3104,
        R=0.026,
        αₕ=0.0014,
        ε=0.002,
        Kᵣ=2R,
        μ=αₕ * Kᵣ / (π * R^2),
        ψ₀=π / 2,
        ν=0.3,
        k=complexfloat(Ω * sin(ψ₀)))
        x = complexfloat(k)
        y = promote(Ω, R, αₕ, ε, Kᵣ, μ, ψ₀, ν)
        p² = SVector(((k^2 - Ω^2 - 1), (k^2 - Ω^2 + 1)))
        new{eltype(x),eltype(y)}(y..., x, p²)
    end
end

wavenumber(K::PoroelasticK) = K.k
GammaKernel(K::PoroelasticK) = GammaKernel(wavenumber(K))

poles²(K::PoroelasticK) = K.p²
poles⁺(K::PoroelasticK) = forceabove(sqrt.(poles²(K)), Line())

function poles_poly(K::PoroelasticK)
    α₁², α₂² = K.p²
    return Polynomial((α₁² * α₂², 0, -(α₁² + α₂²), 0, 1))
end

function roots_poly(K::PoroelasticK)
    @unpack k, αₕ, Kᵣ, Ω, ε, μ = K
    ζ = Ω * (1 + αₕ * Kᵣ) / ε
    γ = GammaKernel(K)
    x = poles_poly(K)
    return (γ^2 * x^2) - (μ^2 * x^2 + 2 * μ * ζ * x + ζ^2)
end

function evaluate(K::PoroelasticK, z)
    @unpack k, αₕ, Kᵣ, Ω, ε, μ = K
    z² = z^2
    α₁², α₂² = poles²(K)
    γ = sqrt(-im * (z + k)) * sqrt(im * (z - k))
    ζ = Ω * (1 + αₕ * Kᵣ) / ε
    return γ + μ + ζ / ((z² - α₁²) * (z² - α₂²))
end

isolate_inf(K::PoroelasticK) = false, GammaKernel(K)

struct PoroelasticI{T,S} <: WienerHopfKernel
    K::PoroelasticK{T,S}
end
PoroelasticK(I::PoroelasticI) = I.K
wavenumber(I::PoroelasticI) = wavenumber(PoroelasticK(I))
GammaKernel(I::PoroelasticI) = GammaKernel(PoroelasticK(I))
roots_poly(I::PoroelasticI) = roots_poly(PoroelasticK(I))
poles(::PoroelasticI{T,S}) where {T,S} = T[]

function evaluate(I::PoroelasticI, z)
    @unpack k, αₕ, Kᵣ, Ω, ε, μ = PoroelasticK(I)
    z² = z^2
    α₁², α₂² = poles²(PoroelasticK(I))
    γ = sqrt(-im * (z + k)) * sqrt(im * (z - k))
    ζ = Ω * (1 + αₕ * Kᵣ) / ε
    return ((z² - α₁²) * (z² - α₂²)) * (γ + μ) + ζ
end

isolate_inf(I::PoroelasticI) = false, GammaKernel(I) * WienerHopfPair(Polynomial((10im, 1.0))^2, Polynomial((-10im, 1.0))^2)
