struct RobinKernel{T} <: WienerHopfKernel
    γ::GammaKernel{T}
    μ::T
end

function RobinKernel(k::T, μ::S) where {T<:Number,S}
    k, μ = promote(complex(float(k)), μ)
    RobinKernel(GammaKernel(k), μ)
end

GammaKernel(K::RobinKernel) = K.γ
compliance(K::RobinKernel) = K.μ
wavenumber(K::RobinKernel) = wavenumber(GammaKernel(K))

evaluate(K::RobinKernel, z) = evaluate(GammaKernel(K), z) + compliance(K)

roots_poly(K::RobinKernel) = Polynomial((-compliance(K)^2 - wavenumber(K)^2, 0, 1))

isolate_inf(K::RobinKernel, args...) = false, GammaKernel(K)

struct NormalisedRobinKernel{T} <: WienerHopfKernel
    K::T
end
NormalisedRobinKernel(k::T, μ::S) where {T<:Number,S} = NormalisedRobinKernel(RobinKernel(k, μ))

RobinKernel(K::NormalisedRobinKernel) = K.K

for op in [:GammaKernel, :compliance, :wavenumber, :roots_poly]
    @eval $op(K::NormalisedRobinKernel) = $op(RobinKernel(K))
end

evaluate(K::NormalisedRobinKernel, z) = 1 + compliance(K) / (GammaKernel(K))(z)

roots_poly(K::NormalisedRobinKernel, opts...) = roots_poly(RobinKernel(K), opts...)

poles(K::RobinKernel{T}, atol=eps()) where {T} = T[]
poles(K::NormalisedRobinKernel, atol=eps()) = poles(RobinKernel(K), atol)
