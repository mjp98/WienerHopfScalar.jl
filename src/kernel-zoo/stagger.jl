struct StaggerKernel{S<:Bool,R<:Complex,T<:Real} <: WienerHopfKernel
    γ::GammaKernel{R}
    h::T
end
GammaKernel(K::StaggerKernel) = K.γ
height(K::StaggerKernel) = K.h

evaluate(K::StaggerKernel{true},z) = 1 + exp(height(K)*GammaKernel(K,z))
evaluate(K::StaggerKernel{false},z) = 1 - exp(height(K)*GammaKernel(K,z))

#roots(K::StaggerKernel{true}) =
poles(K::StaggerKernel{S,R}) where {S,R} = R[]


# Use rational space for factorisation by default
