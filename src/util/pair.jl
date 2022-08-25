struct WienerHopfPair{P,M} <: WienerHopfKernel
    p::P
    m::M
end
getindex(K::WienerHopfPair, u::Bool) = u ? K.p : K.m
#evaluate(K::WienerHopfPair, u::Bool) = K[u]
evaluate(K::WienerHopfPair, α, u) = K[u](α)
evaluate(K::WienerHopfPair, α) = evaluate(K, α, true) * evaluate(K, α, false)
factorise(K::WienerHopfPair, args...) = K
factors(K::WienerHopfKernel) = K[true], K[false]
# *(a::WienerHopfPair,b::WienerHopfPair) = WienerHopfPair(a[true]*b[true],a[false]*b[false])
