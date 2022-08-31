
function evaluate(K::T, z) where {T<:WienerHopfKernel}
    return error("Implement evaluate(K,z) for K::$T")
end

function evaluate(K::T, z, u) where {T<:WienerHopfKernel}
    return  error("Factorisation not known for type $T. Try factorise(K,sp) or implement explicitly.")
end

function leftlimit(K::T) where {T<:WienerHopfKernel}
    return error("Define leftlimit(K) = K(-Inf)")
end

function rightlimit(K::T) where {T<:WienerHopfKernel}
    return error("Define rightlimit(K) = K(+Inf)")
end

defaultspace(::WienerHopfKernel) = Chebyshev(Line{-1 / 4}(0.0))
defaultscale(::WienerHopfKernel) = 1.0
defaultpoint(u::Bool) = 10 *(1im * lsign(u) + u)
defaultpoint(K::WienerHopfKernel, u::Bool) = defaultscale(K) * defaultpoint(u)
