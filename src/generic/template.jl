
function evaluate(K::T, z) where {T<:WienerHopfKernel}
    return error("Implement evaluate(K,z) for K::$T")
end
function evaluate(K::T, z, u) where {T<:WienerHopfKernel}
    return  error("Factorisation not known for type $T. Try factorise(K,sp) or implement explicitly.")
end
# function isolate_poleroot(K::T, sp) where {T<:WienerHopfKernel}
#     return @error "Implement isolate_poleroot(K,sp) for K::$T"
# end
# function isolate_inf(K::T) where {T<:WienerHopfKernel}
#     return @error "Implement isolate_inf(K,sp) for K::$T"
# end

leftlimit(K::T) where {T<:WienerHopfKernel} = error("Define leftlimit(K) = K(-Inf)")
rightlimit(K::T) where {T<:WienerHopfKernel} = error("Define rightlimit(K) = K(+Inf)")



# # Generic factorization

# factorise(K::WienerHopfKernel) = factorise(K, defaultspace(K))

# function factorise(K::WienerHopfKernel, sp, args...)
#     L = isolate_inf(K)
#     R = isolate_poleroot(K, Line())
#     return logfactorise(K / (L * R), sp, args...) * L * R
# end

# # Defaults

# defaultspace(::WienerHopfKernel) = Chebyshev(Line{-1 / 4}(0.0))
# defaultscale(K::WienerHopfKernel) = 1
# defaultpoint(K::WienerHopfKernel, u::Bool) = 10im * defaultscale(K) * lsign(u) + π - 10
