# roots
function filter_roots!(z::AbstractArray{T}, f::Function, atol=100 * eps(T)) where T
    return filter!(z -> abs(f(z)) < atol, z)
end

function candidate_roots(K::WienerHopfKernel)
    p = roots_poly(K)
    coeffs = collect(p.coeffs)
    return PolynomialRoots.roots(coeffs)
end

function roots(K::WienerHopfKernel, atol=100 * eps())
    return filter_roots!(candidate_roots(K), K, atol)
end

# poles
function filter_poles!(z::AbstractArray{T}, f::Function, atol=100 * eps(T)) where T
    return filter!(z -> abs(f(z)) > inv(atol) || isnan(abs(f(z))), z)
end

function candidate_poles(K::WienerHopfKernel)
    p = poles_poly(K)
    coeffs = collect(p.coeffs)
    return PolynomialRoots.roots(coeffs)
end

function poles(K::WienerHopfKernel, atol=100 * eps())
    return filter_poles!(candidate_poles(K), K, atol)
end

import Base.:^

# !!! Type piracy
function ^(x::Polynomial{N,T}, n::Int) where {N,T}
    nmax = 20
    n == 0 && return Polynomial{1,T}(NTuple{1,T}(one(T)))
    n == 1 && return x
    n == 2 && return x * x
    iseven(n) && return (x^2)^(n ÷ 2)
    n < nmax && return x * (x^(n - 1))
    return @error "implement ^($(typeof(x)),n) for n > $nmax using optimal multiplication chain"
end
