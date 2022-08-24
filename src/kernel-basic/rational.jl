factorise(K::ScalarRational, sp::Space) = factorise(K, domain(sp))
function factorise(K::ScalarRational, d::Domain)
    p₊, p₋ = factorise(poles(K), d)
    r₊, r₋ = factorise(roots(K), d)
    K₊ = ScalarRational(p₋, r₋)
    K₋ = ScalarRational(p₊, r₊)
    return SplitRationalKernel(K₊, K₋)
end

struct SplitRationalKernel{T} <: WienerHopfKernel
    plus::ScalarRational{T}
    minus::ScalarRational{T}
end

evaluate(K::SplitRationalKernel, z) = evaluate(K, z, true) * evaluate(K, z, false)

(K::SplitRationalKernel)(u::Bool) = u ? K.plus : K.minus

function evaluate(K::SplitRationalKernel, z, u)
    if u
        return K.plus(z)
    else
        return K.minus(z)
    end
end

poles(K::SplitRationalKernel, u) = poles(K(u))
roots(K::SplitRationalKernel, u) = roots(K(u))

function isolate_poleroot(K::AbstractScalarFunction, d::Union{Space,Domain}, z₋=-10.0im, z₊=10.0im)
    L = factorise(ScalarRational(poles(K), roots(K)), d)
    return isolate_poleroot(L, z₋, z₊)
end
function isolate_poleroot(K::SplitRationalKernel, z₋, z₊)
    K₋ = isolate_poleroot(K.minus, z₊)
    K₊ = isolate_poleroot(K.plus, z₋)
    return SplitRationalKernel(K₊, K₋)
end
function isolate_poleroot(K::AbstractScalarFunction, z)
    return isolate_poleroot(poles(K), roots(K), z)
end
function isolate_poleroot(poles::AbstractVector{T}, roots::AbstractVector{T}, z) where {T}
    p = copy(poles)
    r = copy(roots)
    n = length(p) - length(r)
    n < 0 && append!(p, z * ones(T, -n))
    n > 0 && append!(r, z * ones(T, n))
    return ScalarRational(p, r)
end
