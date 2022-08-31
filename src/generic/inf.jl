isolate_inf(K::WienerHopfKernel; kwargs...) = isolate_inf(leftlimit(K), rightlimit(K); kwargs...)

function parse_inf(l, r; tol=1e-14)
    abs(l) < tol && return @error "K ≈ 0 as z → -∞"
    abs(r) < tol && return @error "K ≈ 0 as z → +∞"
end

struct IsolatedInfLog{T<:Number} <: WienerHopfKernel
    c₊::T
    c₋::T
    z₊::T
    z₋::T
end

cphase(K::IsolatedInfLog) = K.c₋, K.c₊

function half2zero_phase(f::IsolatedInfLog, z, u)
    if u
        return -(log(-im * (z - f.z₋)) - log(-1im)) / (-i2π)
    else
        return  (log( im * (z - f.z₊)) - log( 1im)) / (-i2π)
    end
end

function one2zero_phase(f::IsolatedInfLog, z)
    return half2zero_phase(f, z, true) + half2zero_phase(f, z, false)
end

function evaluate(f::IsolatedInfLog, z)
    c₋, c₊ = cphase(f)
    x = c₊ + (c₋ - c₊) * one2zero_phase(f, z)
    return exp(i2π * x)
end

function evaluate(f::IsolatedInfLog, z, u)
    c₋, c₊ = cphase(f)
    x = c₊ / 2 + (c₋ - c₊) * half2zero_phase(f, z,  u)
    return exp(i2π * x)
end

function isolate_inf(l, r; z₊=10im, z₋=-10im, tol=1e-14)
    parse_inf(l, r; tol)
    isolated = make_inf(l, r, z₊, z₋)
    return isolated
end

function make_inf(l, r, z₊=10im, z₋=-10im)
    c₋, c₊ = cphase(l), cphase(r)
    return IsolatedInfLog(promote(c₊, c₋, z₊, z₋)...)
end

# function isolate_infpow(l, r; z₊=10im, z₋=-10im, tol=1e-14)
#     parse_inf(l, r; tol)
#     isolated = make_infpow(l, r, z₊, z₋)
#     return isolated
# end

# function make_infpow(l, r, z₊=10im, z₋=-10im)
#     c₋, c₊ = cphase(l), cphase(r)
#     scale = cispi(c₊)
#     Δ = c₊ - c₋
#     K₊ = ScalarPower(scale, z₋,  Δ,  im)
#     K₋ = ScalarPower(scale, z₊, -Δ, -im)
#     return WienerHopfPair(K₊, K₋)
# end

# isolate_infpow(K::WienerHopfKernel; kwargs...) = isolate_infpow(leftlimit(K), rightlimit(K); kwargs...)
