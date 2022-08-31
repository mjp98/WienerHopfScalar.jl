function isolate_inf(
    K::WienerHopfKernel;
    innerpoint = defaultpoint(K, true),
    outerpoint = defaultpoint(K,false)
    )

    return isolate_inf(leftlimit(K), rightlimit(K); innerpoint, outerpoint)
end

function isolate_inf(
    l,
    r;
    innerpoint=defaultpoint(true),
    outerpoint=defaultpoint(false),
    tol=1e-14
    )

    check_asymptotes(l, r; tol)

    return split_asymptotes(l, r, innerpoint, outerpoint)
end

function check_asymptotes(l, r; tol=1e-14)
    abs(l) < tol && return @error "K ≈ 0 as z → -∞"
    abs(r) < tol && return @error "K ≈ 0 as z → +∞"
    abs(l) > inv(tol) && return @error "abs(K) >> 1 as z → -∞"
    abs(r) > inv(tol) && return @error "abs(K) >> 1 as z → +∞"
end

function split_asymptotes(
    l,
    r,
    innerpoint=defaultpoint(true),
    outerpoint=defaultpoint(false)
    )
    return SplitAsymptotes(cphase(l),cphase(r),innerpoint,outerpoint)
end


struct SplitAsymptotes{T<:Number} <: WienerHopfKernel
    lphase::T
    rphase::T
    innerpoint::T
    outerpoint::T
    function SplitAsymptotes(a,b,c,d)
        x = promote(a,b,c,d)
        new{eltype(x)}(x...)
    end
end

innerpoint(K::SplitAsymptotes) = K.innerpoint
outerpoint(K::SplitAsymptotes) = K.outerpoint
lphase(K::SplitAsymptotes) = K.lphase
rphase(K::SplitAsymptotes) = K.rphase

function half2zero(f::SplitAsymptotes, z::Number, u::Bool)
    if u
        return  (log(-im * (z - outerpoint(f))) + iπ/2) / i2π
    else
        return -(log( im * (z - innerpoint(f))) - iπ/2) / i2π
    end
end

function one2zero(f::SplitAsymptotes, z::Number)
    return half2zero(f, z, true) + half2zero(f, z, false)
end

function evaluate(f::SplitAsymptotes, z::Number)
    theta = rphase(f) + (lphase(f) - rphase(f)) * one2zero(f, z)
    return exp(i2π * theta)
end

function evaluate(f::SplitAsymptotes, z::Number, u::Bool)
    theta = rphase(f)/2 + (lphase(f) - rphase(f)) * half2zero(f, z, u)
    return exp(i2π * theta)
end
