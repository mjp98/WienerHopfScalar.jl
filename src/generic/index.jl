function isolate_index(
    n::Int,
    innerpoint=defaultpoint(true),
    outerpoint=defaultpoint(false)
)

    n > 0 && return ScalarRational(innerpoint * ones(n), outerpoint * ones(n))
    n < 0 && return ScalarRational(outerpoint * ones(-n), innerpoint * ones(-n))

    return factorise(ScalarConstant(1))
end

function logbranchdata(f::Fun; ε=1e-6)
    F = ApproxFun.setcanonicaldomain(copy(f))
    Fr, Fi = reim(F)
    x = ApproxFun.roots(Fi)
    filter!(z -> abs(z) < 0.9999, x)  # remove endpoints
    filter!(z -> Fr(z) < 0, x)         # log branch on negative real axis
    sort!(x)                           # sort

    any(abs(z) < ε for z in x) && @warn "possible branch cut at origin"

    n = zeros(Int, length(x) + 1)
    for (i, z) in enumerate(x)
        n[i+1] = n[i] + lsign(Fi(leftvalue(z, ε)) < 0)
    end
    x = [-Inf; x]
    return BranchData(x, n)
end

struct BranchData{T}
    jumps::StructVector{BranchJump{T}}
end

BranchData(a, b) = BranchData(StructArray{BranchJump{eltype(a)}}((a, b)))

jumps(d::BranchData) = d.jumps
windingnumber(d::BranchData) = last(jumps(d).n)

struct BranchJump{T<:Number}
    x::T
    n::Int
end

struct ContinuousLog{T<:Fun,S<:BranchData}
    f::T
    data::S
end
@forward ContinuousLog.data jumps

function evaluate(f::ContinuousLog, z)
    fz = f.f(z)
    for (; x, n) in jumps(f)
        z < x && return log(fz) - n * i2π
    end
    return log(fz)
end

(f::ContinuousLog)(z) = evaluate(f, z)
