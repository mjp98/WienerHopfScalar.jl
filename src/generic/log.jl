# @forward LogFun.f (isabove, coefficients, ncoefficients, space, domain)
# @forward SplitNormalised.logf (isabove, coefficients, ncoefficients, space, domain, expcauchy)

expcauchy(f::Fun, z) = exp(complex(cauchy(f, z)))

abstract type AbstractLogSplit end
struct LogFun{T<:Fun} <: AbstractLogSplit
    f::T
end
expcauchy(L::LogFun, z) = expcauchy(L.f, z)
isabove(L::LogFun, z) = isabove(L.f, z)

LogFun(f::Function, sp...) = LogFun(Fun(log ∘ f, sp...))

struct SplitNormalised{T,S<:AbstractLogSplit} <: WienerHopfKernel
    f::T
    logf::S
end

evaluate(K::SplitNormalised, z) = K.f(z)
logf(K::SplitNormalised) = K.logf
expcauchy(K::SplitNormalised, z) = expcauchy(logf(K), z)
isabove(K::SplitNormalised, z) = isabove(logf(K), z)

function evaluate(f::SplitNormalised, z, u)
    above = isabove(f, z)
    (above && u) && return expcauchy(f, z)
    (above && !u) && return f(z) / expcauchy(f, z)
    (!above && u) && return f(z) * expcauchy(f, z)
    (!above && !u) && return inv(expcauchy(f, z))
end

for op in [:coefficients, :ncoefficients, :space, :domain]
    @eval $op(L::LogFun) = $op(L.f)
    @eval $op(L::SplitNormalised) = $op(logf(L))
end

function logfactorise(
    f,
    sp...;
    innerpoint=defaultpoint(true),
    outerpoint=defaultpoint(false)
)

    F = ApproxFun.setcanonicaldomain(Fun(f, sp...))

    branchdata = logbranchdata(F)

    winding_number = windingnumber(branchdata)

    @info "winding_number(K,Γ) = $winding_number"

    if !iszero(winding_number)
        @info "winding_number non-zero: isolating..."
        w = isolate_index(winding_number, innerpoint, outerpoint)
        return logfactorise(f / w, sp...) * factorise(w, first(sp))
    end

    if any(n != 0 for n in jumps(branchdata).n)
        continuouslogf = ContinuousLog(F, branchdata)
        c = Fun(continuouslogf, space(F), ncoefficients(F))
        L = LogFun(Fun(first(sp), coefficients(c)))
        return SplitNormalised(f, L)
    end

    return SplitNormalised(f, LogFun(f, sp...))
end
