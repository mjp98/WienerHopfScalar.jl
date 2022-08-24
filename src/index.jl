function isolate_index(n::Int, z₊=10.0im, z₋=-10.0im)
    # can we make this an already factorised object?
    n > 0 && return ScalarRational(z₊ * ones(n), z₋ * ones(n))
    n < 0 && return ScalarRational(z₋ * ones(-n), z₊ * ones(-n))
    return factorise(ScalarConstant(1))
end

# Refactor this into smaller functions

struct BranchJump{T<:Number}
    x::T
    n::Int
end

struct BranchData{T}
    jumps::StructVector{BranchJump{T}}
end

BranchData(a, b) = StructArray{BranchJump{eltype(a)}}((a, b))

jumps(d::BranchData) = d.jumps
windingnumber(d::BranchData) = sum(jumps(d).n)

function logbranchdata(f::Fun; tol=1e-6)
    if domain(f) !== canonicaldomain(f)
        canonicalf = setcanonicaldomain(copy(f))
        return logbranchdata(canonicalf; tol, ε)
    end
    realF, imagF = reim(F)
    x = ApproxFun.roots(imagF)
    filter!(z -> abs(z) < 1 - sqrt(tol), x)  # remove endpoints
    filter!(z -> realF(z) < 0, x)            # log branch on negative real axis

    any(abs(z) < tol for z in x) && @warn "possible branch cut at origin"

    sort!(x)                                 # sort
    Δ = [lsign(imagF(leftvalue(z, tol)) < 0) for z in x]
    return BranchData(x, cumsum(Δ))
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

# function continuouslog(logf::Fun,data::BranchData)
#     if any(n != 0 for n in data.jumps.n)
#         @warn "Branch cut crossed, fixing..."
#         F = ApproxFun.setcanonicaldomain(Fun(f, sp...))
#         function logf(z)
#             for (x, n) in zip(xbranch, nbranch)
#                 z < x && return log(F(z)) - n * i2π
#             end
#             return log(F(z))
#         end
#         c = Fun(logf, space(F), ncoefficients(F))
#     end

# end



function windingdata(f, sp...; ε=1e-6, tol = 1e-6)
    F = ApproxFun.setcanonicaldomain(copy(Fun(f, sp...)))
    Fr, Fi = reim(F)
    xbranch = sort(ApproxFun.roots(Fi))
    # remove endpoints.
    filter!(z -> abs(z) < 0.9999, xbranch)

    # log branch on negative real axis
    filter!(z -> Fr(z) < 0, xbranch)
    # jump in branch across cut


   # x = logbranchjumps(f; tol)

    nbranch = zeros(Int, length(xbranch) + 1)
    for (i, x) in enumerate(xbranch)
        abs(x) < ε && @warn "possible branch cut at origin"
        # Δ = log(f(rightvalue(x,ε))) - log(f(leftvalue(x,ε)))
        #nbranch[i+1] = nbranch[i] + Int(round(Δ/(2π*im)))
        nbranch[i+1] = nbranch[i] + lsign(Fi(leftvalue(x, ε)) < 0)
    end
    windingnumber = last(nbranch)
    return windingnumber, xbranch, nbranch, F
end
