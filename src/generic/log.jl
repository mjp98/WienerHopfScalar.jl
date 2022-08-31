# expcauchy(f::Fun, z) =  exp(complex(cauchy(f, z)))
# abstract type AbstractLogSplit end



# """

# """
# struct LogFun{T} <: AbstractLogSplit
#     f::T
# end
# LogFun(f::Function, sp, args...) = LogFun(Fun(z->log(f(z)), sp, args...))
# expcauchy(f::LogFun, z) = expccauchy(f,f,z)
# @forward LogFun.f (isabove, coefficients, ncoefficients, space, domain)

# """

# """
# struct LogSplitKernel{T,S<:AbstractLogSplit} <: WienerHopfKernel
#     f::T
#     splitlog::S
# end
# splitlog(K::LogSplitKernel) = K.splitlog

# # function evaluate(K:::LogSplitKernel,z)
# #     return evaluate()
# @forward LogSplitKernel.splitlog (isabove, coefficients, ncoefficients, space, domain, expcauchy)

# function evaluate(f::LogSplitKernel, z, u)
#     above = isabove(f, z)
#     (above && u) && return expcauchy(f, z)
#     (above && !u) && return f(z) / expcauchy(f, z)
#     (!above && u) && return f(z) * expcauchy(f, z)
#     (!above && !u) && return inv(expcauchy(f, z))
# end


# function logfactorise(f, sp, args...; z₋=-10.0im, z₊=10.0im + 10)
#     windingnumber, xbranch, nbranch, F = windingdata(f, sp,args...)
#     @info "winding_number(K,Γ) = $windingnumber"
#     if windingnumber != 0
#         @info "isolating winding_number"
#         w = isolate_index(windingnumber, z₊, z₋)
#         return logfactorise(f / w, sp, args...) * factorise(w, sp)
#     end
#     if any(n != 0 for n in nbranch)
#         @warn "Branch cut crossed, fixing..."
#         F = ApproxFun.setcanonicaldomain(Fun(f, sp,args...))
#         function logf(z)
#             for i in zip(xbranch, nbranch)
#                 x, n = i
#                 z < x && return log(F(z)) - n * i2π
#             end
#             return log(F(z))
#         end
#         c = Fun(logf, space(F), ncoefficients(F))
#         L = LogFun(Fun(sp, coefficients(c)))
#         return LogSplitKernel(f, L)
#     end
#     return LogSplitKernel(f, LogFun(f, sp))
# end



# function continuouslog(f,sp...)


# end


abstract type AbstractLogSplit end
struct LogFun{T<:Fun} <: AbstractLogSplit
    f::T
end
expcauchy(L::LogFun, z) = exp(complex(cauchy(L.f, z)))
isabove(L::LogFun, z) = isabove(L.f, z)

LogFun(f::Function, sp...) = LogFun(Fun(log ∘ f, sp...))

struct LogSplitKernel{T,S<:AbstractLogSplit} <: WienerHopfKernel
    f::T
    splitlog::S
end

evaluate(K::LogSplitKernel, z) = K.f(z)
splitlog(K::LogSplitKernel) = K.splitlog
expcauchy(K::LogSplitKernel, z) = expcauchy(splitlog(K), z)
isabove(K::LogSplitKernel, z) = isabove(splitlog(K), z)

function evaluate(f::LogSplitKernel, z, u)
    above = isabove(f, z)
    (above && u) && return expcauchy(f, z)
    (above && !u) && return f(z) / expcauchy(f, z)
    (!above && u) && return f(z) * expcauchy(f, z)
    (!above && !u) && return inv(expcauchy(f, z))
end

for op in [:coefficients, :ncoefficients, :space, :domain]
    @eval $op(L::LogFun) = $op(L.f)
    @eval $op(L::LogSplitKernel) = $op(splitlog(L))
end

function logfactorise(
    f,
    sp...;
    innerpoint=defaultpoint(true),
    outerpoint=defaultpoint(false)
    )

    windingnumber, xbranch, nbranch, F = windingdata(f, sp...)
    @info "winding_number(K,Γ) = $windingnumber"
    if windingnumber != 0
        @info "isolating winding_number"
        w = isolate_index(windingnumber, innerpoint, outerpoint)
        return logfactorise(f / w, sp...) * factorise(w, sp[1])
    end
    if any(n != 0 for n in nbranch)
        @info "Branch cut crossed, fixing..."
        F = ApproxFun.setcanonicaldomain(Fun(f, sp...))
        function logf(z)
            for i in zip(xbranch, nbranch)
                x, n = i
                z < x && return log(F(z)) - n * i2π
            end
            return log(F(z))
        end
        c = Fun(logf, space(F), ncoefficients(F))
        L = LogFun(Fun(first(sp), coefficients(c)))
        return LogSplitKernel(f, L)
    end
    return LogSplitKernel(f, LogFun(f, sp...))
end
