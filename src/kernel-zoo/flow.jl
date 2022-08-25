
# using ApproxFun
# import Base: *, /, -

# mysqrt(z, k) = sqrt(complex(z - k))

# direction(θ::Real) = exp(im * pi * (1 - θ))
# direction2(θ::Real) = exp(im * pi * (1 - θ) / 2)

# function mysqrt(z, k, θ)
#     d = direction(θ)
#     d2 = direction2(θ)
#     return sqrt(d * (z - k)) / d2
# end

# # -----------------------------
# # Sqrt - shifted with branch cut at angle θπ
# # -----------------------------

# struct Sqrt{R,T,S,U} <: Function
#     k::T
#     θ::S
#     c::U
#     function Sqrt{R}(k, θ, c=1) where {R}
#         w = complex(float(k))
#         new{R,typeof(w),typeof(θ),typeof(c)}(w, θ, c)
#     end
# end
# center(K::Sqrt) = K.k
# scale(K::Sqrt) = K.c
# angle(K::Sqrt) = K.θ


# direction2(f::Sqrt) = direction2(angle(f))

# (f::Sqrt{true})(z) = scale(f) * mysqrt(z, center(f), angle(f))
# (f::Sqrt{false})(z) = scale(f) / mysqrt(z, center(f), angle(f))

# *(a::Number, b::Sqrt{T}) where {T} = Sqrt{T}(center(b), angle(b), a * scale(b))
# /(a::Number, b::Sqrt{T}) where {T} = Sqrt{!T}(center(b), angle(b), a / scale(b))

# *(b::Sqrt{T}, a::Number) where {T} = Sqrt{T}(center(b), angle(b), a * scale(b))
# /(b::Sqrt{T}, a::Number) where {T} = Sqrt{T}(center(b), angle(b), scale(b) / a)

# -(b::Sqrt{T}) where {T} = Sqrt{T}(center(b), angle(b), -scale(b))

# cutdomain(f::Sqrt) = Ray{angle(f)}(center(f))
# cutspace(f::Sqrt) = JacobiWeight(-0.5, -0.5, cutdomain(f))

# Δ(f::Sqrt{false}) = Fun(cutspace(f), [2im, -2im]) * direction2(f) * scale(f)
# Δ(f::Sqrt{true}) = Fun(cutspace(f), [-2im, -2im]) / direction2(f) * scale(f)

# Σ(f::Sqrt{true}) = Fun(cutdomain(f), [0])
# Σ(f::Sqrt{false}) = Fun(cutdomain(f), [0])

# ε = 1e-12
# z = 0.2 + 1im

# S = Sqrt{true}(0.2, 0.5)
# u = Δ(S)
# u(0.2 + im)
# S(z - ε) - S(z + ε)

# roots(::Sqrt{T,S}) where {T} = S[]

# struct ZeroOverSqrt{T,S} <: Function
#     s::Sqrt{true}
#     x::S
#     c::S
# end
# (f::ZeroOverSqrt{true})(z) = scale(f) * (f.x - z) / f.s(z)
# (f::ZeroOverSqrt{false})(z) = scale(f) * f.s(z) / (f.x - z)

# *(a::Number, b::ZeroOverSqrt{T}) where {T} = ZeroOverSqrt{T}(b.s, b.x, a * scale(b))
# /(a::Number, b::ZeroOverSqrt{T}) where {T} = ZeroOverSqrt{!T}(b.s, b.x, a / scale(b))


# roots(f::ZeroOverSqrt{true,S}) = S[f.x]
# roots(::ZeroOverSqrt{false,S}) = S[]

# cutdomain(f::ZeroOverSqrt{T}) where {T} = cutdomain(f.s)

# struct FlowKernel{T}
#     γ::GammaKernel{T}
#     p::Polynomial{T}
#     function FlowKernel(k, k₀, k₁)
#         k, k₀, k₁ = promote(k, k₀, k₁)
#         new{eltype(k)}(GammaKernel(k), Polynomial((k₀, k₁)))
#     end
# end
# (K::FlowKernel)(z) = K.γ(z) / K.p(z)

# GammaKernel(K::FlowKernel) = K.γ

# factor(K::GammaKernel, u) = Sqrt{true}(lsign(u) * wavenumber(K), lsign(u) * 0.5)

# function factor(K::FlowKernel, u)
#     k₀, k₁ = coeffs(K.p)
#     γ = GammaKernel(K)
#     γu = factor(γ, u)
#     if imag(-k₀ / k₁) > 0
#         !u && return ZeroOverSqrt{false}(γu, -k₀ / k₁, -1 / k₁)
#         u && return γu
#     else
#         !u && return γu
#         u && return ZeroOverSqrt{false}(γu, -k₀ / k₁, -1 / k₁)
#     end
# end

# factorise(K::FlowKernel) = WienerHopfPair(factor(K, true), factor(K, false))
