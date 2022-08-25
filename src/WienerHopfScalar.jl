module WienerHopfScalar

using Reexport

@reexport using ApproxFun
using HolomorphicFun
using Lazy
using MacroTools
using PolynomialRoots
@reexport using SingularIntegralEquations
using StaticArrays
using StaticUnivariatePolynomials
using StructArrays
using UnPack

import ApproxFun: domain, space, ncoefficients, coefficients, setcanonicaldomain, evaluate
import Base: +, -, *, /, inv, literal_pow, getindex, angle, exponent, ==, convert
import HolomorphicFun: poles, roots
import MacroTools: postwalk, @capture
import SpecialFunctions: gamma

export WienerHopfKernel
export @wienerhopf
export @nopoles
export @noroots
export WienerHopfPair

export isolate_inf, isolate_poleroot
export logfactorise, factorise, factors

export GammaKernel
export NobleKernel
export RobinKernel, NormalisedRobinKernel
export RawlinsKernel
export PoroelasticK, PoroelasticI

export wavenumber, compliance

abstract type WienerHopfKernel <: AbstractScalarFunction end

include("macros/wienerhopf.jl")
include("macros/nopoles.jl")
include("util/util.jl")

include("kernel-basic/nan.jl")
include("kernel-basic/constant.jl")
include("kernel-basic/rational.jl")


include("index.jl")
include("inf.jl")
include("log.jl")


defaultspace(::WienerHopfKernel) = Chebyshev(Line{-1 / 4}(0.0))
defaultscale(K::WienerHopfKernel) = 1
defaultpoint(K::WienerHopfKernel, u::Bool) = 10im * defaultscale(K) * lsign(u) + randn(Float64) - 10

function factorise(K::WienerHopfKernel,sp::Space)
    L = isolate_inf(K)
    R = isolate_poleroot(K, Line())
    return (logfactorise(K / (L * R), sp) * L) * R
end
factorise(K::WienerHopfKernel) = factorise(K, defaultspace(K))

include("generic/template.jl")
include("kernel-zoo/zoo.jl")

end
