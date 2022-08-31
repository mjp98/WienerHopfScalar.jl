module WienerHopfScalar


using Reexport

@reexport using ApproxFun
@reexport using SingularIntegralEquations

using HolomorphicFun
using Lazy
using MacroTools
using PolynomialRoots
using StaticArrays
using StaticUnivariatePolynomials
using StructArrays
using UnPack

import ApproxFun: domain, space, ncoefficients, coefficients, setcanonicaldomain, evaluate, factors
import Base: +, -, *, /, inv, literal_pow, getindex, angle, exponent, ==, convert
import HolomorphicFun: poles, roots
import MacroTools: postwalk, @capture
import SpecialFunctions: gamma

export WienerHopfKernel
export @wienerhopf
export @nopoles
export @noroots
export WienerHopfPair
export WienerHopfNaN

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

include("generic/generic.jl")

include("kernel-zoo/zoo.jl")

end
