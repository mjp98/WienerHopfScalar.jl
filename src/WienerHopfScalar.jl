module WienerHopfScalar

using ApproxFun
using HolomorphicFun
using MacroTools
using PolynomialRoots
using SingularIntegralEquations
using StaticArrays
using StaticUnivariatePolynomials
using UnPack

import ApproxFun: domain, space, ncoefficients, coefficients, setcanonicaldomain, evaluate
import Base: +, -, *, /, inv, literal_pow, getindex, angle, exponent, ==, convert
import HolomorphicFun: poles, roots
import MacroTools: postwalk, @capture
import SpecialFunctions: gamma

export WienerHopfKernel
export @wienerhopf

abstract type WienerHopfKernel <: AbstractScalarFunction end

include("macros/wienerhopf.jl")

end
