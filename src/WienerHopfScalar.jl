module WienerHopfScalar

using HolomorphicFun
using MacroTools

import MacroTools: postwalk, @capture

export WienerHopfKernel
export @wienerhopf

abstract type WienerHopfKernel <: AbstractScalarFunction end

include("macros/wienerhopf.jl")

end
