complexfloat(z) = complex(float(z))
const cfloat = complexfloat
const iπ = im * π
const i2π = 2iπ

cphase(z) = log(complex(z)) / i2π

leftvalue(z, ε) = abs(z) < ε ? -ε : z - ε * abs(z)
rightvalue(z, ε) = abs(z) < ε ? ε : z + ε * abs(z)
lsign(n::Bool) = n ? 1 : -1

include("factorise.jl")
include("isabove.jl")
include("pair.jl")
include("roots.jl")
