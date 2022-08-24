function factorise(x::AbstractVector, d::Domain)
    p, m = copy(x), copy(x)
    filter!(z -> isabove(d, z), p)
    filter!(z -> !isabove(d, z), m)
    return p, m
end
