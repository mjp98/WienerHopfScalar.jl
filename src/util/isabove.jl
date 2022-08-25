isabove(x::Fun, z) = isabove(space(x), z)
isabove(x::Space, z) = isabove(domain(x), z)
isabove(d::Line, z) = isaboveline(d.center, angle(d), z)
function isaboveline(c, θ, z)
    w = cis(-θ) * (z - c)
    return imag(w) > 0 || (imag(w) == 0 && real(w) > 0)
end

forceabove(z, d) = isabove(d, z) ? z : -z
forceabove(z, d, u) = (u && isabove(d, z)) ? z : -z
