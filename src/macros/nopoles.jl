macro noroots(ex)
    return :(candidate_roots(x::$(esc(ex)){T},args...) where T = T[])
end

macro nopoles(ex)
    return :(candidate_poles(x::$(esc(ex)){T},args...) where T = T[])
end
