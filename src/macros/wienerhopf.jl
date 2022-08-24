
isplus(x::Symbol) = last(string(x)) ∈ ['₊', '⁺']
isminus(x::Symbol) = last(string(x)) ∈ ['₋', '⁻']

del_lastchar!(x::Symbol) = Symbol(chop(string(x)))

macro wienerhopf(ex)
    ex = postwalk(parse_pm_subscripts, ex)
    return esc(ex)
end

function parse_pm_subscripts(x)
    if @capture(x, f_(xs__))
        if isminus(f)
            return :($(del_lastchar!(f))($(xs...), false))
        elseif isplus(f)
            return :($(del_lastchar!(f))($(xs...), true))
        end
    end
    return x
end
