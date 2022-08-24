function factorise(K::ScalarConstant, args...)
    K₊ = ScalarConstant(sqrt(complex(K.value)))
    WienerHopfPair(K₊, K₊)
end

value(K::WienerHopfPair{<:ScalarConstant,<:ScalarConstant}) = K(1)

*(x::Number, K::WienerHopfKernel) = *(factorise(ScalarConstant(x)), K)
*(K::WienerHopfKernel, x::Number) = *(x, K)
