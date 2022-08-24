
# Noble 1.11
struct NobleKernel <: WienerHopfKernel end

function evaluate(::NobleKernel,α)
	πα = π*α
	return πα*coth(πα)
end

function evaluate(K::NobleKernel,α,u)
	if u
		return sqrt(π)*gamma(1-im*α)/gamma(0.5-im*α)
	else
		return K(-α,true)
	end
end

factorise(K::NobleKernel,args...) = K
