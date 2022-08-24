struct WienerHopfNaN <: WienerHopfKernel end
evaluate(::WienerHopfNaN,args...) = NaN
