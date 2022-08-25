using WienerHopfScalar
using Test

sp = Chebyshev(Line{-1/4}(0.0))
K = RobinKernel(randn(ComplexF64),randn(ComplexF64))
L = isolate_inf(K,sp)
@test L == GammaKernel(wavenumber(K))
@test WienerHopfScalar.roots(K) == WienerHopfScalar.roots(NormalisedRobinKernel(K))

k = 1

μ = 0.5
K = RobinKernel(k,μ)
L = isolate_inf(K)
R = isolate_poleroot(K,Line())
Kl = logfactorise(K/(R*L),sp)

K2 = factorise(K,sp)
