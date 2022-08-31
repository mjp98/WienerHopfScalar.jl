# Generic factorization

function factorise(K::WienerHopfKernel,sp::Space)
    L = isolate_inf(K)
    R = isolate_poleroot(K, Line())
    return (logfactorise(K / (L * R), sp) * L) * R
end
factorise(K::WienerHopfKernel) = factorise(K, defaultspace(K))
