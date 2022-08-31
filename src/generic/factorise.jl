# Generic factorization

function factorise(K::WienerHopfKernel, sp::Space)
    A = isolate_inf(K)
    R = isolate_poleroot(K, Line()) # Split contour is real axis.
    normal_kernel = K/(A*R)

    return (logfactorise(K / (A * R), sp) * A) * R
end
factorise(K::WienerHopfKernel) = factorise(K, defaultspace(K))
