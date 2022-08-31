# Generic factorization

function factorise(kernel::WienerHopfKernel, sp::Space)
    asymptotes = isolate_inf(kernel)
    rational = isolate_poleroot(kernel, Line()) # Split contour is real axis.
    
    normalised_kernel = kernel / (asymptotes * rational)

    return (logfactorise(normalised_kernel, sp) * asymptotes) * rational
end
factorise(K::WienerHopfKernel) = factorise(K, defaultspace(K))
