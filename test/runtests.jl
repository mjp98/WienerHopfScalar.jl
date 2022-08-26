using WienerHopfScalar
using ApproxFun
using Test

import WienerHopfScalar: forceabove, defaultscale, leftlimit, rightlimit, poles, roots, isabove

@testset "WienerHopfScalar.jl" begin

    @testset "remove ₊₋⁺⁻ macro" begin

        import WienerHopfScalar: isplus, isminus, del_lastchar!

        @testset "identify last character" begin
            @test isminus(:K₋)
            @test isplus(:K₊)
            @test !isminus(:K₊)
            @test !isplus(:K₋)
        end

        @testset "delete last character" begin
            @test del_lastchar!(:K₋) == :K
        end

        function K(z, u)
            u ? sin(z) : cos(z)
        end
        function K(z)
            K(z, true) * K(z, false)
        end
        z = randn(ComplexF64)
        k = randn(ComplexF64)

        @testset "subscript variations value" begin
            @test @wienerhopf K₊(z) == K(z, true)
            @test @wienerhopf K⁺(z) == K(z, true)
            @test @wienerhopf K₋(z) == K(z, false)
            @test @wienerhopf K⁻(z) == K(z, false)
        end

        @testset "expression including an unchanged value" begin
            @test @wienerhopf K₊(z) + K₋(z) + k == K(z, true) + K(z, false) + k
        end

        @testset "macroexpand" begin
            @test (@macroexpand @wienerhopf K₊(z) + K₋(z) + k) == :(K(z, true) + K(z, false) + k)
        end

    end


    @testset "generic" begin
        struct TestKernel{T} <:WienerHopfKernel end
        K = TestKernel{Float64}()
        @test_throws ErrorException evaluate(K,0.2)
        @test_throws ErrorException evaluate(K,0.2, true)
        @test_throws ErrorException leftlimit(K)
        @test_throws ErrorException rightlimit(K)
        @test defaultscale(K) == 1
        @test isabove(Line(),WienerHopfScalar.defaultpoint(K,true))
    end


    @testset "nopoles macro" begin
        @test begin x = @nopoles TestKernel; r = TestKernel{Float64}(); isempty(poles(r)); end
        @test begin x = @noroots TestKernel; r = TestKernel{Float64}(); isempty(roots(r)); end

    end

    @testset "isolate_inf" begin
        ∞₋ = -1e15
        ∞₊ =  1e15

        a = randn(ComplexF64)
        b = randn(ComplexF64)

        z = randn(ComplexF64)

        import WienerHopfScalar: one2zero_phase

        f = isolate_inf(a,b;z₊ = 1im, z₋ = -10im)

        @test one2zero_phase(f,∞₊) ≈ 0 atol = 1e-12
        @test one2zero_phase(f,∞₋) ≈ 1 atol = 1e-12
        @test f(∞₋) ≈ a
        @test f(∞₊) ≈ b
        @test f(z,true)*f(z,false) ≈ f(z)  atol = 100*eps()
        @test ncoefficients(Fun(f,Line())) < 2048
    end

    @testset "isabove" begin
        import WienerHopfScalar: isabove
		# Real line
		@test  isabove(Fun(tanh,Line()), 1im)
		@test !isabove(Fun(tanh,Line()),-1im)
		@test  isabove(Fun(tanh,Line{-1/4}()), 1)
		@test !isabove(Fun(tanh,Line{-1/4}()),-1)
        @test  isabove(Fun(tanh,Chebyshev(Line())), 1im) == isabove(Line(), 1im)
        @test  isabove(Chebyshev(Line()), 1im) == isabove(Line(), 1im)

        @test forceabove(1im,Line()) == 1im
        @test forceabove(-1im,Line()) == 1im
        @test forceabove(1im,Line(),true) == 1im
        @test forceabove(1im,Line(),false) == -1im
	end

    @testset "index" begin
        K = WienerHopfScalar.isolate_index(0)
        z = randn(ComplexF64)
        @test K(z) == 1
        @test K(z,true) == 1
        @test K(z,false) == 1
    end

    @testset "ConstantFunction" begin
        z = -1
        F = factorise(WienerHopfScalar.ScalarConstant(z))
        @test F(z) == z
        @test F(z,true) == F(z,false)
        @test F(z,true)*F(z,false) == F(z)
        @test factorise(F) == F
        @test WienerHopfScalar.value(F) == z

        u = randn(ComplexF64,2)
        U = [WienerHopfScalar.ScalarConstant(v) for v in u]
        @test prod(U) == WienerHopfScalar.ScalarConstant(prod(u))

        w = randn(ComplexF64)
        γ = GammaKernel(randn(ComplexF64))
        @test (2*γ)(w) ≈ 2*γ(w)
        @test (γ*2)(w) ≈ 2*γ(w)
    end

    @testset "split" begin
        p,m = factorise([im,-im,1+im],Line())
        @test Set(p) == Set([im,1+im])
        @test Set(m) == Set([-im])
    end

    @testset "GammaKernel" begin

        k = 1
        γ = GammaKernel(k)
        @test real(γ( 2k))>0
        @test real(γ(-2k))>0
        @test @wienerhopf γ₊(-k) == 0
        @test @wienerhopf γ₋( k) == 0
        @test factorise(γ) == γ
        z = randn(ComplexF64)
        @test @wienerhopf γ₊(z)*γ₋(z) ≈ γ(z)
        @test (γ^4)(z) == γ(z)^4

        @test γ[true](0.2) ≈ γ(0.2,true)
        @test defaultscale(γ) == k
        @test factors(γ)[1] == γ[true]
    end

    @testset "NobleKernel" begin
        K = NobleKernel()
        z = randn(ComplexF64)
        @test K(z,true)*K(z,false) ≈ K(z)
        @test K(z,true) == K(-z,false)
        @test factorise(K) == K
    end

    @testset "robin kernel" begin
        sp = Chebyshev(Line{-1/4}(0.0))
        k = randn(ComplexF64)
        μ = randn(ComplexF64)
        K = RobinKernel(k,μ)
        L = isolate_inf(K,sp)
        @test L == GammaKernel(wavenumber(K))
        @test WienerHopfScalar.roots(K) == WienerHopfScalar.roots(NormalisedRobinKernel(K))
        @test GammaKernel(K,0.2) == GammaKernel(k)(0.2)


        @test WienerHopfScalar.roots_poly(K) == WienerHopfScalar.roots_poly(NormalisedRobinKernel(k,μ))



        k = 1
        μtest = [0.5,2,-2im,2im+2,-2im-2,-2im+2,-2+1im,-2+0.1im]
        for μ ∈ μtest

            K = RobinKernel(k,μ)
            L = isolate_inf(K)
            R = isolate_poleroot(K,Line())
            Kl = logfactorise(K/(R*L),sp)
            K1 = Kl*L*R
            K2 = factorise(K,sp)

            z = randn(ComplexF64)
            @test K2(z) ≈ K1(z)
            z = randn(ComplexF64)
            @test K2(z,true)*K2(z,false) ≈ K2(z)
            @test ncoefficients(Kl) < 512
            @test ncoefficients(Fun(z->K2(z,true),-0.5k..0.5k)) < 512
            @test ncoefficients(Fun(z->K2(z,false),-0.5k..0.5k)) < 512
            r = WienerHopfScalar.roots(K)
            r .+= 1e-12
            if !isempty(r)
                r₊,r₋ = factorise(r,Line())
                for (_,x) in enumerate(r₊)
                    @test isapprox(K2(x,false),0;atol = 1e-6)
                end
                for (_,x) in enumerate(r₋)
                    @test isapprox(K2(x,true),0;atol = 1e-6)
                end
            end
        end
    end

    @testset "log branch" begin
        K = factorise(RobinKernel(1,2im+0.3),Chebyshev(Line{-1/4}(0.0)))
        z = randn(ComplexF64)
        @test K(z,true)*K(z,false) ≈ K(z)
        L = NormalisedRobinKernel(RobinKernel(1,2im+0.3))
	    R = isolate_poleroot(L,Line())
        sp = Chebyshev(Line{-1/4}(0.0))
        H = logfactorise(L/R,sp)
	    @test ncoefficients(H) < 512
    end

    @testset "NaN kernel" begin
        K = WienerHopfNaN()
        @test isnan(K(0.2))
        @test isnan(K(0.2,true))
    end

    @testset "Rawlins kernel" begin
        k = randn(Float64)
        K = RawlinsKernel(k,2,3)
        @test WienerHopfScalar.wavenumber(K) == complex(k)

        k = 1

        # factorise(RawlinsKernel(4,-20im-20,2),Chebyshev(Line{-1/4}(0.0))) gives negative winding number
        μtest = [0.5,2,-2im,2im+2,2im-2,-2im-2,-2im+2]
        for μ ∈ μtest
            sp = Chebyshev(Line{-1/4}(0.0))
            K = RawlinsKernel(k,μ,-2k)
            L = isolate_inf(K)
            R = isolate_poleroot(K,Line())
            Kl = logfactorise((K/(R*L)),sp)
            K1 = Kl*L*R
            K2 = factorise(K,sp)

            # Add test for poles and roots, and check that they are handled properly...
            z = randn(ComplexF64)
            @test K2(z) ≈ K1(z)
            @test ncoefficients(Kl) < 512
            @test ncoefficients(Fun(z->K2(z,true),-0.5k..0.5k)) < 512
            @test ncoefficients(Fun(z->K2(z,false),-0.5k..0.5k)) < 512
            z = randn(ComplexF64)
            @test K2(z,true)*K2(z,false) ≈ K2(z)

            r = WienerHopfScalar.roots(K)
            r .+= 1e-12
            if !isempty(r)
                r₊,r₋ = factorise(r,Line())
                for (_,x) in enumerate(r₊)
                    @test isapprox(K2(x,false),0;atol = 1e-6)
                end
                for (_,x) in enumerate(r₋)
                    @test isapprox(K2(x,true),0;atol = 1e-6)
                end
            end
        end
    end


    @testset "Poroelastic kernel" begin
        sp = Chebyshev(Line{-1/4}(0.0))

        K = PoroelasticK()
        k = abs(wavenumber(K))
        H = factorise(K)
        z = randn(ComplexF64)
        @test H(z) ≈ K(z)
        @test H(z,true)*H(z,false) ≈ H(z)
        @test ncoefficients(Fun(z->H(z,true),-0.5abs(k)..0.5abs(k))) < 512
        @test ncoefficients(Fun(z->H(z,false),-0.5abs(k)..0.5abs(k))) < 512


        I = PoroelasticI(K)
        k = abs(wavenumber(I))
        H = factorise(I)
        z = randn(ComplexF64)
        @test H(z) ≈ I(z)
        @test H(z,true)*H(z,false) ≈ H(z)
        @test ncoefficients(Fun(z->H(z,true),-0.5abs(k)..0.5abs(k))) < 512
        @test ncoefficients(Fun(z->H(z,false),-0.5abs(k)..0.5abs(k))) < 512

        @test all(isabove(Line(),z) for z in WienerHopfScalar.poles⁺(K))



    end

    @testset "rightvalue" begin
        @test WienerHopfScalar.leftvalue(0.2,0.01) == WienerHopfScalar.rightvalue(0.2,-0.01)
    end


end
