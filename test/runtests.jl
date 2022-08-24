using WienerHopfScalar
using Test

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
end
