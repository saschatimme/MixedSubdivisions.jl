using TropicalHomotopyContinuation
using Test

@testset "TropicalHomotopyContinuation" begin
    A₁ = [0 1 2 1; 0 0 0 1]
    A₂ = [1 0 1 2; 0 1 1 1]
    C = Cayley(A₁, A₂)

    @test C.A == [0 1 2 1 1 0 1 2
                  0 0 0 1 0 1 1 1
                  1 1 1 1 0 0 0 0
                  0 0 0 0 1 1 1 1]

    @test C.mᵢs == [4, 4]

    @test_throws ErrorException Cayley([1 0; 0 1], [1 0; 0 0; 0 2])
end
