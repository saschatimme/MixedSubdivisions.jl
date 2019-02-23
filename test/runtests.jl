using TropicalHomotopyContinuation
using Test

@testset "TropicalHomotopyContinuation" begin
    A₁ = [0 1 2 1; 0 0 0 1]
    A₂ = [1 0 1 2; 0 1 1 1]
    C = cayley(A₁, A₂)

    @test C == [0 1 2 1 1 0 1 2
                0 0 0 1 0 1 1 1
                1 1 1 1 0 0 0 0
                0 0 0 0 1 1 1 1]

    @test_throws ErrorException cayley([1 0; 0 1], [1 0; 0 0; 0 2])
end
