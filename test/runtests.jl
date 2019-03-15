using TropicalHomotopyContinuation
const THC = TropicalHomotopyContinuation
import PolynomialTestSystems: equations, cyclic, ipp2
using Test

@testset "TropicalHomotopyContinuation" begin
    A₁ = [0 0 1 1; 0 2 0 1]
    A₂ = [0 0 1 2; 0 1 1 0]

    A = cayley(A₁, A₂)

    @test A == [0  0  1  1  0  0  1  2
                0  2  0  1  0  1  1  0
                1  1  1  1  0  0  0  0
                0  0  0  0  1  1  1  1]

    @test_throws ErrorException cayley([1 0; 0 1], [1 0; 0 0; 0 2])

    w₁ = [0, 0, 0, -2]
    w₂ = [0, -3, -4, -8]
    w = [w₁; w₂]

    v₁ = [0, 0, 0, -1]
    v₂ = copy(w₂)
    v = [v₁; v₂]


    mixed_cell_indices = [(2, 3), (1, 3)]
    indexing = THC.CayleyIndexing(size.((A₁, A₂), 2))
    ord = THC.DotOrdering(w)
    cell = MixedCell(mixed_cell_indices, A, indexing)
    @test cell.volume == 3
    @test cell.circuit_table == [1 2; 0 0; 0 0; 1 -1; 0 0; 1 2; 0 0; -2 -1]
    ineq = THC.first_violated_inequality(cell, v, ord)
    @test ineq.config_index == 1
    @test ineq.col_index == 4

    @test THC.exchange_column(cell, THC.exchange_first, ineq) == MixedCell([(4, 3), (1, 3)], A, indexing)
    @test THC.exchange_column(cell, THC.exchange_second, ineq) == MixedCell([(2, 4), (1, 3)], A, indexing)

    ind_back = THC.reverse_index(ineq, cell, THC.exchange_second)
    cell2 = THC.exchange_column(cell, THC.exchange_second, ineq)
    @test cell2.volume == 2
    @test cell == THC.exchange_column(cell2, THC.exchange_second, ind_back)

    ind_back = THC.reverse_index(ineq, cell, THC.exchange_first)
    cell2 = THC.exchange_column(cell, THC.exchange_first, ineq)
    @test cell2.volume == 1
    @test cell == THC.exchange_column(cell2, THC.exchange_first, ind_back)

	@test mixed_volume(equations(cyclic(5))) == 70
	@test mixed_volume(equations(cyclic(7))) == 924
	@test mixed_volume(equations(ipp2())) == 288
end
