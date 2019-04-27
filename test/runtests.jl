using MixedSubdivisions
const THC = MixedSubdivisions
import PolynomialTestSystems: equations, cyclic, ipp2
using Test

@testset "MixedSubdivisions" begin
	@testset "Basic" begin
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
	    ord = THC.DotOrdering(Int32.(w))
	    cell = MixedCellTable(mixed_cell_indices, A, indexing)
	    @test cell.volume == 3
	    @test cell.circuit_table == [1 2; 3 0; 0 0; 1 -1; 0 3; 1 2; 0 0; -2 -1]
	    ineq = THC.first_violated_inequality(cell, v, ord)
	    @test ineq.config_index == 1
	    @test ineq.col_index == 4

	    @test THC.exchange_column(cell, THC.exchange_first, ineq) == MixedCellTable([(4, 3), (1, 3)], A, indexing)
	    @test THC.exchange_column(cell, THC.exchange_second, ineq) == MixedCellTable([(2, 4), (1, 3)], A, indexing)

	    ind_back = THC.reverse_index(ineq, cell, THC.exchange_second)
	    cell2 = THC.exchange_column(cell, THC.exchange_second, ineq)
	    @test cell2.volume == 2
	    @test cell == THC.exchange_column(cell2, THC.exchange_second, ind_back)

	    ind_back = THC.reverse_index(ineq, cell, THC.exchange_first)
	    cell2 = THC.exchange_column(cell, THC.exchange_first, ineq)
	    @test cell2.volume == 1
	    @test cell == THC.exchange_column(cell2, THC.exchange_first, ind_back)
	end

	@testset "Mixed Volume" begin
		@test mixed_volume(equations(cyclic(5))) == 70
		@test mixed_volume(equations(cyclic(7))) == 924
		@test mixed_volume(equations(cyclic(10))) == 35940
		@test mixed_volume(equations(cyclic(11))) == 184756
		@test mixed_volume(equations(ipp2())) == 288

		@test mixed_volume(equations(cyclic(5)), algorithm=:total_degree) == 70
		@test mixed_volume(equations(ipp2()), algorithm=:total_degree) == 288
	end

	@testset "Mixed Cells" begin
		A₁ = [0 0 1 1; 0 2 0 1]
		A₂ = [0 0 1 2; 0 1 1 0]

		A = [A₁, A₂]

		w₁ = [0, 0, 0, 2]
		w₂ = [0, 3, 4, 8]
		w = [w₁, w₂]

		v₁ = [0, 0, 0, 1]
		v₂ = copy(w₂)
		v = [v₁, v₂]

		cells_v = mixed_cells(A, v)
		@test length(cells_v) == 3
		@test sum(volume, cells_v) == 4
		@test sort(volume.(cells_v)) == [1, 1, 2]
		@test all(cells_v) do c
		    all(sort.(map(A_i -> A_i' * normal(c), A) .+ v)) do r
		        isapprox(r[1], r[2], atol=1e-12)
		    end
		end

		cells_w = mixed_cells(A, w)
		@test length(cells_w) == 2
		@test sum(volume, cells_w) == 4
		@test sort(volume.(cells_w)) == [1, 3]
		@test all(cells_w) do c
		    all(sort.(map(A_i -> A_i' * normal(c), A) .+ w)) do r
		        isapprox(r[1], r[2], atol=1e-12)
		    end
		end
	end
end
