using MixedSubdivisions
const MS = MixedSubdivisions
import MultivariatePolynomials
const MP = MultivariatePolynomials
import PolynomialTestSystems: equations, cyclic, ipp2, cyclooctane
using Test

@testset "MixedSubdivisions" begin
    @testset "Basic" begin
        A₁ = [0 0 1 1; 0 2 0 1]
        A₂ = [0 0 1 2; 0 1 1 0]

        A = MS.cayley(A₁, A₂)

        @test A == [
            0 0 1 1 0 0 1 2
            0 2 0 1 0 1 1 0
            1 1 1 1 0 0 0 0
            0 0 0 0 1 1 1 1
        ]

        @test_throws ErrorException MS.cayley([1 0; 0 1], [1 0; 0 0; 0 2])

        w₁ = [0, 0, 0, -2]
        w₂ = [0, -3, -4, -8]
        w = [w₁; w₂]

        v₁ = [0, 0, 0, -1]
        v₂ = copy(w₂)
        v = [v₁; v₂]


        mixed_cell_indices = [(2, 3), (1, 3)]
        indexing = MS.CayleyIndexing(size.((A₁, A₂), 2))
        ord = MS.DotOrdering(Int32.(w))
        cell = MS.MixedCellTable(mixed_cell_indices, A, indexing)
        @test cell.volume == 3
        @test cell.circuit_table == [1 2; 3 0; 0 0; 1 -1; 0 3; 1 2; 0 0; -2 -1]
        ineq = MS.first_violated_inequality(cell, v, ord, typemax(Int32))
        @test ineq.config_index == 1
        @test ineq.col_index == 4

        @test MS.exchange_column(cell, MS.exchange_first, ineq, 8) == MS.MixedCellTable(
            [(4, 3), (1, 3)],
            A,
            indexing
        )
        @test MS.exchange_column(cell, MS.exchange_second, ineq, 8) == MS.MixedCellTable(
            [(2, 4), (1, 3)],
            A,
            indexing,
        )

        ind_back = MS.reverse_index(ineq, cell, MS.exchange_second)
        cell2 = MS.exchange_column(cell, MS.exchange_second, ineq, 8)
        @test cell2.volume == 2
        @test cell == MS.exchange_column(cell2, MS.exchange_second, ind_back, 8)

        ind_back = MS.reverse_index(ineq, cell, MS.exchange_first)
        cell2 = MS.exchange_column(cell, MS.exchange_first, ineq, 8)
        @test cell2.volume == 1
        @test cell == MS.exchange_column(cell2, MS.exchange_first, ind_back, 8)
    end

    @testset "Mixed Volume" begin
        @test mixed_volume(equations(cyclic(5))) == 70
        @test mixed_volume(equations(cyclic(7))) == 924
        @test mixed_volume(equations(cyclic(10))) == 35940
        @test mixed_volume(equations(cyclic(11))) == 184756
        @test mixed_volume(equations(ipp2())) == 288

        @test mixed_volume(equations(cyclic(5)), algorithm = :total_degree) == 70
        @test mixed_volume(equations(ipp2()), algorithm = :total_degree) == 288
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
                isapprox(r[1], r[2], atol = 1e-12)
            end
        end

        cells_w = mixed_cells(A, w)
        @test length(cells_w) == 2
        @test sum(volume, cells_w) == 4
        @test sort(volume.(cells_w)) == [1, 3]
        @test all(cells_w) do c
            all(sort.(map(A_i -> A_i' * normal(c), A) .+ w)) do r
                isapprox(r[1], r[2], atol = 1e-12)
            end
        end
    end

    @testset "Fine mixed cells" begin
        f = equations(cyclic(7))
        cells, lift = fine_mixed_cells(f)
        @test sum(c -> c.volume, cells) == 924
        @test lift isa Vector{Vector{Int32}}

        @test isnothing(fine_mixed_cells(f, max_tries = 0))
    end

    @testset "Overflow error messages" begin
        f = equations(cyclooctane())
        @test_throws ArgumentError MS.mixed_volume(f)
        F = [f; randn(2, 18) * [MP.variables(f); 1]]
        A = support(F)
        lifting = map(Ai -> rand(-Int32(2^20):Int32(2^20), size(Ai, 2)), A)
        cells = mixed_cells(A, lifting)
        @test_deprecated is_fully_mixed_cell(first(cells), A, lifting)
        @test all(is_fine, cells)
        @test sum(volume, cells) == 32768
        @test sum(volume, first(fine_mixed_cells(F))) == 32768
    end

    @testset "Bugfix#1" begin
        lift = Array{Int32,1}[
            [-312, 347, 823, 901],
            [-1021, -296, -421, -491, -9, 168, -905, -79, -610, -898],
            [489, 385, -540, 612, 792, -986, 697, 695, 842, -731],
            [-126, 677, -699, 96, -270, 186, 263, 75, -765, 68],
            [-230, 173, 142, 531, 218, 668, -305, 177, -710, -94],
            [283, -967, 530, 226, -452, -557, 501, -828],
            [-656, -75, -709, -588],
            [407, 920, 330, 265],
            [157, 2, 524, 417, -449, -897, -19, -603, -209, 257],
            [-683, 763, 1023, 983, -863, -727, 211, -466, 998, -644],
            [-916, 661, 401, 808, 979, -793, -815, 155, -1018, -409],
            [1015, -274, -71, 955, 100, 433, 249, -552],
            [-523, -225, 997, -784, -884],
            [457, -579, 749],
            [624, -74, -563, 861, 86],
            [
             -428,
             75,
             274,
             283,
             691,
             -782,
             513,
             -915,
             -189,
             279,
             793,
             624,
             -741,
             -698,
             -94,
             945,
             -58,
             -164,
            ],
            [
             533,
             -156,
             -488,
             539,
             718,
             724,
             -32,
             -845,
             441,
             -273,
             137,
             511,
             -1005,
             -856,
             -788,
             824,
             279,
             -12,
            ],
        ]
        f = equations(cyclooctane())
        F = [f; randn(2, 18) * [MP.variables(f); 1]]
        A = support(F)
        cells = mixed_cells(A, lift)
        @test !isempty(cells)
        sum(volume, cells) == 32768
    end
end
