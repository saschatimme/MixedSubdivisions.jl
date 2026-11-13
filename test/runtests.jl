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
        for I in cell.indexing
            MS.is_valid_inquality(cell, I) || continue
            @test MS.inequality_dot(cell, I, ord.w) == MS.inequality(cell, I)' * ord.w
        end
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
        @test sum(c -> volume(c), MixedCellIterator(A, w; algorithm = :total_degree)) == 4
        @test all(cells_w) do c
            all(sort.(map(A_i -> A_i' * normal(c), A) .+ w)) do r
                isapprox(r[1], r[2], atol = 1e-12)
            end
        end
    end

    @testset "MixedCellIterator" begin
        A = support(equations(cyclic(5)))
        lift = map(Ai -> rand(-Int32(2^15):Int32(2^15), size(Ai, 2)), A)
        expected = sort(indices.(mixed_cells(A, lift)))
        iter = MixedCellIterator(A, lift)
        # each iteration returns its own cell
        @test sort(indices.(collect(iter))) == expected
        # a new iteration after stopping early starts from scratch
        first(iter)
        @test sort(indices.(collect(iter))) == expected
        @test sum(volume, MixedCellIterator(A, lift; algorithm = :total_degree)) == 70
    end

    @testset "is_fine and extreme liftings" begin
        # A lifted point lies 1 above the face of a cell: all cells are fine,
        # but a floating point check reported one as not fine
        A = [[6 0 4 8; 7 0 3 8], [3 10 5 1; 5 6 5 4]]
        lift = [[862884, 393514, -929749, 1038738], [355334, -832450, -991389, -804171]]
        cells = mixed_cells(A, lift)
        @test sum(volume, cells) == 56
        @test all(is_fine, cells)

        # a lifting value whose negation overflows Int32
        A = [[0 0 1 1; 0 2 0 1], [0 0 1 2; 0 1 1 0]]
        @test_throws OverflowError mixed_cells(A, [[0, 0, 0, typemin(Int32)], [0, 3, 4, 8]])
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
        @test sum(volume, cells) == 32768
    end
    @testset "Bugfix#30" begin
        As = Matrix{Int64}[
            [1 0 1 0 1 1 0; 0 1 1 0 1 1 0; 0 0 0 1 -1 0 0],
            [1 1 0; 1 1 0; -1 0 0],
            [0 1 0; 1 1 0; 0 0 0]
        ]
        weights = [
            [2, 4, 1, 4, 5, 4, 3],
            [1, 4, 3],
            [5, 5, 2]
        ]
        @test mixed_volume(As, algorithm=:regeneration) == 3
        @test mixed_volume(As, algorithm=:total_degree) == 3
        cells = mixed_cells(As, weights)
        @test length(cells) == 3
        @test volume.(cells) == [1, 1, 1]
        @test sort(normal.(cells)) == [[-4.0, 1.0, -3.0], [2.0, -3.0, -4.0], [3.0, -3.0, -2.0]]
    end

    @testset "Overflow bounds" begin
        # Underestimated bounds let the Int32 fast paths wrap around silently.
        # Expected values verified by brute-force enumeration of all cells.
        A = [
            [0 24 12 9 4 9; 1 5 2 0 12 11; 15 6 26 24 2 23],
            [0 25 22 0; 9 23 24 24; 14 25 11 7],
            [16 20 18 18; 12 24 20 15; 6 28 0 30],
        ]
        lift = [
            [786261, 716467, 533220, -891724, 877712, 124907],
            [-69575, 347349, 297421, -368259],
            [-187312, -460139, -867041, 868538],
        ]
        cells = mixed_cells(A, lift)
        @test length(cells) == 11
        @test sum(volume, cells) == 21123

        A = [
            [8 1 6 8 5 7 4 5 7; 1 7 10 8 6 5 9 10 7; 8 2 10 1 8 5 2 10 10; 4 5 6 8 6 10 7 8 7; 6 7 5 6 2 3 8 4 9],
            [0 1 5 1 2 10 0 2 2; 9 3 9 9 8 1 0 0 6; 10 9 1 8 10 6 10 3 4; 9 5 4 9 9 3 6 7 8; 0 10 2 2 4 9 5 0 3],
            [2 6 0 6 8 6; 5 1 1 9 9 6; 8 7 10 2 4 4; 10 2 2 10 6 9; 7 3 6 7 10 10],
            [0 2 8 0 1 6 1 5 4; 10 0 10 0 8 0 4 6 9; 6 1 1 4 0 7 5 9 4; 1 2 5 5 9 3 5 7 8; 4 0 5 10 0 8 8 4 9],
            [5 2 8 7 10 6 7; 4 2 9 0 10 10 9; 3 6 7 0 7 1 2; 5 0 3 8 8 10 1; 4 7 4 3 10 7 1],
        ]
        @test mixed_volume(A, algorithm = :regeneration) == 1119513
        @test mixed_volume(A, algorithm = :total_degree) == 1119513
    end

    @testset "Int64 circuit tables" begin
        # These need circuit tables beyond Int32, which used to throw an OverflowError.
        # Expected values verified by brute-force enumeration of all cells.
        A = [
            [28 34 27 80 59 98; 41 28 45 103 66 14; 22 36 117 108 58 11; 40 98 81 99 52 43],
            [75 65 93 91 46 46; 37 1 78 119 88 7; 115 65 68 17 41 67; 70 62 74 118 9 50],
            [106 35 2 39 33 79; 31 89 26 15 59 107; 27 79 76 6 30 16; 38 82 47 38 108 39],
            [84 87 8 33 90 47; 3 11 55 61 26 31; 3 44 104 81 39 17; 41 45 74 111 18 113],
        ]
        @test mixed_volume(A, algorithm = :regeneration) == 404441984
        @test mixed_volume(A, algorithm = :total_degree) == 404441984
        @test mixed_volume(A, threads = true) == 404441984

        # lifting too large for Int32 dot products
        A = [
            [388 91 22 62 332; 273 258 70 105 263; 73 319 87 109 358],
            [98 277 214 204 394; 144 279 124 146 8; 98 318 106 20 15],
            [111 106 283 184 247; 308 76 11 253 381; 241 184 249 282 243],
        ]
        lift = [
            [-191983741, -83874253, -174357837, 178798983, -316331603],
            [417735239, -536306775, -505198336, 355914830, 333007717],
            [-424101032, 158185279, -335004225, -238058356, -183751679],
        ]
        @test sum(volume, mixed_cells(A, lift)) == 62259741
        @test sum(volume, mixed_cells(A, lift; threads = true)) == 62259741
    end

    @testset "Multithreading" begin
        # CI runs the tests with several threads, see run_tests.yml
        f = equations(cyclic(10))
        @test mixed_volume(f; threads = true) == 35940
        @test mixed_volume(support(f); threads = true) == 35940
        @test_throws ArgumentError mixed_volume(f; threads = true, algorithm = :total_degree)

        A = support(equations(cyclic(9)))
        lift = map(Ai -> rand(-Int32(2^15):Int32(2^15), size(Ai, 2)), A)
        key(c) = (c.indices, c.volume, c.normal)
        cells = mixed_cells(A, lift; threads = true)
        @test sort(key.(cells)) == sort(key.(mixed_cells(A, lift)))

        cells, lift = fine_mixed_cells(f; threads = true, show_progress = false)
        @test sum(volume, cells) == 35940
        @test all(is_fine, cells)
    end


end