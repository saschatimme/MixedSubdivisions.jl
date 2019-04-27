using MixedSubdivisions
using PolynomialTestSystems
const THC = MixedSubdivisions
using Profile


f = equations(cyclic(5))
A = THC.support(f)
lift = map(A) do A_i
    rand(0:100, size(A_i, 2))
end

iter = MixedCellIterator(A, lift)
@time THC.iterate_mixed_cells(iter)
@time THC.iterate_mixed_cells(iter)
Profile.clear_malloc_data()
THC.iterate_mixed_cells(iter)
