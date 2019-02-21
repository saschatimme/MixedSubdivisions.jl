using TropicalHomotopyContinuation
const THC = TropicalHomotopyContinuation

A₁ = [0 0 1 1
      0 2 0 1]
w₁ = [0, 0, 0, -2]
v₁ = [0, 0, 0, -1]

A₂ = [0 0 1 2
      0 1 1 0]

w₂ = [0, -3, -4, -8]
v₂ = copy(w₂)

C = THC.cayley((A₁, A₂))
w = [w₁; w₂]
v = [v₁; v₂]

mixed_cells = [[(2, 3), (1, 3)], [(3, 4), (3, 4)]]

subdivison = THC.MixedSubdivision([A₁, A₂], mixed_cells)
ord = THC.DotOrdering(w)

violated_ind = THC.first_violated_inequality(subdivison.mixed_cells[1], v, ord)

THC.mixed_cell_split(subdivison.mixed_cells[1], violated_ind)


subdivison.mixed_cells[1]



indexing = THC.ConfigurationIndexing([4, 3])
for (i, j, offset) in indexing
      @show i, j, offset
end
collect(indexing)




C = Cayley(A₁, A₂)

D = THC.mixed_cell_submatrix(C, mixed_cells[1])

det(D)

using LinearAlgebra

LUD = lu(D)


inv(D) * C.A


THC.MixedCell(C, mixed_cells[1])
