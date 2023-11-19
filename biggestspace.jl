using LinearAlgebra
using SparseArrays
using BasicBSpline
using Test
function Base.print_array(io::IO, S::SparseArrays.AbstractSparseMatrixCSCInclAdjointAndTranspose)
    if max(size(S)...) < 52
        Base.print_matrix(io, S)
    else
        _show_with_braille_patterns(io, S)
    end
end

Pd1 = BSplineSpace{3}(knotvector"1111 2  3 44")
Pd2 = BSplineSpace{4}(knotvector" 3 2 3 26 55")
Pd3 = BSplineSpace{4}(knotvector"2322 3 26 55")

Pd1 = BSplineSpace{3}(knotvector"44 3  2 121 ")
Pd2 = BSplineSpace{4}(knotvector"55 62 3 2121")
Pd3 = BSplineSpace{4}(knotvector"55 62 3 2321")
n1 = dim(Pd1)
n2 = dim(Pd2)
n3 = dim(Pd3)
A12 = changebasis_I(Pd1, Pd2)
A13 = changebasis_R(Pd1, Pd3)
A23 = changebasis_R(Pd2, Pd3)

A12 * A23[:,1:n2] - A13[:,1:n2]

A12


plot(Pd1)
plot!(Pd2)
plot!(Pd3)

savefig("Pd123.png")



Pd1 = BSplineSpace{3}(knotvector"44 3  2 121 ")
Pd2 = BSplineSpace{4}(knotvector"55 62 3 2121")
changebasis_I(Pd1, Pd2)

function _changebasis_nonzero_R(P1, P2)
end



