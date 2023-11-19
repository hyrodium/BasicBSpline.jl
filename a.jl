using Revise
using LinearAlgebra
using StaticArrays
using SparseArrays
using Test
using BasicBSpline
using BenchmarkTools
using Plots
import BasicBSpline._lower_R
import BasicBSpline._vec
import BasicBSpline._generate_flags
import BasicBSpline._changebasis_I
import BasicBSpline._changebasis_R

P = BSplineSpace{3}(knotvector"4 1  2 1 14")
P′= BSplineSpace{4}(knotvector"5 2  3 2 25")
P ⊑ P′
P ⊆ P′
isnondegenerate_R(P)
isnondegenerate_R(P′)

changebasis_I(P, P′)
changebasis_I(P, P′)-changebasis_R(P, P′)

changebasis_I(P, P)
changebasis_I(P, P)-changebasis_R(P, P)

P = BSplineSpace{3}(knotvector"13 1  2 1 14")
P′= BSplineSpace{4}(knotvector" 5 2  3 2 25")
P ⊑ P′
P ⊆ P′

changebasis_I(P, P′)
changebasis_I(P, P′)-changebasis_R(P, P′)

changebasis_I(P, P)
changebasis_I(P, P)-changebasis_R(P, P)

# common
U = Float64
k = knotvector(P)
k′ = knotvector(P′)
p = degree(P)
p′ = degree(P′)
n = dim(P)
n′ = dim(P′)

### on R
Aᵖ = changebasis_R(P, P′)
Aᵖ⁻¹ = _changebasis_R(_lower_R(P), _lower_R(P′))

K′ = [k′[i+p′] - k′[i] for i in 1:n′+1]
K = U[ifelse(k[i+p] ≠ k[i], U(1 / (k[i+p] - k[i])), zero(U)) for i in 1:n+1]

i = 2
v1 = [Aᵖ[i,j] - Aᵖ[i,j-1] for j in 2:n′]
v2 = [(p/p′)*(K′[j])*(K[i]*Aᵖ⁻¹[i,j]-K[i+1]*Aᵖ⁻¹[i+1,j]) for j in 2:n′]
norm(v1-v2)

### on I
Aᵖ = changebasis_R(P, P′)
Aᵖ⁻¹ = _changebasis_R(_lower_R(P), _lower_R(P′))[2:n, 2:n′]

K′ = [k′[i+p′] - k′[i] for i in 1:n′+1]
K = U[ifelse(k[i+p] ≠ k[i], U(1 / (k[i+p] - k[i])), zero(U)) for i in 1:n+1]

i = n-1
v1 = [Aᵖ[i,j] - Aᵖ[i,j-1] for j in 2:n′]
v2 = [(p/p′)*(K′[j])*(K[i]*Aᵖ⁻¹[i-1,j-1]-K[i+1]*Aᵖ⁻¹[i,j-1]) for j in 2:n′]
norm(v1-v2)








function ranges(P::BSplineSpace{p,T,<:AbstractKnotVector{T}}, P′::BSplineSpace{p′,T′,<:AbstractKnotVector{T′}}) where {p,p′,T,T′}
    k = knotvector(P)
    k′ = knotvector(P′)
    n = dim(P)
    j_begin = 1
    j_end = 1
    ranges = fill(1:0, n)
    for i in 1:n
        # Skip for degenerated basis
        isdegenerate_I(P,i) && continue
        k_ = view(k, i:i+p+1)
        kᵢ_I = max(k[i], k[1+p])
        kᵢ₊ₚ₊₁_I = min(k[i+p+1], k[end-p])
        j_begin = findlast(==(kᵢ_I), _vec(k′)) - count(≤(kᵢ_I), _vec(k_)) - p′ + p + 1
        j_end = findfirst(==(kᵢ₊ₚ₊₁_I), _vec(k′)) + count(≥(kᵢ₊ₚ₊₁_I), _vec(k_)) - p′ - 1
        j_range = j_begin:j_end
        ranges[i] = j_range
    end
    return ranges
end

using BasicBSpline
p = 3
P = BSplineSpace{p}(KnotVector(1:8))
P′ = BSplineSpace{p}(KnotVector([2,3,3,4,5,6,7,8]))

k = knotvector(P)
k′ = knotvector(P′)
p = degree(P)
p′ = degree(P′)

i = 1
k_ = view(k, i:i+p+1)
kᵢ_I = max(k[i], k[1+p])
kᵢ₊ₚ₊₁_I = min(k[i+p+1], k[end-p])
j_begin = findlast(==(kᵢ_I), _vec(k′)) - count(≤(kᵢ_I), _vec(k_)) - p′ + p + 1
j_end = findfirst(==(kᵢ₊ₚ₊₁_I), _vec(k′)) + count(≥(kᵢ₊ₚ₊₁_I), _vec(k_)) - p′ - 1
j_range = j_begin:j_end






k_ = view(k, i:i+p+1)
kᵢ₊ₚ₊₁_I = min(k[i+p+1], k[end-p])

findfirst(==(kᵢ₊ₚ₊₁_I), _vec(k′)) + count(≥(kᵢ₊ₚ₊₁_I), _vec(k_)) - 1 - p′ - 1




k′

[12345]
[23345678]

ranges[i] = j_range




plot(P1, xlims=extrema(domain(P1)))
plot!(P2)


ranges(P1, P2)
dim(P1)
dim(P2)

