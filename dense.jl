using BasicBSpline
using Plots

k1 = knotvector" 1 1 11  1 1 1 1"
k2 = knotvector"1 1 1 1  11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)


k1 = KnotVector(rand(8)) + 4*KnotVector([0,1])
k2 = k1 + KnotVector(rand(3))
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)

BasicBSpline._changebasis_I(P1,P2)
BasicBSpline._changebasis_R(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)


BasicBSpline._changebasis_R(P1,P2) ≈ BasicBSpline._changebasis_I_new(P1,P2)

BasicBSpline._changebasis_R(P1,P2) ≈ BasicBSpline._changebasis_I_new(P1,P2)

# case 1
k1 = knotvector" 1 1 11   1 1 1 1"
k2 = knotvector"1 1 1 1 3 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)


# works!
k1 = knotvector" 1 1 11 2 1 1 1 1"
k2 = knotvector"1 1 1 1 3 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I(P1,P2)

# works!
P1 = BSplineSpace{3}(KnotVector(1:8))
P2 = BSplineSpace{3}(KnotVector([2,3,3,4,5,6,7,8]))
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)

# works!
k3 = knotvector" 1 1 11   4"
k4 = knotvector"1 1 1 1 3 4"
P3 = BSplineSpace{3}(k3)
P4 = BSplineSpace{3}(k4)
BasicBSpline._changebasis_I_old(P3,P4)
BasicBSpline._changebasis_I_new(P3,P4)

# works!
k1 = knotvector" 1 1 11   1 1 1 1"
k2 = knotvector"1 1 1 1 3 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)

# works!
k1 = knotvector"33"
k2 = knotvector"33"
P1 = BSplineSpace{2}(k1)
P2 = BSplineSpace{2}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)

# works!
k1 = knotvector"213"
k2 = knotvector" 33"
P1 = BSplineSpace{2}(k1)
P2 = BSplineSpace{2}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)

# works!
k1 = knotvector"1 1 111111"
k2 = knotvector" 1 1 21111"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)



# case 1
k1 = knotvector" 1 1 11   1 1 1 1"
k2 = knotvector"1 1 1 1 2 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
P1 ⊆ P2
P1 ⊑ P2
BasicBSpline._changebasis_I_old(P1,P2)
BasicBSpline._changebasis_I_new(P1,P2)



plot(P1)
plot!(P2)

k3 = knotvector" 1 1 11   4"
k4 = knotvector"1 1 1 1 3 4"
P3 = BSplineSpace{3}(k3)
P4 = BSplineSpace{3}(k4)
BasicBSpline._changebasis_I_old(P3,P4)
BasicBSpline._changebasis_I_new(P3,P4)


i = 1
findfirst(!iszero, Aᵖ[i, :])
findlast(!iszero, Aᵖ[i, :])




k1 = KnotVector(rand(8)) + 4*KnotVector([0,1])
k1 = knotvector"1111 2 1111"
P1 = BSplineSpace{3}(k1)
P2 = expandspace_R(P1, Val(1), KnotVector(4*rand(3)))

P1 ⊆ P2
isdegenerate(P1)
isdegenerate(P2)

BasicBSpline._find_j_begin_end_R(P1, P2, 1, 1, 1)
BasicBSpline._find_j_begin_end_R(P1, P2, 2, 1, 1)
Matrix(changebasis(P1,P2))


##### impl
P = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
P′ = BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
P = P1
P′ = P2
i = 1
j_begin = 1
j_end = 1
p = degree(P)
p′ = degree(P′)

n′ = dim(P′)
k = knotvector(P)
k′ = knotvector(P′)

_, j_end = BasicBSpline._find_j_begin_end_R(P, P′, i, 1, 1)

t_begin = k[i]
m = p′-p
for ii in 0:p+1
    if k[i+ii] == t_begin
        m += 1
    else
        break
    end
end
m

k′[j_end+p′+1:-1:1]

for j in j_end+p′+1:-1:1
    if k′[j] == t_begin
        @show j
        m -= 1
    end
    if m == 0
        j_begin = j
        @show j_begin
        break
    end
end

BasicBSpline._find_j_begin_end_R(P, P′, i, 1, 1)


k[i:i+p+1]




BSplineSpace{2}(KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))





