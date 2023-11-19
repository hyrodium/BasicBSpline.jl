using Revise
using BasicBSpline
using IntervalSets
using Test
using Plots
import BasicBSpline._lower_I


sin(3)
println("hoge")

Pa1 = BSplineSpace{3}(knotvector"1 1 11  11 1 111")
Pa2 = BSplineSpace{3}(knotvector"  2 11 111 1 111")
Pa3 = BSplineSpace{3}(knotvector"  1111 111 1 111")

Pb1 = BSplineSpace{3}(knotvector"1111 1111")
Pb2 = BSplineSpace{3}(knotvector" 2 2 1111")
Pb3 = BSplineSpace{3}(knotvector"1111 1 21")
Pb4 = BSplineSpace{3}(knotvector" 2 2 4")
Pb5 = BSplineSpace{3}(knotvector"   4 4")
Pb6 = BSplineSpace{4}(knotvector"   5 5")
Pb7 = BSplineSpace{4}(knotvector"  1415")

Pc1 = BSplineSpace{2}(knotvector"3     3")
Pc2 = BSplineSpace{2}(knotvector"3  1 23")
Pc3 = BSplineSpace{2}(knotvector"3 12 33")
Pc4 = BSplineSpace{3}(knotvector"4 23 44")
Pc5 = BSplineSpace{4}(knotvector"5 36 55")

Pa1 ⊑ Pa1
Pa1 ⊑ Pa2
Pa1 ⊑ Pa3

Pb1 ⊑ Pb1
Pb1 ⊑ Pb2
Pb1 ⊑ Pb3
Pb1 ⊑ Pb4
Pb1 ⊑ Pb5
Pb2 ⊑ Pb1
Pb2 ⊑ Pb2
Pb2 ⊑ Pb3
Pb2 ⊑ Pb4
Pb2 ⊑ Pb5
Pb3 ⊑ Pb1
Pb3 ⊑ Pb2
Pb3 ⊑ Pb3
Pb3 ⊑ Pb4
Pb3 ⊑ Pb5
Pb4 ⊑ Pb1
Pb4 ⊑ Pb2
Pb4 ⊑ Pb3
Pb4 ⊑ Pb4
Pb4 ⊑ Pb5
Pb5 ⊑ Pb1
Pb5 ⊑ Pb2
Pb5 ⊑ Pb3
Pb5 ⊑ Pb4
Pb5 ⊑ Pb5
Pb1 ⊑ Pb6
Pb2 ⊑ Pb6
Pb3 ⊑ Pb6
Pb4 ⊑ Pb6
Pb5 ⊑ Pb6
Pb1 ⊑ Pb7
Pb2 ⊑ Pb7
Pb3 ⊑ Pb7
Pb4 ⊑ Pb7
Pb5 ⊑ Pb7
Pb6 ⊑ Pb7

Pc1 ⊑ Pc2
Pc2 ⊑ Pc3
Pc3 ⊑ Pc4
Pc4 ⊑ Pc5

changebasis_I(Pc1, Pc2)


P = Pc1
P′= Pc2
function get_ranges(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
    n = dim(P)
    k = knotvector(P)
    k′ = knotvector(P′)
    ranges = fill(1:0, n)
    A = changebasis_I(P, P′)
    for i in 1:n
        isdegenerate_I(P, i) && continue
        flags = abs.(A[i, :]) .> 1e-13
        ranges[i] = findfirst(flags):findlast(flags)
        _k = k[i:i+p+1]
        @show _k
    end
    @show (p, p′)
    @show k′
    return ranges
end

get_ranges(Pc4, Pc5)
get_ranges(Pc5, Pc5)


Matrix(changebasis_I(Pc4, Pc5))
Matrix(changebasis_I(_lower_I(Pc4), _lower_I(Pc5)))

Matrix(changebasis_I(Pb4, Pb5))
Matrix(changebasis_I(_lower_I(Pb4), _lower_I(Pb5)))

Pb4
Pb5


Matrix(changebasis(Pc4, Pc5))




# 整数に限定してvisualization
# 


k = knotvector"10345"

function str(k::AbstractKnotVector)
    l = length(k)
    kᵢ = k[1]
    s = " "^(kᵢ-1)
    n = 0
    for i in 1:l
        if kᵢ ≠ k[i]
            s *= "$n"*" "^(k[i]-kᵢ-1)
            n = 1
        else
            n += 1
        end
        kᵢ = k[i]
    end
    s *= "$n"
    return s
end

str(k)



function show_ranges(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
    @assert P ⊑ P′
    p₊ = p′-p
    d = domain(P)
    n = dim(P)
    k = knotvector(P)
    k′ = knotvector(P′)
    ranges = fill(1:0, n)
    A = changebasis_I(P, P′)
    for i in 1:n
        isdegenerate_I(P, i) && continue
        flags = abs.(A[i, :]) .> 1e-13
        j_begin = findfirst(flags)
        j_end = findlast(flags)
        ranges[i] = j_begin:j_end
        
        kₛ = k[i:i+p+1]
        k₊ = clamp_knotvector(kₛ + unique(kₛ[[kᵢ in d for kᵢ in kₛ]]) * p₊, d)
        kᵣ = clamp_knotvector(k′,d)[j_begin:j_end+p′+1]

        if k₊ == kᵣ
            # println("kₛ   : ", str(kₛ))
            # println("k₊=kᵣ: ", str(k₊))
            # println("k′   : ", str(k′))
            # println()
        else
            println("i = $(i)")
            println("kₛ: ", str(kₛ))
            println("k₊: ", str(k₊))
            println("kᵣ: ", str(kᵣ))
            println("k′: ", str(k′))
            println()
        end
    end
    @show (p, p′)
    @show k′
    return ranges
end

function clamp_knotvector(k::AbstractKnotVector, d::ClosedInterval)
    _k = copy(KnotVector(k))
    l = length(k)
    v = _k.vector
    for i in 1:l
        v[i] = clamp(v[i],d)
    end
    return _k
end

clamp_knotvector(knotvector"111111", 2..4)

# show_ranges(Pc4, Pc5)
show_ranges(Pb4, Pb7)
changebasis(Pb4, Pb7)

# show_ranges(Pa1, Pa1)
# show_ranges(Pa1, Pa2)
show_ranges(Pa1, Pa3)
# show_ranges(Pb1, Pb1)
show_ranges(Pb1, Pb2)
show_ranges(Pb1, Pb3)

# show_ranges(Pb1, Pb4)
# show_ranges(Pb1, Pb5)
show_ranges(Pb2, Pb1)
# show_ranges(Pb2, Pb2)
show_ranges(Pb2, Pb3)
show_ranges(Pb2, Pb4)
show_ranges(Pb2, Pb5)
show_ranges(Pb3, Pb1)
show_ranges(Pb3, Pb2)
# show_ranges(Pb3, Pb3)
show_ranges(Pb3, Pb4)
show_ranges(Pb3, Pb5)
show_ranges(Pb4, Pb1)
show_ranges(Pb4, Pb2)
show_ranges(Pb4, Pb3)
# show_ranges(Pb4, Pb4)
show_ranges(Pb4, Pb5)
show_ranges(Pb5, Pb1)
show_ranges(Pb5, Pb2)
show_ranges(Pb5, Pb3)
show_ranges(Pb5, Pb4)
# show_ranges(Pb5, Pb5)
# show_ranges(Pb1, Pb6)
show_ranges(Pb2, Pb6)
show_ranges(Pb3, Pb6)
show_ranges(Pb4, Pb6)
# show_ranges(Pb5, Pb6)
show_ranges(Pb1, Pb7)
show_ranges(Pb2, Pb7)
show_ranges(Pb3, Pb7)
show_ranges(Pb4, Pb7)
# show_ranges(Pb5, Pb7)
# show_ranges(Pb6, Pb7)
# show_ranges(Pc1, Pc2)
# show_ranges(Pc2, Pc3)
# show_ranges(Pc3, Pc4)
# show_ranges(Pc4, Pc5)




show_ranges(Pa1, Pa3)
show_ranges(Pb1, Pb2)
show_ranges(Pb1, Pb3)
show_ranges(Pb2, Pb1)
show_ranges(Pb2, Pb3)
show_ranges(Pb2, Pb4)
show_ranges(Pb2, Pb5)
show_ranges(Pb3, Pb1)
show_ranges(Pb3, Pb2)
show_ranges(Pb3, Pb4)
show_ranges(Pb3, Pb5)
show_ranges(Pb4, Pb1)
show_ranges(Pb4, Pb2)
show_ranges(Pb4, Pb3)
show_ranges(Pb4, Pb5)
show_ranges(Pb5, Pb1)
show_ranges(Pb5, Pb2)
show_ranges(Pb5, Pb3)
show_ranges(Pb5, Pb4)
show_ranges(Pb2, Pb6)
show_ranges(Pb3, Pb6)
show_ranges(Pb4, Pb6)
show_ranges(Pb1, Pb7)
show_ranges(Pb2, Pb7)
show_ranges(Pb3, Pb7)
show_ranges(Pb4, Pb7)

plot(Pb4, xlims=extrema(domain(Pb4)))
plot!(Pb7)

domain(Pb4)
Pb7

knotvector(Pb4)
knotvector(Pb7)


Px1 = BSplineSpace{3}(knotvector"11111111")
Px2 = BSplineSpace{3}(knotvector"   41111")
changebasis_I(Px1, Px2)
changebasis_I(Px2, Px1)
show_ranges(Px1, Px2)
show_ranges(Px2, Px1)


Px1 = BSplineSpace{3}(knotvector"11111111")
Px2 = BSplineSpace{3}(knotvector"  131111")
changebasis_I(Px1, Px2)
changebasis_I(Px2, Px1)
show_ranges(Px1, Px2)
show_ranges(Px2, Px1)


Px1 = BSplineSpace{3}(knotvector"11111111")
Px2 = BSplineSpace{3}(knotvector" 2111111")
changebasis_I(Px1, Px2)
changebasis_I(Px2, Px1)
show_ranges(Px1, Px2)
show_ranges(Px2, Px1)


function similarities_on_boundaries(P1::BSplineSpace{p}, P2::BSplineSpace{p}) where p
    @assert domain(P1) == domain(P2)
    k1 = knotvector(P1)
    k2 = knotvector(P2)
    s1 = 0
    s2 = 0
    for i in reverse(1:p)
        if k1[1+i] == k2[1+i]
            s1 += 1
        else
            break
        end
    end
    for i in reverse(1:p)
        if k1[end-i] == k2[end-i]
            s2 += 1
        else
            break
        end
    end
    return (s1,s2)
end

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector"   4 4")
similarities_on_boundaries(Px1, Px2)
Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector"  31 112 ")
similarities_on_boundaries(Px1, Px2)
Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
similarities_on_boundaries(Px1, Px2)



Px1 = BSplineSpace{3}(knotvector" 121 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
changebasis_I(Px1, Px2)
similarities_on_boundaries(Px1, Px2)





function breadth_on_boundaries(P1::BSplineSpace{p}, P2::BSplineSpace{p}) where p
    @assert domain(P1) == domain(P2)
    k1 = knotvector(P1)
    k2 = knotvector(P2)
    s1 = 0
    s2 = 0
    for i in reverse(1:p)
        if k1[1+i] == k2[1+i]
            s1 += 1
        else
            break
        end
    end
    for i in reverse(1:p)
        if k1[end-i] == k2[end-i]
            s2 += 1
        else
            break
        end
    end
    return (p-s1+1, p-s2+1)
end

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector"   4 4")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector"  31 112 ")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)

Px1 = BSplineSpace{3}(knotvector" 121 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)





Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{4}(knotvector"   515")
changebasis_I(Px1, Px2)
# breadth_on_boundaries(Px1, Px2)  # (3, 3)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{4}(knotvector"2111 113")
changebasis_I(Px1, Px2)
# breadth_on_boundaries(Px1, Px2)  # (4,4)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{4}(knotvector"2111 311")
changebasis_I(Px1, Px2)
# breadth_on_boundaries(Px1, Px2)  # (4,4)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{4}(knotvector"2111 221")
changebasis_I(Px1, Px2)
# breadth_on_boundaries(Px1, Px2)  # (4,4)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector"  31 112 ")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)

Px1 = BSplineSpace{3}(knotvector"1111 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)

Px1 = BSplineSpace{3}(knotvector" 121 1111")
Px2 = BSplineSpace{3}(knotvector" 211 1111")
breadth_on_boundaries(Px1, Px2)
changebasis_I(Px1, Px2)




changebasis_I(Px1, Px1) * changebasis_I(Px1, Px1)





