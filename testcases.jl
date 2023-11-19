using BasicBSpline
using BenchmarkTools
using Plots
using Test

Q3 = BSplineSpace{4, Int64, KnotVector{Int64}}(KnotVector([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 7]))
Q4 = BSplineSpace{5, Int64, KnotVector{Int64}}(KnotVector([1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 5, 6]))
changebasis_I(Q3, Q4)
BasicBSpline.__changebasis_I_old(Q3, Q4)

changebasis_I(Q3, Q4)


changebasis_I(Q3, Q4) - BasicBSpline.__changebasis_I_old(Q3, Q4)











P1 = BSplineSpace{3}(knotvector"1111 132")
P2 = BSplineSpace{3}(knotvector"11121132")
isdegenerate_I(P1)
isdegenerate_I(P2)
isdegenerate_R(P1)
isdegenerate_R(P2)
P1 ⊆ P2
P1 ⊑ P2

A_R = changebasis_R(P1, P2)
A_I = changebasis_I(P1, P2)

@benchmark changebasis_R(P1, P2)
@benchmark changebasis_I(P1, P2)

A_R
for j in 1:dim(P2)
    @test any(Base.isstored.(Ref(A_R), 1:dim(P1), j)) == isnondegenerate_R(P2, j)
end




domain(P1)
domain(P2)

plot(P)

domain(P)

P1 

[3][4]

Q1 = BSplineSpace{1, Int64, KnotVector{Int64}}(KnotVector([2, 2, 4, 4, 6, 6]))
Q2 = BSplineSpace{3, Int64, KnotVector{Int64}}(KnotVector([1, 1, 1, 2, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7]))
A = changebasis_I(Q1, Q2)


Q3 = BSplineSpace{4, Int64, KnotVector{Int64}}(KnotVector([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 7]))
Q4 = BSplineSpace{5, Int64, KnotVector{Int64}}(KnotVector([1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 5, 6]))
changebasis_I(Q3, Q4)

changebasis_I(BasicBSpline._lower_I(Q3), BasicBSpline._lower_I(Q4))

exactdim_I(Q4)

isdegenerate_I(Q4)
plot(Q4)
domain(Q4)

BasicBSpline.__changebasis_I_old(Q3,Q4)
changebasis_I(Q3, Q4)





Q3 = BSplineSpace{4, Int64, KnotVector{Int64}}(KnotVector([1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 7]))
Q4 = BSplineSpace{5, Int64, KnotVector{Int64}}(KnotVector([1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 5, 6]))
changebasis_I(Q3, Q4)
BasicBSpline.__changebasis_I_old(Q3, Q4)


changebasis_I(Q3, Q4) - BasicBSpline.__changebasis_I_old(Q3, Q4)


function test_changebasis_I(P, P′; check_zero=true)
    # The test case must hold P ⊑ P′
    @test P ⊑ P′
    # Test type stability
    A = @inferred changebasis_I(P,P′)
    # Test output type
    @test A isa SparseMatrixCSC
    # Zeros must not be stored
    if check_zero
        @test !any(iszero.(A.nzval))
    end
    # Test the size of A
    n = dim(P)
    n′ = dim(P′)
    @test size(A) == (n,n′)
    # Elements must not be stored for degenerate row/col
    for j in 1:n′
        @test any(Base.isstored.(Ref(A), 1:n, j)) == isnondegenerate_I(P′, j)
    end
    for i in 1:n
        @test any(Base.isstored.(Ref(A), i, 1:n′)) == isnondegenerate_I(P, i)
    end
    # B_{(i,p,k)} = ∑ⱼ A_{i,j} B_{(j,p′,k′)}
    ts = range(domain(P), length=21)[2:end-1]
    for t in ts
        @test norm(bsplinebasis.(P,1:n,t) - A*bsplinebasis.(P′,1:n′,t), Inf) < ε
    end
end

test_changebasis_I(Q3, Q4)


ε = 1e-13
using LinearAlgebra
