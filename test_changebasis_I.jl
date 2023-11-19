using BasicBSpline
using Random
using Plots
using Test
using SparseArrays
using LinearAlgebra
using IntervalSets

# Useful functions
function _isless(i1::ClosedInterval, i2::ClosedInterval)
    return extrema(i1) < extrema(i2)
end
function test_changebasis_I(P1,P2)
    ε = 1e-13
    @test P1 ⊑ P2
    A = @inferred changebasis_I(P1,P2)
    @test A isa SparseMatrixCSC
    # @test !any(iszero.(A.nzval))
    n1 = dim(P1)
    n2 = dim(P2)
    @test size(A) == (n1,n2)
    d = domain(P1)
    ts = range(extrema(d)..., length=21)[2:end-1]
    for t in ts
        @test norm(bsplinebasis.(P1,1:n1,t) - A*bsplinebasis.(P2,1:n2,t), Inf) < ε
    end
end

# Set seed
Random.seed!(11)

# Generate randomized test cases
p = 5
chars = '0' .+ (0:p+1)
spaces = BSplineSpace[]
for _ in 1:1000
    s = join(rand(chars,5))
    k = sum(KnotVector(findall(==('0'+i), s))*i for i in 1:9)
    P = BSplineSpace{rand(0:p)}(k)
    dim(P) < p+1 && continue
    width(domain(P)) == 0 && continue
    push!(spaces, P)
end
unique!(spaces)
domains = domain.(spaces)
perm = sortperm(domains, lt=_isless)
sorted_spaces = spaces[perm]
sorted_domains = domains[perm]
ranges = UnitRange{Int}[]
c = 1
for _ in 1:100
    r = searchsorted(sorted_domains, sorted_domains[c], lt=_isless)
    push!(ranges, r)
    c = maximum(r) + 1
    maximum(r) == length(spaces) && break
end
testcases = Tuple{BSplineSpace, BSplineSpace}[]
for range in ranges
    sps = unique(sorted_spaces[range])
    for P1 in sps
        for P2 in sps
            if P1 ⊑ P2
                push!(testcases, (P1, P2))
            end
        end
    end
end
ranges
testcases

unique(sorted_domains)

for (i, (P1, P2)) in enumerate(testcases)
    println(i)
    test_changebasis_I(P1, P2)
end

# i=1から片付ける
# 以下の2つに分類して後者を狭めていく
    # 計算できるケース
    # 計算できないケース

# Pick up for i=1
i = 1
j_ranges = Vector{Union{Nothing, UnitRange{Int}}}(nothing, length(testcases))
for (c, (P1, P2)) in enumerate(testcases)
    if isnondegenerate_I(P1, i)
        j_ranges[c] = BasicBSpline._find_j_range_I(P1,P2,i,1,dim(P2))
    end
end

c = 49
j_ranges[c]
P1, P2 = testcases[c]
plot(P1)
plot!(P2)
k1 = knotvector(P1)
k2 = knotvector(P2)
p1 = degree(P1)
p2 = degree(P2)
k1[i:i+p1+1]
k2[i:i+p2+1]
k1_sub = KnotVector(clamp.(collect(k1[i:i+p1+1]), Ref(domain(P1))))
j_begin = 1
j_end = 1
KnotVector(clamp.(collect(k1[j_begin:j_end+p1+1]), Ref(domain(P1))))

j_end = findfirst(j_end -> k1_sub ⊆ KnotVector(clamp.(collect(k2[1:j_end+p1+1]), Ref(domain(P1)))), 1:dim(P2))
j_begin = findlast(j_begin -> k1_sub ⊆ KnotVector(clamp.(collect(k2[j_begin:j_end+p2+1]), Ref(domain(P1)))), 1:dim(P2))
j_begin:j_end
j_ranges[c]

for c in 1:length(testcases)
    isnothing(j_ranges[c]) && continue
    println(c, ",", j_ranges[c])
    P1, P2 = testcases[c]
    plot(P1)
    plot!(P2)
    k1 = knotvector(P1)
    k2 = knotvector(P2)
    p1 = degree(P1)
    p2 = degree(P2)
    k1[i:i+p1+1]
    k2[i:i+p2+1]
    k1_sub = KnotVector(clamp.(collect(k1[i:i+p1+1]), Ref(domain(P1))))
    j_begin = 1
    j_end = 1
    KnotVector(clamp.(collect(k1[j_begin:j_end+p1+1]), Ref(domain(P1))))

    j_end = findfirst(j_end -> k1_sub ⊆ KnotVector(clamp.(collect(k2[1:j_end+p1+1]), Ref(domain(P1)))), 1:dim(P2))
    j_begin = findlast(j_begin -> k1_sub ⊆ KnotVector(clamp.(collect(k2[j_begin:j_end+p2+1]), Ref(domain(P1)))), 1:dim(P2))
    @test j_begin:j_end == j_ranges[c]
end

P1, P2 = testcases[661]
Matrix(changebasis_I(P1, P2))

plot(P1)
plot!(P2)
domain(P1)

# Pick up for i=2
i = 2
j_ranges = Vector{Union{Nothing, UnitRange{Int}}}(nothing, length(validcases))
for (c, (P1, P2)) in enumerate(validcases)
    if isnondegenerate_I(P1, i)
        j_ranges[c] = BasicBSpline._find_j_range_I(P1,P2,i,1,dim(P2))
    end
end

c = 2
j_ranges[c]
P1, P2 = validcases[c]
plot(P1)
plot!(P2)
k1 = knotvector(P1)
k2 = knotvector(P2)
p1 = degree(P1)
p2 = degree(P2)
k1[i:i+p1+1]
k2[1:2+p2+1]
j_begin = 1
j_end = 2
k1_sub = KnotVector(clamp.(collect(k1[i:i+p1+1]), Ref(domain(P1))))
KnotVector(clamp.(collect(k1[j_begin:j_end+p1+1]), Ref(domain(P1))))
k1[i:i+p1+1]
k2[j_begin:j_end+p2+1]

j_end = findfirst(j_end -> k1_sub ⊆ KnotVector(clamp.(collect(k2[1:j_end+p1+1]), Ref(domain(P1)))), 1:dim(P2))
j_begin = findlast(j_begin -> k1_sub ⊆ KnotVector(clamp.(collect(k2[j_begin:j_end+p2+1]), Ref(domain(P1)))), 1:dim(P2))
j_begin:j_end
j_ranges[c]








collect(k1[i:i+p1+1])


clamp(3,1..5)









validcases
c = 1
j_ranges[c]
P1, P2 = validcases[c]
plot(P1)
plot!(P2)
k1 = knotvector(P1)
k2 = knotvector(P2)
p1 = degree(P1)
p2 = degree(P2)
k1[i:i+p+1]
k2[1:5]



A = changebasis_I(P1,P2)
n1 = dim(P1)
n2 = dim(P2)
d = domain(P1)
ts = range(extrema(d)..., length=20)
for t in ts
    println(norm(bsplinebasis.(P1,1:n1,t) - A*bsplinebasis.(P2,1:n2,t), Inf))
end


t = 3.9999
bsplinebasis.(P1,1:n1,t)
A*bsplinebasis.(P2,1:n2,t)

plot(P2)

P2

j_ranges[c]

P2



invalidcases
testcases

P1, P2 = invalidcases[1]
domain(P1)

plot(P1)
plot!(P2)
isdegenerate_I(P1)
isdegenerate_I(P2)

BasicBSpline._changebasis_I_old(P1,P2)
isdegenerate_I.(P2,1:dim(P2))

BasicBSpline._changebasis_I(P1,P2)




# 現状で計算できるケース
    # i=1から片付ける
    # 以下の2つに分類して後者を狭めていく
        # 計算できるケース
        # 計算できないケース
# エラーになるケース
    # とりあえずcasesを列挙して集めとく
    # 


s1 = join(rand(chars,5))
s2 = join(rand(chars,5))
k1 = sum(KnotVector(findall(==('0'+i), s1))*i for i in 1:9)
k2 = sum(KnotVector(findall(==('0'+i), s2))*i for i in 1:9)
P1 = BSplineSpace{p}(k1)
P2 = BSplineSpace{p}(k2)
dim(P1) < 1 && continue
dim(P2) < 1 && continue
P1 ⊑ P2 && push!(testcases, (P1, P2))
P2 ⊑ P1 && push!(testcases, (P2, P1))

testcases

testcases

P1, P2 = testcases[11]

P1
P2

i = 25
testcases[i][1]
testcases[i][2]

dim(P2)

plot(P1)
plot!(P2)

changebasis_I(P1, P2)

changebasis_I(P1, P2)

changebasis_I(P2, P2)

P1
P2


