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

# Set seed
Random.seed!(11)

# Generate randomized test cases
p = 5
chars = '0' .+ (0:p+1)
spaces = BSplineSpace[]
for _ in 1:1500
    s = join(rand(chars,12))
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

function new_issqsubset(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
    p₊ = p′ - p
    p₊ < 0 && return false

    k = knotvector(P)
    k′ = knotvector(P′)
    !(k[1+p] == k′[1+p′] < k′[end-p′] == k[end-p]) && return false

    l = length(k)
    l′ = length(k′)
    inner_knotvector = view(k, p+2:l-p-1)
    inner_knotvector′ = view(k′, p′+2:l′-p′-1)

    _P = BSplineSpace{p}(inner_knotvector)
    _P′ = BSplineSpace{p′}(inner_knotvector′)
    return _P ⊆ _P′
end

for range in ranges
    sps = unique(sorted_spaces[range])
    for P1 in sps
        for P2 in sps
            @test (P1 ⊑ P2) == new_issqsubset(P1, P2)
        end
    end
end



for P1 in spaces
    for P2 in spaces
        @test (P1 ⊑ P2) == new_issqsubset(P1, P2)
    end
end

