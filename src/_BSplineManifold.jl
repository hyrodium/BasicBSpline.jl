# B-spline manifold
function ⊗(X::Array{Float64},Y::Array{Float64})::Array{Float64}
    m = size(X)
    n = size(Y)
    reshape(reshape(X,length(X)) * reshape(Y,length(Y))', m..., n...)
end

function tensorprod(X::Array{T,1}) where T <: Array{Float64}
    n=length(X)
    # X[1] ⊗ … ⊗ X[n]
    @inbounds Y = X[1]
    for i ∈ 2:n
        @inbounds Y = Y ⊗ X[i]
    end
    return Y
end

struct BSplineManifold
    bsplinespaces::Array{BSplineSpace,1}
    controlpoints::Array{Float64}
    function BSplineManifold(bsplinespaces::Array{BSplineSpace,1}, controlpoints::Array{Float64})
        if collect(size(controlpoints)[1:end-1]) ≠ dim.(bsplinespaces)
            error("dimension does not match")
        else
            new(bsplinespaces, controlpoints)
        end
    end
end

# function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
#     if length(𝒫s)==length(t)==1
#         return BSplineBasis(𝒫s[1],t[1])
#     elseif length(𝒫s)==length(t)==2
#         return BSplineBasis(𝒫s[1],t[1])*BSplineBasis(𝒫s[2],t[2])'
#     else
#         error("dimension does not match")
#     end
# end

@doc raw"""
Multi-dimentional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
    d = length(t)
    Bs = [BSplineBasis(𝒫s[i],t[i]) for i ∈ 1:d]
    return tensorprod(Bs)
end

# function Mapping(M::BSplineManifold, t::Array{Float64,1})
#     𝒫s = M.bsplinespaces
#     𝒂 = M.controlpoints
#     d=length(𝒫s)
#     d̂=size(𝒂)[end]
#     return [sum(BSplineBasis(𝒫s,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
# end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function Mapping(M::BSplineManifold, t::Array{Float64,1})
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d = length(Ps)
    d̂ = size(𝒂)[end]
    return [sum(BSplineBasis(Ps,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
end
