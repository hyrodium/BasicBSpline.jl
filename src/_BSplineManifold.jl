# B-spline manifold

abstract type AbstractBSplineManifold end

function ⊗(X::AbstractArray{<:Real}, Y::AbstractArray{<:Real})
    # TODO: remove this function and use Tensor.jl
    m = size(X)
    n = size(Y)
    reshape(reshape(X, length(X)) * reshape(Y, length(Y))', m..., n...)
end

function tensorprod(X::AbstractArray)
    n = length(X)
    # X[1] ⊗ … ⊗ X[n]
    Y = X[1]
    for i in 2:n
        Y = Y ⊗ X[i]
    end
    return Y
end

"""
B-spline manifold for general polynomial degree
"""
struct BSplineManifold{T} <: AbstractBSplineManifold
    bsplinespaces::Vector{BSplineSpace}
    controlpoints::Array{T} where T<:Point
    function BSplineManifold(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
        Ps = BSplineSpace.(Ps)
        if collect(size(a)) ≠ dim.(Ps)
            throw(DimensionMismatch())
        else
            P = convert(Vector{BSplineSpace}, Ps)
            a′ = float(a)
            new{T}(P, a′)
        end
    end
end

"""
convert AbstractBSplineManifold to BSplineManifold
"""
function BSplineManifold(M::AbstractBSplineManifold)
    BSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

@doc raw"""
Multi-dimensional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function bsplinebasis(Ps::AbstractVector{BSplineSpace}, t::AbstractVector{<:Real})
    d = length(t)
    Bs = [bsplinebasis(Ps[i], t[i]) for i in 1:d]
    return tensorprod(Bs)
end

@doc raw"""
Multi-dimensional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function bsplinebasis(I::CartesianIndex, Ps::AbstractVector{BSplineSpace}, t::AbstractVector{<:Real})
    d = length(Ps)
    Bs = prod(bsplinebasis(I[i], Ps[i], t[i]) for i in 1:d)
    return tensorprod(Bs)
end

function bsplinesupport(I::CartesianIndex, Ps::AbstractVector{BSplineSpace})
    d = length(Ps)
    return [bsplinesupport(I[i], Ps[i]) for i in 1:d]
end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function (M::BSplineManifold)(t::AbstractVector{<:Real})
    Ps = M.bsplinespaces
    a = M.controlpoints
    d = length(Ps)
    d̂ = size(a)[end]
    N = prod(dim.(Ps))

    B = bsplinebasis(Ps, t)
    B_flat = reshape(B,N)
    a_flat = reshape(a,N,d̂)
    return [sum(B_flat .* a_flat[:,i]) for i in 1:d̂]
end

@doc raw"""
Calculate the dimension of B-spline manifold.
"""
dim

dim(M::AbstractBSplineManifold) = length(M.bsplinespaces)

function bsplinespaces(M::BSplineManifold)
    return M.bsplinespaces
end

function controlpoints(M::BSplineManifold)
    return M.controlpoints
end
