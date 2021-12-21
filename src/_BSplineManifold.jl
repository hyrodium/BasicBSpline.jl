# B-spline manifold

abstract type AbstractBSplineManifold{Dim} end

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
    @warn "The method `bsplinebasis(I::CartesianIndex, Ps::AbstractVector{BSplineSpace}, t::AbstractVector{<:Real})` will be deprecated."
    d = length(Ps)
    Bs = prod(bsplinebasis(I[i], Ps[i], t[i]) for i in 1:d)
    return tensorprod(Bs)
end

function bsplinesupport(I::CartesianIndex, Ps::AbstractVector{BSplineSpace})
    @warn "The method `bsplinesupport(I::CartesianIndex, Ps::AbstractVector{BSplineSpace})` will be deprecated."
    d = length(Ps)
    return [bsplinesupport(I[i], Ps[i]) for i in 1:d]
end

@doc raw"""
Calculate the dimension of B-spline manifold.
"""
dim

dim(M::AbstractBSplineManifold) = length(M.bsplinespaces)
