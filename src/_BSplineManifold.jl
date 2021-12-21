# B-spline manifold

# TODO add more type parameter Deg, T
abstract type AbstractBSplineManifold{Dim} end

# struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim,Deg,T}
struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim}
    bsplinespaces::S
    controlpoints::Array{T,Dim₊₁}
    function BSplineManifold(Ps::S,a::Array{T,Dim₊₁}) where {S<:Tuple,Dim₊₁,T<:Real}
        if !all(isa.(Ps,AbstractBSplineSpace))
            # TODO: update error message
            error("invalid")
        end
        if size(a)[1:Dim₊₁-1] != dim.(Ps)
            # TODO: update error message
            error("invalid")
        end
        d = length(Ps)
        p = degree.(Ps)
        new{d,p,T,S,d+1}(Ps,a)
    end
end

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
