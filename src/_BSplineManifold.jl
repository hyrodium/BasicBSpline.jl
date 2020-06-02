# B-spline manifold
abstract type AbstractBSplineManifold end

function ⊗(X::Array{Float64},Y::Array{Float64})::Array{Float64}
    m = size(X)
    n = size(Y)
    reshape(reshape(X,length(X)) * reshape(Y,length(Y))', m..., n...)
end

function tensorprod(X::Array{T,1}) where T <: Array{Float64}
    n=length(X)
    # X[1] ⊗ … ⊗ X[n]
    Y = X[1]
    for i ∈ 2:n
        Y = Y ⊗ X[i]
    end
    return Y
end

"""
B-spline manifold for general polynomial degree
"""
struct BSplineManifold <: AbstractBSplineManifold
    bsplinespaces::Array{BSplineSpace,1}
    controlpoints::Array{Float64}
    function BSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        Ps = BSplineSpace.(Ps)
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        else
            P = convert(Array{BSplineSpace,1}, Ps)
            a = convert(Array{Float64}, a)
            new(P, a)
        end
    end
end

function BSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::Array{Array{T,1}} where T<:Real)
    d̂ = length(a[1])
    A = reshape(transpose(hcat(reshape(a,prod(size(a)))...)), size(a)..., d̂)
    return BSplineManifold(Ps, A)
end


"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    bsplinespaces::Array{T,1} where T <: FastBSplineSpace
    controlpoints::Array{Float64}
    function FastBSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        Ps = FastBSplineSpace.(Ps)
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        else
            P = convert(Array{FastBSplineSpace,1}, Ps)
            a = convert(Array{Float64}, a)
            new(P, a)
        end
    end
end

function FastBSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::Array{Array{T,1}} where T<:Real)
    d̂ = length(a[1])
    A = reshape(transpose(hcat(reshape(a,prod(size(a)))...)), size(a)..., d̂)
    return FastBSplineManifold(Ps, A)
end


"""
convert FastBSplineManifold to BSplineManifold
"""
function BSplineManifold(M::AbstractBSplineManifold)
    BSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

"""
convert BSplineManifold to FastBSplineManifold
"""
function FastBSplineManifold(M::AbstractBSplineManifold)
    FastBSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

@doc raw"""
Multi-dimentional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function bsplinebasis(Ps::Array{BSplineSpace,1},t::Array{<:Real,1})
    d = length(t)
    Bs = [bsplinebasis(Ps[i],t[i]) for i ∈ 1:d]
    return tensorprod(Bs)
end

@doc raw"""
Multi-dimentional B-spline basis function.
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```
"""
function bsplinebasis(I::CartesianIndex, Ps::Array{BSplineSpace,1},t::Array{<:Real,1})
    d = length(Ps)
    Bs = prod(bsplinebasis(I[i],Ps[i],t[i]) for i ∈ 1:d)
    return tensorprod(Bs)
end

function bsplinesupport(I::CartesianIndex, Ps::Array{BSplineSpace,1})
    d = length(Ps)
    return [bsplinesupport(I[i],Ps[i]) for i ∈ 1:d]
end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function mapping(M::BSplineManifold, t::Array{<:Real,1})
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d = length(Ps)
    d̂ = size(𝒂)[end]
    return [sum(bsplinebasis(Ps,t).*𝒂[..,i]) for i ∈ 1:d̂]
end

function mapping(M::FastBSplineManifold, t::Array{<:Real,1})
    return mapping(BSplineManifold(M),t)
end


@doc raw"""
Calculate the dimention of B-spline manifold.
"""
function dim(M::AbstractBSplineManifold)
    length(M.bsplinespaces)
end
