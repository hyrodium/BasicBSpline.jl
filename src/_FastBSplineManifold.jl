# B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    # TODO: avoid abstract type FastBSplineSpace
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

struct BSplineCurve{p1} <: AbstractBSplineManifold
    knotvector1::Array{Float64,1}
    controlpoints::Array{Float64,2}
    function BSplineCurve(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        p1, k1 = degree(Ps[1]), knots(Ps[1]).vector
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        elseif p1 > MAX_DEGREE
            return BSplineManifold(Ps, a)
        else
            a = convert(Array{Float64}, a)
            new{p1}(k1,a)
        end
    end
end

struct BSplineSurface{p1,p2} <: AbstractBSplineManifold
    knotvector1::Array{Float64,1}
    knotvector2::Array{Float64,1}
    controlpoints::Array{Float64,3}
    function BSplineSurface(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        p1, k1 = degree(Ps[1]), knots(Ps[1]).vector
        p2, k2 = degree(Ps[2]), knots(Ps[2]).vector
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        elseif p1 > MAX_DEGREE || p2 > MAX_DEGREE
            return BSplineManifold(Ps, a)
        else
            a = convert(Array{Float64}, a)
            new{p1,p2}(k1,k2,a)
        end
    end
end

struct BSplineSolid{p1,p2,p3} <: AbstractBSplineManifold
    knotvector1::Array{Float64,1}
    knotvector2::Array{Float64,1}
    knotvector3::Array{Float64,1}
    controlpoints::Array{Float64,4}
    function BSplineSolid(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        p1, k1 = degree(Ps[1]), knots(Ps[1]).vector
        p2, k2 = degree(Ps[2]), knots(Ps[2]).vector
        p3, k3 = degree(Ps[3]), knots(Ps[3]).vector
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        elseif p1 > MAX_DEGREE || p2 > MAX_DEGREE || p3 > MAX_DEGREE
            return BSplineManifold(Ps, a)
        else
            a = convert(Array{Float64}, a)
            new{p1,p2,p3}(k1,k2,k3,a)
        end
    end
end

for fname in (:BSplineCurve, :BSplineSurface, :BSplineSolid, :FastBSplineManifold)
    @eval function $fname(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::Array{Array{T,1}} where T<:Real)
        d̂ = length(a[1])
        A = reshape(transpose(hcat(reshape(a,prod(size(a)))...)), size(a)..., d̂)
        return $fname(Ps, A)
    end
end

"""
convert AbstractBSplineManifold to FastBSplineManifold
"""
function FastBSplineManifold(M::AbstractBSplineManifold)
    return FastBSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function mapping(M::FastBSplineManifold, t::Array{<:Real,1})
    # TODO: faster
    return mapping(BSplineManifold(M),t)
end

@doc raw"""
Calculate the dimension of B-spline manifold.
"""
function dim(M::FastBSplineManifold)
    return length(M.bsplinespaces)
end

dim(M::BSplineCurve) = 1
dim(M::BSplineSurface) = 2
dim(M::BSplineSolid) = 3
