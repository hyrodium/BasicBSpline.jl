# B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    bsplinespaces::Array{T,1} where {T<:FastBSplineSpace}
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
            new{p1}(k1, a)
        end
    end
    function BSplineCurve{q1}(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real}) where {q1}
        return BSplineCurve(Ps, a)
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
            new{p1,p2}(k1, k2, a)
        end
    end
    function BSplineSurface{q1,q2}(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real}) where {q1} where {q2}
        return BSplineSurface(Ps, a)
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
            new{p1,p2,p3}(k1, k2, k3, a)
        end
    end
    function BSplineSolid{q1,q2,q3}(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real}) where {q1} where {q2} where {q3}
        return BSplineSolid(Ps, a)
    end
end

for fname in (:BSplineCurve, :BSplineSurface, :BSplineSolid, :FastBSplineManifold, :BSplineManifold)
    @eval function $fname(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::Array{Array{T,1}} where {T<:Real})
        d̂ = length(a[1])
        A = reshape(transpose(hcat(reshape(a, prod(size(a)))...)), size(a)..., d̂)
        return $fname(Ps, A)
    end
end

function BSplineManifold(M::BSplineCurve{p1}) where {p1}
    k1 = Knots(M.knotvector1)
    P1 = BSplineSpace(p1, k1)
    a = M.controlpoints
    return BSplineManifold([P1], a)
end

function BSplineManifold(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    k1 = Knots(M.knotvector1)
    P1 = BSplineSpace(p1, k1)
    k2 = Knots(M.knotvector2)
    P2 = BSplineSpace(p2, k2)
    a = M.controlpoints
    return BSplineManifold([P1, P2], a)
end

function BSplineManifold(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    k1 = Knots(M.knotvector1)
    P1 = BSplineSpace(p1, k1)
    k2 = Knots(M.knotvector2)
    P2 = BSplineSpace(p2, k2)
    k3 = Knots(M.knotvector3)
    P3 = BSplineSpace(p3, k3)
    a = M.controlpoints
    return BSplineManifold([P1, P2, P3], a)
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
    return mapping(BSplineManifold(M), t)
end

function mapping(M::BSplineCurve, t::Array{<:Real,1})
    # TODO: faster
    return mapping(BSplineManifold(M), t)
end

function mapping(M::BSplineSurface, t::Array{<:Real,1})
    # TODO: faster
    return mapping(BSplineManifold(M), t)
end

function mapping(M::BSplineSolid, t::Array{<:Real,1})
    # TODO: faster
    return mapping(BSplineManifold(M), t)
end

dim(M::FastBSplineManifold) = length(M.bsplinespaces)
dim(M::BSplineCurve) = 1
dim(M::BSplineSurface) = 2
dim(M::BSplineSolid) = 3

function bsplinespaces(M::FastBSplineManifold)
    return M.bsplinespaces
end

function bsplinespaces(M::BSplineCurve{p1}) where {p1}
    return [FastBSplineSpace(p1, Knots(M.knotvector1))]
end

function bsplinespaces(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    return [FastBSplineSpace(p1, Knots(M.knotvector1)), FastBSplineSpace(p2, Knots(M.knotvector2))]
end

function bsplinespaces(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    return [FastBSplineSpace(p1, Knots(M.knotvector1)), FastBSplineSpace(p2, Knots(M.knotvector2)), FastBSplineSpace(p3, Knots(M.knotvector3))]
end

function controlpoints(M::FastBSplineManifold)
    return M.controlpoints
end

function controlpoints(M::BSplineCurve{p1}) where {p1}
    return M.controlpoints
end

function controlpoints(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    return M.controlpoints
end

function controlpoints(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    return M.controlpoints
end
