# Faster B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    bsplinespaces::Array{T,1} where {T<:FastBSplineSpace}
    controlpoints::Array{Float64}
    function FastBSplineManifold(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
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
    function BSplineCurve(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
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
    function BSplineCurve{q1}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1}
        return BSplineCurve(Ps, a)
    end
end

struct BSplineSurface{p1,p2} <: AbstractBSplineManifold
    knotvector1::Array{Float64,1}
    knotvector2::Array{Float64,1}
    controlpoints::Array{Float64,3}
    function BSplineSurface(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
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
    function BSplineSurface{q1,q2}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1} where {q2}
        return BSplineSurface(Ps, a)
    end
end

struct BSplineSolid{p1,p2,p3} <: AbstractBSplineManifold
    knotvector1::Array{Float64,1}
    knotvector2::Array{Float64,1}
    knotvector3::Array{Float64,1}
    controlpoints::Array{Float64,4}
    function BSplineSolid(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
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
    function BSplineSolid{q1,q2,q3}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1} where {q2} where {q3}
        return BSplineSolid(Ps, a)
    end
end

for fname in (:BSplineCurve, :BSplineSurface, :BSplineSolid, :FastBSplineManifold, :BSplineManifold)
    @eval function $fname(Ps::AbstractVector{<:AbstractBSplineSpace}, a::Array{Array{T,1}} where {T<:Real})
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
# function mapping(M::FastBSplineManifold, t::Array{<:Real,1})
#     # TODO: faster
#     return mapping(BSplineManifold(M), t)
# end

function mapping(M::BSplineCurve{p1}, t::Array{<:Real,1}) where {p1}
    P1, = bsplinespaces(M)
    a = controlpoints(M)
    n1 = dim(P1)
    t1, = t
    v1 = M.knotvector1
    j1 = _knotindex(v1, t1)
    b1 = _bsplinebasis(P1,t1,j1)
    d̂ = size(a)[end]

    S1 = Array{Float64}(undef, d̂)
    f1 = j1-p1
    l1 = j1
    for xx in 1:d̂
        S1[xx] = b1[1]*a[f1,xx]
        for i1 in f1+1:l1
            S1[xx] += b1[i1-f1+1]*a[i1,xx]
        end
    end

    return S1
end

function mapping(M::BSplineSurface{p1,p2}, t::Array{<:Real,1}) where {p1} where {p2}
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    n1, n2 = dim(P1), dim(P2)
    t1, t2 = t
    v1, v2 = M.knotvector1, M.knotvector2
    j1, j2 = _knotindex(v1, t1), _knotindex(v2, t2)
    b1 = _bsplinebasis(P1,t1,j1)
    b2 = _bsplinebasis(P2,t2,j2)
    d̂ = size(a)[end]

    S1 = Array{Float64}(undef, d̂)
    S2 = Array{Float64}(undef, d̂)
    f1, f2 = j1-p1, j2-p2
    l1, l2 = j1, j2
    for xx in 1:d̂
        S2[xx] = b2[1]*a[f1,f2,xx]
        for i2 in f2+1:l2
            S2[xx] += b2[i2-f2+1]*a[f1,i2,xx]
        end
        S1[xx] = b1[1]*S2[xx]
        for i1 in f1+1:l1
            S2[xx] = b2[1]*a[i1,f2,xx]
            for i2 in f2+1:l2
                S2[xx] += b2[i2-f2+1]*a[i1,i2,xx]
            end
            S1[xx] += b1[i1-f1+1]*S2[xx]
        end
    end

    return S1
end

function mapping(M::BSplineSolid{p1,p2,p3}, t::Array{<:Real,1}) where {p1} where {p2} where {p3}
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    t1, t2, t3 = t
    v1, v2, v3 = M.knotvector1, M.knotvector2, M.knotvector3
    j1, j2, j3 = _knotindex(v1, t1), _knotindex(v2, t2), _knotindex(v3, t3)
    b1 = _bsplinebasis(P1,t1,j1)
    b2 = _bsplinebasis(P2,t2,j2)
    b3 = _bsplinebasis(P3,t3,j3)
    d̂ = size(a)[end]

    S1 = Array{Float64}(undef, d̂)
    S2 = Array{Float64}(undef, d̂)
    S3 = Array{Float64}(undef, d̂)
    f1, f2, f3 = j1-p1, j2-p2, j3-p3
    l1, l2, l3 = j1, j2, j3
    for xx in 1:d̂
        S3[xx] = b3[1]*a[f1,f2,f3,xx]
        for i3 in f3+1:l3
            S3[xx] += b3[i3-f3+1]*a[f1,f2,i3,xx]
        end
        S2[xx] = b1[1]*b2[1]*S3[xx]
        for i2 in f2+1:l2
            S3[xx] = b3[1]*a[f1,i2,f3,xx]
            for i3 in f3+1:l3
                S3[xx] += b3[i3-f3+1]*a[f1,i2,i3,xx]
            end
            S2[xx] += b1[1]*b2[i2-f2+1]*S3[xx]
        end
        S1[xx] = S2[xx]
        for i1 in f1+1:l1
            S3[xx] = b3[1]*a[i1,f2,f3,xx]
            for i3 in f3+1:l3
                S3[xx] += b3[i3-f3+1]*a[i1,f2,i3,xx]
            end
            S2[xx] = b1[i1-f1+1]*b2[1]*S3[xx]
            for i2 in f2+1:l2
                S3[xx] = b3[1]*a[i1,i2,f3,xx]
                for i3 in f3+1:l3
                    S3[xx] += b3[i3-f3+1]*a[i1,i2,i3,xx]
                end
                S2[xx] += b1[i1-f1+1]*b2[i2-f2+1]*S3[xx]
            end
            S1[xx] += S2[xx]
        end
    end

    return S1
end

dim(M::FastBSplineManifold) = length(M.bsplinespaces)
dim(M::BSplineCurve) = 1
dim(M::BSplineSurface) = 2
dim(M::BSplineSolid) = 3

function bsplinespaces(M::FastBSplineManifold)
    return M.bsplinespaces
end

function bsplinespaces(M::BSplineCurve{p1}) where {p1}
    return (FastBSplineSpace(p1, Knots(M.knotvector1)), )
end

function bsplinespaces(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    return (FastBSplineSpace(p1, Knots(M.knotvector1)), FastBSplineSpace(p2, Knots(M.knotvector2)))
end

function bsplinespaces(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    return (FastBSplineSpace(p1, Knots(M.knotvector1)), FastBSplineSpace(p2, Knots(M.knotvector2)), FastBSplineSpace(p3, Knots(M.knotvector3)))
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
