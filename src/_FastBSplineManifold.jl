# Faster B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold{T} <: AbstractBSplineManifold
    bsplinespaces::Vector{<:FastBSplineSpace}
    controlpoints::Array{T} where T<:Point
    function FastBSplineManifold(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{T}) where T
        Ps = FastBSplineSpace.(Ps)
        if collect(size(a)) ≠ dim.(Ps)
            throw(DimensionMismatch())
        else
            P = convert(Vector{FastBSplineSpace}, Ps)
            new{T}(P, float(a))
        end
    end
end

struct BSplineCurve{p1,T} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    controlpoints::Array{T,1} where T<:Point
    function BSplineCurve(P1::AbstractBSplineSpace, a::AbstractArray{T,1}) where T
        p1 = degree(P1)
        if size(a) ≠ (dim(P1),)
            throw(DimensionMismatch())
        else
            new{p1,T}(P1, float(a))
        end
    end
end
function BSplineCurve(P::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Point,1})
    if length(P) ≠ 1
        throw(DimensionMismatch())
    else
        BSplineCurve(P[1], a)
    end
end

struct BSplineSurface{p1,p2,T} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    bsplinespace2::FastBSplineSpace{p2}
    controlpoints::Array{T,2} where T<:Point
    function BSplineSurface(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, a::AbstractArray{T,2}) where T
        p1, p2 = degree(P1), degree(P2)
        if size(a) ≠ (dim(P1), dim(P2))
            throw(DimensionMismatch())
        else
            new{p1,p2,T}(P1, P2, float(a))
        end
    end
end
function BSplineSurface(P::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Point,2})
    if length(P) ≠ 2
        throw(DimensionMismatch())
    else
        BSplineSurface(P[1], P[2], a)
    end
end

struct BSplineSolid{p1,p2,p3,T} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    bsplinespace2::FastBSplineSpace{p2}
    bsplinespace3::FastBSplineSpace{p3}
    controlpoints::Array{T,3} where T<:Point
    function BSplineSolid(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace, a::AbstractArray{T,3}) where T
        p1, p2, p3 = degree(P1), degree(P2), degree(P3)
        if size(a) ≠ (dim(P1), dim(P2), dim(P3))
            throw(DimensionMismatch())
        else
            new{p1,p2,p3,T}(P1, P2, P3, float(a))
        end
    end
end
function BSplineSolid(P::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Point,3})
    if length(P) ≠ 3
        throw(DimensionMismatch())
    else
        BSplineSolid(P[1], P[2], P[3], a)
    end
end

function BSplineManifold(M::BSplineCurve{p1}) where {p1}
    P1 = M.bsplinespace1
    a = M.controlpoints
    return BSplineManifold([P1], a)
end

function BSplineManifold(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    k1 = M.bsplinespace1
    P2 = M.bsplinespace2
    a = M.controlpoints
    return BSplineManifold([P1, P2], a)
end

function BSplineManifold(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    P1 = M.bsplinespace1
    P2 = M.bsplinespace2
    P3 = M.bsplinespace3
    a = M.controlpoints
    return BSplineManifold([P1, P2, P3], a)
end

"""
convert AbstractBSplineManifold to FastBSplineManifold
"""
function FastBSplineManifold(M::AbstractBSplineManifold)
    return FastBSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

function (M::FastBSplineManifold)(t::AbstractVector{<:Real})
    # TODO: faster
    return BSplineManifold(M)(t)
end

function (M::BSplineCurve{p1})(t1::Real) where {p1}
    P1, = bsplinespaces(M)
    a = controlpoints(M)
    n1 = dim(P1)
    k1 = P1.knots
    j1 = _knotindex(P1, t1)
    b1 = _bsplinebasis(P1,t1,j1)

    f1 = j1-p1
    l1 = j1
    S1 = b1[1]*a[f1]
    for i1 in f1+1:l1
        S1 += b1[i1-f1+1]*a[i1]
    end

    return S1
end

function (M::BSplineSurface{p1,p2})(t1::Real,t2::Real) where {p1} where {p2}
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    n1, n2 = dim(P1), dim(P2)
    k1, k2 = P1.knots, P2.knots
    j1, j2 = _knotindex(P1, t1), _knotindex(P2, t2)
    b1 = _bsplinebasis(P1,t1,j1)
    b2 = _bsplinebasis(P2,t2,j2)

    f1, f2 = j1-p1, j2-p2
    l1, l2 = j1, j2
    S2 = b2[1]*a[f1,f2]
    for i2 in f2+1:l2
        S2 += b2[i2-f2+1]*a[f1,i2]
    end
    S1 = b1[1]*S2
    for i1 in f1+1:l1
        S2 = b2[1]*a[i1,f2]
        for i2 in f2+1:l2
            S2 += b2[i2-f2+1]*a[i1,i2]
        end
        S1 += b1[i1-f1+1]*S2
    end

    return S1
end

function (M::BSplineSolid{p1,p2,p3})(t1::Real,t2::Real,t3::Real) where {p1} where {p2} where {p3}
    P1, P2, P3 = M.bsplinespace1, M.bsplinespace2, M.bsplinespace3
    a = controlpoints(M)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    k1, k2, k3 = P1.knots, P2.knots, P3.knots
    j1, j2, j3 = _knotindex(P1, t1), _knotindex(P2, t2), _knotindex(P3, t3)
    b1 = _bsplinebasis(P1,t1,j1)
    b2 = _bsplinebasis(P2,t2,j2)
    b3 = _bsplinebasis(P3,t3,j3)

    f1, f2, f3 = j1-p1, j2-p2, j3-p3
    l1, l2, l3 = j1, j2, j3
    S3 = b3[1]*a[f1,f2,f3]
    for i3 in f3+1:l3
        S3 += b3[i3-f3+1]*a[f1,f2,i3]
    end
    S2 = b1[1]*b2[1]*S3
    for i2 in f2+1:l2
        S3 = b3[1]*a[f1,i2,f3]
        for i3 in f3+1:l3
            S3 += b3[i3-f3+1]*a[f1,i2,i3]
        end
        S2 += b1[1]*b2[i2-f2+1]*S3
    end
    S1 = S2
    for i1 in f1+1:l1
        S3 = b3[1]*a[i1,f2,f3]
        for i3 in f3+1:l3
            S3 += b3[i3-f3+1]*a[i1,f2,i3]
        end
        S2 = b1[i1-f1+1]*b2[1]*S3
        for i2 in f2+1:l2
            S3 = b3[1]*a[i1,i2,f3]
            for i3 in f3+1:l3
                S3 += b3[i3-f3+1]*a[i1,i2,i3]
            end
            S2 += b1[i1-f1+1]*b2[i2-f2+1]*S3
        end
        S1 += S2
    end

    return S1
end

# function (M::Union{BSplineCurve,BSplineSurface,BSplineSolid})(t::AbstractVector{<:Real})
#     return M(t...)
# end

function (M::BSplineCurve)(t::AbstractVector{<:Real})
    return M(t...)
end

function (M::BSplineSurface)(t::AbstractVector{<:Real})
    return M(t...)
end

function (M::BSplineSolid)(t::AbstractVector{<:Real})
    return M(t...)
end

dim(M::FastBSplineManifold) = length(M.bsplinespaces)
dim(M::BSplineCurve) = 1
dim(M::BSplineSurface) = 2
dim(M::BSplineSolid) = 3

function bsplinespaces(M::FastBSplineManifold)
    return M.bsplinespaces
end

function bsplinespaces(M::BSplineCurve{p1}) where {p1}
    return (M.bsplinespace1, )
end

function bsplinespaces(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    return (M.bsplinespace1, M.bsplinespace2)
end

function bsplinespaces(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    return (M.bsplinespace1, M.bsplinespace2, M.bsplinespace3)
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
