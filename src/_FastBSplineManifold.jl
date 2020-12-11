# Faster B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    bsplinespaces::Array{<:FastBSplineSpace,1}
    controlpoints::Array{Float64}
    function FastBSplineManifold(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
        Ps = FastBSplineSpace.(Ps)
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        else
            P = convert(Array{FastBSplineSpace,1}, Ps)
            a′ = convert(Array{Float64}, a)
            new(P, a′)
        end
    end
end

struct BSplineCurve{p1} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    controlpoints::Array{Float64,2}
    function BSplineCurve(P1::AbstractBSplineSpace, a::AbstractArray{<:Real})
        p1 = degree(P1)
        if size(a)[1:end-1] ≠ (dim(P1),)
            error("dimension does not match")
        else
            a′ = convert(Array{Float64}, a)
            new{p1}(P1, a′)
        end
    end
end
function BSplineCurve(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
    if collect(size(a)[1:end-1]) ≠ dim.(Ps)
        error("dimension does not match")
    elseif length(Ps) ≠ 1
        error("dimension does not match")
    else
        a′ = convert(Array{Float64}, a)
        return BSplineCurve(Ps[1],a′)
    end
end
function BSplineCurve{q1}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1}
    return BSplineCurve(Ps, a)
end

struct BSplineSurface{p1,p2} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    bsplinespace2::FastBSplineSpace{p2}
    controlpoints::Array{Float64,3}
    function BSplineSurface(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, a::AbstractArray{<:Real})
        p1, p2 = degree(P1), degree(P2)
        if size(a)[1:end-1] ≠ (dim(P1), dim(P2))
            error("dimension does not match")
        else
            a′ = convert(Array{Float64}, a)
            new{p1,p2}(P1, P2, a′)
        end
    end
end
function BSplineSurface(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
    if collect(size(a)[1:end-1]) ≠ dim.(Ps)
        error("dimension does not match")
    elseif length(Ps) ≠ 2
        error("dimension does not match")
    else
        a′ = convert(Array{Float64}, a)
        return BSplineSurface(Ps[1], Ps[2], a′)
    end
end
function BSplineSurface{q1,q2}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1} where {q2}
    return BSplineSurface(Ps, a)
end

struct BSplineSolid{p1,p2,p3} <: AbstractBSplineManifold
    bsplinespace1::FastBSplineSpace{p1}
    bsplinespace2::FastBSplineSpace{p2}
    bsplinespace3::FastBSplineSpace{p3}
    controlpoints::Array{Float64,4}
    function BSplineSolid(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace, a::AbstractArray{<:Real})
        p1, p2, p3 = degree(P1), degree(P2), degree(P3)
        if size(a)[1:end-1] ≠ (dim(P1), dim(P2), dim(P3))
            error("dimension does not match")
        else
            a′ = convert(Array{Float64}, a)
            new{p1,p2,p3}(P1, P2, P3, a′)
        end
    end
end
function BSplineSolid(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real})
    if collect(size(a)[1:end-1]) ≠ dim.(Ps)
        error("dimension does not match")
    elseif length(Ps) ≠ 3
        error("dimension does not match")
    else
        a′ = convert(Array{Float64}, a)
        BSplineSolid(Ps[1],Ps[2],Ps[3],a′)
    end
end
function BSplineSolid{q1,q2,q3}(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:Real}) where {q1} where {q2} where {q3}
    return BSplineSolid(Ps, a)
end

for fname in (:BSplineCurve, :BSplineSurface, :BSplineSolid, :FastBSplineManifold, :BSplineManifold)
    @eval function $fname(Ps::AbstractVector{<:AbstractBSplineSpace}, a::AbstractArray{<:AbstractVector{<:Real}})
        d̂ = length(a[1])
        A = reshape(transpose(hcat(reshape(a, prod(size(a)))...)), size(a)..., d̂)
        return $fname(Ps, A)
    end
end

function BSplineManifold(M::BSplineCurve{p1}) where {p1}
    k1 = M.knots1
    P1 = BSplineSpace(p1, k1)
    a = M.controlpoints
    return BSplineManifold([P1], a)
end

function BSplineManifold(M::BSplineSurface{p1,p2}) where {p1} where {p2}
    k1 = M.knots1
    P1 = BSplineSpace(p1, k1)
    k2 = M.knots2
    P2 = BSplineSpace(p2, k2)
    a = M.controlpoints
    return BSplineManifold([P1, P2], a)
end

function BSplineManifold(M::BSplineSolid{p1,p2,p3}) where {p1} where {p2} where {p3}
    k1 = M.knots1
    P1 = BSplineSpace(p1, k1)
    k2 = M.knots2
    P2 = BSplineSpace(p2, k2)
    k3 = M.knots3
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

function (M::BSplineSurface{p1,p2})(t1::Real,t2::Real) where {p1} where {p2}
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    n1, n2 = dim(P1), dim(P2)
    k1, k2 = P1.knots, P2.knots
    j1, j2 = _knotindex(P1, t1), _knotindex(P2, t2)
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

function (M::BSplineSolid{p1,p2,p3})(t1::Real,t2::Real,t3::Real) where {p1} where {p2} where {p3}
    P1, P2, P3 = M.bsplinespace1, M.bsplinespace2, M.bsplinespace3
    a = controlpoints(M)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    k1, k2, k3 = P1.knots, P2.knots, P3.knots
    j1, j2, j3 = _knotindex(P1, t1), _knotindex(P2, t2), _knotindex(P3, t3)
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

function (M::Union{BSplineCurve,BSplineSurface,BSplineSolid})(t::AbstractVector{<:Real})
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
