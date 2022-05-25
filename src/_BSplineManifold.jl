# B-spline manifold

abstract type AbstractManifold{Dim} end

# Broadcast like a scalar
Base.Broadcast.broadcastable(M::AbstractManifold) = Ref(M)

abstract type AbstractBSplineManifold{Dim, Deg} <: AbstractManifold{Dim} end

dim(::AbstractBSplineManifold{Dim}) where Dim = Dim

@doc raw"""
Construct Rational B-spline manifold from given control points and B-spline spaces.

# Examples
```jldoctest
julia> using StaticArrays

julia> P = BSplineSpace{2}(KnotVector(0,0,0,1,1,1))
BSplineSpace{2, Int64}(KnotVector([0, 0, 0, 1, 1, 1]))

julia> a = [SVector(1,0), SVector(1,1), SVector(0,1)]
3-element Vector{SVector{2, Int64}}:
 [1, 0]
 [1, 1]
 [0, 1]

julia> M = BSplineManifold(a,(P,));

julia> M(0.4)
2-element SVector{2, Float64} with indices SOneTo(2):
 0.84
 0.64
```
"""
struct BSplineManifold{Dim,Deg,C,S<:Tuple} <: AbstractBSplineManifold{Dim,Deg}
    bsplinespaces::S
    controlpoints::Array{C,Dim}
    function BSplineManifold(a::Array{C,Dim},Ps::S) where {S<:Tuple,C,Dim}
        for P in Ps
            if !(P isa AbstractBSplineSpace)
                throw(TypeError(:BSplineManifold,AbstractBSplineSpace,P))
            end
        end
        if size(a) != dim.(Ps)
            msg = "The size of control points array $(size(a)) and dimensions of B-spline spaces $(dim.(Ps)) must be equal."
            throw(DimensionMismatch(msg))
        end
        Deg = degree.(Ps)
        new{Dim,Deg,C,S}(Ps,a)
    end
end

Base.:(==)(M1::AbstractBSplineManifold, M2::AbstractBSplineManifold) = (bsplinespaces(M1)==bsplinespaces(M2)) & (controlpoints(M1)==controlpoints(M2))

bsplinespaces(M::BSplineManifold) = M.bsplinespaces
controlpoints(M::BSplineManifold) = M.controlpoints

@generated function unsafe_mapping(M::BSplineManifold{1,Deg},t1::Real) where {Deg}
    p1, = Deg
    exs = Expr[]
    for j1 in 1:p1
        push!(exs, :(v += b1[$(1+j1)]*getindex(a,i1+$(j1))))
    end
    Expr(:block,
        :((P1,) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(i1 = intervalindex(P1,t1)),
        :(b1 = bsplinebasisall(P1,i1,t1)),
        :(v = b1[1]*getindex(a,i1)),
        exs...,
        :(return v)
    )
end

@generated function unsafe_mapping(M::BSplineManifold{2,Deg},t1::Real,t2::Real) where {Deg}
    p1, p2 = Deg
    exs = Expr[]
    for j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(v += b1[$(j1)]*b2[$(j2)]*getindex(a,i1+$(j1-1),i2+$(j2-1))))
    end
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :((i1, i2) = (intervalindex(P1,t1), intervalindex(P2,t2))),
        :((b1, b2) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2))),
        :(v = b1[1]*b2[1]*getindex(a,i1,i2)),
        exs...,
        :(return v)
    )
end

@generated function unsafe_mapping(M::BSplineManifold{3,Deg},t1::Real,t2::Real,t3::Real) where {Deg}
    p1, p2, p3 = Deg
    exs = Expr[]
    for j3 in 1:p3+1, j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(v += b1[$(j1)]*b2[$(j2)]*b3[$(j3)]*getindex(a,i1+$(j1-1),i2+$(j2-1),i3+$(j3-1))))
    end
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2, P3) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :((i1, i2, i3) = (intervalindex(P1,t1), intervalindex(P2,t2), intervalindex(P3,t3))),
        :((b1, b2, b3) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2), bsplinebasisall(P3,i3,t3))),
        :(v = b1[1]*b2[1]*b3[1]*getindex(a,i1,i2,i3)),
        exs...,
        :(return v)
    )
end

@inline function (M::AbstractManifold{1})(t1::Real)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    unsafe_mapping(M,t1)
end

@inline function (M::AbstractManifold{2})(t1::Real,t2::Real)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    unsafe_mapping(M,t1,t2)
end

@inline function (M::AbstractManifold{3})(t1::Real,t2::Real,t3::Real)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    t3 in domain(Ps[3]) || throw(DomainError(t3, "The input $(t3) is out of range."))
    unsafe_mapping(M,t1,t2,t3)
end


## currying
# 1dim
@inline function (M::AbstractBSplineManifold{1})(::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end

# 2dim
@inline function (M::AbstractBSplineManifold{2})(::Colon,::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end
@inline function (M::AbstractBSplineManifold{2,p})(t1::Real,::Colon) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    b = sum(a[j1+i1,:]*B1[1+i1] for i1 in 0:p1)
    return BSplineManifold(b,(P2,))
end
@inline function (M::AbstractBSplineManifold{2,p})(::Colon,t2::Real) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[:,j2+i2]*B2[1+i2] for i2 in 0:p2)
    return BSplineManifold(b,(P1,))
end

# 3dim
@inline function (M::AbstractBSplineManifold{3})(::Colon,::Colon,::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end
@inline function (M::AbstractBSplineManifold{3,p})(t1::Real,::Colon,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    b = sum(a[j1+i1,:,:]*B1[1+i1] for i1 in 0:p1)
    return BSplineManifold(b,(P2,P3))
end
@inline function (M::AbstractBSplineManifold{3,p})(::Colon,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[:,j2+i2,:]*B2[1+i2] for i2 in 0:p2)
    return BSplineManifold(b,(P1,P3))
end
@inline function (M::AbstractBSplineManifold{3,p})(::Colon,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j3 = intervalindex(P3,t3)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[:,:,j3+i3]*B3[1+i3] for i3 in 0:p3)
    return BSplineManifold(b,(P1,P2))
end
@inline function (M::AbstractBSplineManifold{3,p})(t1::Real,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    j2 = intervalindex(P2,t2)
    B1 = bsplinebasisall(P1,j1,t1)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[j1+i1,j2+i2,:]*B1[1+i1]*B2[1+i2] for i1 in 0:p1, i2 in 0:p2)
    return BSplineManifold(b,(P3,))
end
@inline function (M::AbstractBSplineManifold{3,p})(t1::Real,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    j3 = intervalindex(P3,t3)
    B1 = bsplinebasisall(P1,j1,t1)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[j1+i1,:,j3+i3]*B1[1+i1]*B3[1+i3] for i1 in 0:p1, i3 in 0:p3)
    return BSplineManifold(b,(P2,))
end
@inline function (M::AbstractBSplineManifold{3,p})(::Colon,t2::Real,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    j3 = intervalindex(P3,t3)
    B2 = bsplinebasisall(P2,j2,t2)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[:,j2+i2,j3+i3]*B2[1+i2]*B3[1+i3] for i2 in 0:p2, i3 in 0:p3)
    return BSplineManifold(b,(P1,))
end

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro
