# Rational B-spline manifold (NURBS)
abstract type AbstractRationalBSplineManifold{Dim, Deg} <: AbstractManifold{Dim} end

"""
Construct Rational B-spline manifold from given control points, weights and B-spline spaces.

# Examples
```jldoctest
julia> using StaticArrays, LinearAlgebra

julia> P = BSplineSpace{2}(KnotVector([0,0,0,1,1,1]))
BSplineSpace{2, Int64}(KnotVector([0, 0, 0, 1, 1, 1]))

julia> w = [1, 1/√2, 1]
3-element Vector{Float64}:
 1.0
 0.7071067811865475
 1.0

julia> a = [SVector(1,0), SVector(1,1), SVector(0,1)]
3-element Vector{SVector{2, Int64}}:
 [1, 0]
 [1, 1]
 [0, 1]

julia> M = RationalBSplineManifold(a,w,P);  # 1/4 arc

julia> M(0.3)
2-element SVector{2, Float64} with indices SOneTo(2):
 0.8973756499953727
 0.4412674277525845

julia> norm(M(0.3))
1.0
```
"""
struct RationalBSplineManifold{Dim,Deg,C,S<:NTuple{Dim, AbstractBSplineSpace},T} <: AbstractRationalBSplineManifold{Dim,Deg}
    bsplinespaces::S
    controlpoints::Array{C,Dim}
    weights::Array{T,Dim}
    function RationalBSplineManifold(a::Array{C,Dim},w::Array{T,Dim},Ps::S) where {S<:NTuple{Dim, AbstractBSplineSpace},C,T<:Real} where Dim
        if size(a) != dim.(Ps)
            msg = "The size of control points array $(size(a)) and dimensions of B-spline spaces $(dim.(Ps)) must be equal."
            throw(DimensionMismatch(msg))
        end
        Deg = degree.(Ps)
        new{Dim,Deg,C,S,T}(Ps,a,w)
    end
end

RationalBSplineManifold(a::Array{C,Dim},w::Array{T,Dim},Ps::Vararg{AbstractBSplineSpace, Dim}) where {C,Dim,T<:Real} = RationalBSplineManifold(a,w,Ps)

Base.:(==)(M1::AbstractRationalBSplineManifold, M2::AbstractRationalBSplineManifold) = (bsplinespaces(M1)==bsplinespaces(M2)) & (controlpoints(M1)==controlpoints(M2)) & (weights(M1)==weights(M2))

controlpoints(M::RationalBSplineManifold) = M.controlpoints
weights(M::RationalBSplineManifold) = M.weights
bsplinespaces(M::RationalBSplineManifold) = M.bsplinespaces

@generated function unbounded_mapping(M::RationalBSplineManifold{1,Deg},t1::Real) where {Deg}
    p1, = Deg
    exs = Expr[]
    for j1 in 1:p1
        push!(exs, :(u += b1[$(1+j1)]*getindex(w,i1+$(j1))))
        push!(exs, :(v += b1[$(1+j1)]*getindex(w,i1+$(j1))*getindex(a,i1+$(j1))))
    end
    Expr(:block,
        :((P1,) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(w = weights(M)),
        :(i1 = intervalindex(P1,t1)),
        :(b1 = bsplinebasisall(P1,i1,t1)),
        :(u = b1[1]*getindex(w,i1)),
        :(v = b1[1]*getindex(w,i1)*getindex(a,i1)),
        exs...,
        :(return v/u)
    )
end

@generated function unbounded_mapping(M::RationalBSplineManifold{2,Deg},t1::Real,t2::Real) where {Deg}
    p1, p2 = Deg
    exs = Expr[]
    for j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(u += b1[$(j1)]*b2[$(j2)]*getindex(w,i1+$(j1-1),i2+$(j2-1))))
        push!(exs, :(v += b1[$(j1)]*b2[$(j2)]*getindex(w,i1+$(j1-1),i2+$(j2-1))*getindex(a,i1+$(j1-1),i2+$(j2-1))))
    end
    deleteat!(exs,1)
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(w = weights(M)),
        :((i1, i2) = (intervalindex(P1,t1), intervalindex(P2,t2))),
        :((b1, b2) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2))),
        :(u = b1[1]*b2[1]*getindex(w,i1,i2)),
        :(v = b1[1]*b2[1]*getindex(w,i1,i2)*getindex(a,i1,i2)),
        exs...,
        :(return v/u)
    )
end

@generated function unbounded_mapping(M::RationalBSplineManifold{3,Deg},t1::Real,t2::Real,t3::Real) where {Deg}
    p1, p2, p3 = Deg
    exs = Expr[]
    for j3 in 1:p3+1, j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(u += b1[$(j1)]*b2[$(j2)]*b3[$(j3)]*getindex(w,i1+$(j1-1),i2+$(j2-1),i3+$(j3-1))))
        push!(exs, :(v += b1[$(j1)]*b2[$(j2)]*b3[$(j3)]*getindex(w,i1+$(j1-1),i2+$(j2-1),i3+$(j3-1))*getindex(a,i1+$(j1-1),i2+$(j2-1),i3+$(j3-1))))
    end
    deleteat!(exs,1)
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2, P3) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(w = weights(M)),
        :((i1, i2, i3) = (intervalindex(P1,t1), intervalindex(P2,t2), intervalindex(P3,t3))),
        :((b1, b2, b3) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2), bsplinebasisall(P3,i3,t3))),
        :(u = b1[1]*b2[1]*b3[1]*getindex(w,i1,i2,i3)),
        :(v = b1[1]*b2[1]*b3[1]*getindex(w,i1,i2,i3)*getindex(a,i1,i2,i3)),
        exs...,
        :(return v/u)
    )
end


## currying
# 1dim
@inline function (M::AbstractRationalBSplineManifold{1})(::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end

# 2dim
@inline function (M::AbstractRationalBSplineManifold{2})(::Colon,::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end
@inline function (M::AbstractRationalBSplineManifold{2,p})(t1::Real,::Colon) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    w′ = sum(w[j1+i1,:]*B1[1+i1] for i1 in 0:p1)
    a′ = sum(a[j1+i1,:].*w[j1+i1,:].*B1[1+i1] for i1 in 0:p1) ./ w′
    return RationalBSplineManifold(a′,w′,(P2,))
end
@inline function (M::AbstractRationalBSplineManifold{2,p})(::Colon,t2::Real) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    w′ = sum(w[:,j2+i2]*B2[1+i2] for i2 in 0:p2)
    a′ = sum(a[:,j2+i2].*w[:,j2+i2].*B2[1+i2] for i2 in 0:p2) ./ w′
    return RationalBSplineManifold(a′,w′,(P1,))
end

# 3dim
@inline function (M::AbstractRationalBSplineManifold{3})(::Colon,::Colon,::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(t1::Real,::Colon,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    w′ = sum(w[j1+i1,:,:]*B1[1+i1] for i1 in 0:p1)
    a′ = sum(a[j1+i1,:,:].*w[j1+i1,:,:].*B1[1+i1] for i1 in 0:p1) ./ w′
    return RationalBSplineManifold(a′,w′,(P2,P3))
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(::Colon,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    w′ = sum(w[:,j2+i2,:]*B2[1+i2] for i2 in 0:p2)
    a′ = sum(a[:,j2+i2,:].*w[:,j2+i2,:].*B2[1+i2] for i2 in 0:p2) ./ w′
    return RationalBSplineManifold(a′,w′,(P1,P3))
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(::Colon,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j3 = intervalindex(P3,t3)
    B3 = bsplinebasisall(P3,j3,t3)
    w′ = sum(w[:,:,j3+i3]*B3[1+i3] for i3 in 0:p3)
    a′ = sum(a[:,:,j3+i3].*w[:,:,j3+i3].*B3[1+i3] for i3 in 0:p3) ./ w′
    return RationalBSplineManifold(a′,w′,(P1,P2))
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(t1::Real,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j1 = intervalindex(P1,t1)
    j2 = intervalindex(P2,t2)
    B1 = bsplinebasisall(P1,j1,t1)
    B2 = bsplinebasisall(P2,j2,t2)
    w′ = sum(w[j1+i1,j2+i2,:]*B1[1+i1]*B2[1+i2] for i1 in 0:p1, i2 in 0:p2)
    a′ = sum(a[j1+i1,j2+i2,:].*w[j1+i1,j2+i2,:].*B1[1+i1].*B2[1+i2] for i1 in 0:p1, i2 in 0:p2) ./ w′
    return RationalBSplineManifold(a′,w′,(P3,))
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(t1::Real,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j1 = intervalindex(P1,t1)
    j3 = intervalindex(P3,t3)
    B1 = bsplinebasisall(P1,j1,t1)
    B3 = bsplinebasisall(P3,j3,t3)
    w′ = sum(w[j1+i1,:,j3+i3]*B1[1+i1]*B3[1+i3] for i1 in 0:p1, i3 in 0:p3)
    a′ = sum(a[j1+i1,:,j3+i3].*w[j1+i1,:,j3+i3].*B1[1+i1].*B3[1+i3] for i1 in 0:p1, i3 in 0:p3) ./ w′
    return RationalBSplineManifold(a′,w′,(P2,))
end
@inline function (M::AbstractRationalBSplineManifold{3,p})(::Colon,t2::Real,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    j2 = intervalindex(P2,t2)
    j3 = intervalindex(P3,t3)
    B2 = bsplinebasisall(P2,j2,t2)
    B3 = bsplinebasisall(P3,j3,t3)
    w′ = sum(w[:,j2+i2,j3+i3]*B2[1+i2]*B3[1+i3] for i2 in 0:p2, i3 in 0:p3)
    a′ = sum(a[:,j2+i2,j3+i3].*w[:,j2+i2,j3+i3].*B2[1+i2].*B3[1+i3] for i2 in 0:p2, i3 in 0:p3) ./ w′
    return RationalBSplineManifold(a′,w′,(P1,))
end
