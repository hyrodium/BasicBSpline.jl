# Rational B-spline manifold (NURBS)
abstract type AbstractRationalBSplineManifold{Dim, Deg} <: AbstractManifold{Dim} end

"""
Construct Rational B-spline manifold from given control points, weights and B-spline spaces.

# Examples
```jldoctest
julia> using StaticArrays, LinearAlgebra

julia> P = BSplineSpace{2}(KnotVector(0,0,0,1,1,1))
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

julia> M = RationalBSplineManifold(a,w,(P,));  # 1/4 arc

julia> M(0.3)
2-element SVector{2, Float64} with indices SOneTo(2):
 0.8973756499953727
 0.4412674277525845

julia> norm(M(0.3))
1.0
```
"""
struct RationalBSplineManifold{Dim,Deg,C,S<:Tuple,T} <: AbstractRationalBSplineManifold{Dim,Deg}
    bsplinespaces::S
    controlpoints::Array{C,Dim}
    weights::Array{T,Dim}
    function RationalBSplineManifold(a::Array{C,Dim},w::Array{T,Dim},Ps::S) where {S<:Tuple,C,Dim,T<:Real}
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
        new{Dim,Deg,C,S,T}(Ps,a,w)
    end
end

controlpoints(M::RationalBSplineManifold) = M.controlpoints
weights(M::RationalBSplineManifold) = M.weights
bsplinespaces(M::RationalBSplineManifold) = M.bsplinespaces

@generated function unsafe_mapping(M::RationalBSplineManifold{1,Deg},t1::Real) where {Deg}
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

@generated function unsafe_mapping(M::RationalBSplineManifold{2,Deg},t1::Real,t2::Real) where {Deg}
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

@generated function unsafe_mapping(M::RationalBSplineManifold{3,Deg},t1::Real,t2::Real,t3::Real) where {Deg}
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
