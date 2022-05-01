# B-spline manifold

abstract type AbstractManifold{Dim} end
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

@inline function (M::AbstractManifold{1})(t1)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    unsafe_mapping(M,t1)
end

@inline function (M::AbstractManifold{2})(t1,t2)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    unsafe_mapping(M,t1,t2)
end

@inline function (M::AbstractManifold{3})(t1,t2,t3)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    t3 in domain(Ps[3]) || throw(DomainError(t3, "The input $(t3) is out of range."))
    unsafe_mapping(M,t1,t2,t3)
end

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro
