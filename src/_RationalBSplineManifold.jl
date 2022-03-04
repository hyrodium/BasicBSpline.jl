# Rational B-spline manifold (NURBS)
abstract type AbstractRationalBSplineManifold{Dim, Deg} <: AbstractManifold{Dim} end

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
        push!(exs, :(v += b1[$(1+j1)]*getindex(a,i1+$(j1))*getindex(w,i1+$(j1))))
    end
    Expr(:block,
        :((P1,) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(w = weights(M)),
        :(i1 = intervalindex(P1,t1)),
        :(b1 = bsplinebasisall(P1,i1,t1)),
        :(u = b1[1]*getindex(w,i1)),
        :(v = b1[1]*getindex(a,i1)*getindex(w,i1)),
        exs...,
        :(return v/u)
    )
end
