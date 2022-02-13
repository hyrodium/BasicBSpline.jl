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
