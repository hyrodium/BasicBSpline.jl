# B-spline manifold

# TODO add more type parameter Deg, T
abstract type AbstractBSplineManifold{Dim} end

# struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim,Deg,T}
struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim}
    bsplinespaces::S
    controlpoints::Array{T,Dim₊₁}
    function BSplineManifold(Ps::S,a::Array{T,Dim₊₁}) where {S<:Tuple,Dim₊₁,T<:Real}
        if !all(isa.(Ps,AbstractBSplineSpace))
            # TODO: update error message
            error("invalid")
        end
        if size(a)[1:Dim₊₁-1] != dim.(Ps)
            # TODO: update error message
            error("invalid")
        end
        d = length(Ps)
        p = degree.(Ps)
        new{d,p,T,S,d+1}(Ps,a)
    end
end
