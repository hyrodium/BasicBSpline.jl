# B-spline manifold

# TODO add more type parameter Deg, T
abstract type AbstractBSplineManifold{Dim} end

dim(M::AbstractBSplineManifold{Dim}) where Dim = Dim

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

bsplinespaces(M::BSplineManifold) = M.bsplinespaces
controlpoints(M::BSplineManifold) = M.controlpoints

@generated function (M::BSplineManifold{1,Deg,S,T})(t1::Real) where {Deg,S,T}
    p1, = Deg
    exs = Expr[]
    for j1 in 1:p1
        push!(exs, :(v .+= b1[$(1+j1)]*view(a,i1+$(j1),:)))
    end
    Expr(:block,
        :((P1,) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(i1 = intervalindex(P1,t1)),
        :(b1 = bsplinebasisall(P1, i1, t1)),
        :(v = b1[1]*view(a,i1,:)),
        exs...,
        :(return v)
    )
end

# TODO use @generated macro
function (M::BSplineManifold{2,(2,2),S,T})(t1,t2) where {S,T}
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    i1, i2 = intervalindex(P1,t1), intervalindex(P2,t2)
    b1, b2 = bsplinebasisall(P1, i1, t1), bsplinebasisall(P2, i2, t2)
    v = b1[1]*b2[1]*view(a,i1,i2,:)
    v .+= b1[2]*b2[1]*view(a,i1+1,i2,:)
    v .+= b1[3]*b2[1]*view(a,i1+2,i2,:)
    v .+= b1[1]*b2[2]*view(a,i1,i2+1,:)
    v .+= b1[2]*b2[2]*view(a,i1+1,i2+1,:)
    v .+= b1[3]*b2[2]*view(a,i1+2,i2+1,:)
    v .+= b1[1]*b2[3]*view(a,i1,i2+2,:)
    v .+= b1[2]*b2[3]*view(a,i1+1,i2+2,:)
    v .+= b1[3]*b2[3]*view(a,i1+2,i2+2,:)
    v
end
