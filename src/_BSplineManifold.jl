# B-spline manifold

abstract type AbstractBSplineManifold{Dim, Deg} end

dim(::AbstractBSplineManifold{Dim}) where Dim = Dim

# struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim,Deg,T}
struct BSplineManifold{Dim,Deg,T,S<:Tuple,Dim₊₁} <: AbstractBSplineManifold{Dim,Deg}
    bsplinespaces::S
    controlpoints::Array{T,Dim₊₁}
    function BSplineManifold(a::Array{T,Dim₊₁},Ps::S) where {S<:Tuple,Dim₊₁,T<:Real}
        for P in Ps
            if !(P isa AbstractBSplineSpace)
                throw(TypeError(:CustomBSplineManifold,AbstractBSplineSpace,P))
            end
        end
        if size(a)[1:Dim₊₁-1] != dim.(Ps)
            msg = "The size of control points array $(size(a)[1:Dim₊₁-1]) and dimensions of B-spline spaces $(dim.(Ps)) must be equal."
            throw(DimensionMismatch(msg))
        end
        Dim = length(Ps)
        Deg = degree.(Ps)
        new{Dim,Deg,T,S,Dim+1}(Ps,a)
    end
end

bsplinespaces(M::BSplineManifold) = M.bsplinespaces
controlpoints(M::BSplineManifold) = M.controlpoints

@generated function _mapping(M::BSplineManifold{1,Deg,S,T},t1) where {Deg,S,T}
    p1, = Deg
    exs = Expr[]
    for j1 in 1:p1
        push!(exs, :(v .+= b1[$(1+j1)]*view(a,i1+$(j1),:)))
    end
    Expr(
        :block,
        :((P1,) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :(i1 = intervalindex(P1,t1)),
        :(b1 = bsplinebasisall(P1,i1,t1)),
        :(v = b1[1]*view(a,i1,:)),
        exs...,
        :(return v)
    )
end

@generated function _mapping(M::BSplineManifold{2,Deg,S,T},t1,t2) where {Deg,S,T}
    p1, p2 = Deg
    exs = Expr[]
    for j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(v .+= b1[$(j1)]*b2[$(j2)]*view(a,i1+$(j1-1),i2+$(j2-1),:)))
    end
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :((i1, i2) = (intervalindex(P1,t1), intervalindex(P2,t2))),
        :((b1, b2) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2))),
        :(v = b1[1]*b2[1]*view(a,i1,i2,:)),
        exs...,
        :(return v)
    )
end

@generated function _mapping(M::BSplineManifold{3,Deg,S,T},t1,t2,t3) where {Deg,S,T}
    p1, p2, p3 = Deg
    exs = Expr[]
    for j3 in 1:p3+1, j2 in 1:p2+1, j1 in 1:p1+1
        push!(exs, :(v .+= b1[$(j1)]*b2[$(j2)]*b3[$(j3)]*view(a,i1+$(j1-1),i2+$(j2-1),i3+$(j3-1),:)))
    end
    deleteat!(exs,1)
    Expr(
        :block,
        :((P1, P2, P3) = bsplinespaces(M)),
        :(a = controlpoints(M)),
        :((i1, i2, i3) = (intervalindex(P1,t1), intervalindex(P2,t2), intervalindex(P3,t3))),
        :((b1, b2, b3) = (bsplinebasisall(P1,i1,t1), bsplinebasisall(P2,i2,t2), bsplinebasisall(P3,i3,t3))),
        :(v = b1[1]*b2[1]*b3[1]*view(a,i1,i2,i3,:)),
        exs...,
        :(return v)
    )
end

function (M::AbstractBSplineManifold{Dim})(ts::Real...) where Dim
    Ps = bsplinespaces(M)
    for i in 1:Dim
        ts[i] in domain(Ps[i]) || throw(DomainError(ts[i], "The input is out of range."))
    end
    _mapping(M,ts...)
end

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro
