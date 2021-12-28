# B-spline manifold with given type as a type of control points

struct CustomBSplineManifold{Dim,Deg,C,S<:Tuple} <: AbstractBSplineManifold{Dim,Deg}
    bsplinespaces::S
    controlpoints::Array{C,Dim}
    function CustomBSplineManifold(a::Array{C,Dim},Ps::S) where {S<:Tuple,C,Dim}
        if !all(isa.(Ps,AbstractBSplineSpace))
            # TODO: update error message
            error("invalid")
        end
        if size(a) != dim.(Ps)
            # TODO: update error message
            error("invalid")
        end
        p = degree.(Ps)
        new{Dim,p,C,S}(Ps,a)
    end
end

bsplinespaces(M::CustomBSplineManifold) = M.bsplinespaces
controlpoints(M::CustomBSplineManifold) = M.controlpoints

@generated function (M::CustomBSplineManifold{1,Deg,S,T})(t1::Real) where {Deg,S,T}
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

@generated function (M::CustomBSplineManifold{2,Deg,S,T})(t1,t2) where {Deg,S,T}
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

@generated function (M::CustomBSplineManifold{3,Deg,S,T})(t1,t2,t3) where {Deg,S,T}
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

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro
