# Rational B-spline manifold (NURBS)
"""
Construct Rational B-spline manifold from given control points, weights and B-spline spaces.

# Examples
```jldoctest
julia> using StaticArrays, LinearAlgebra

julia> P = BSplineSpace{2}(KnotVector([0,0,0,1,1,1]))
BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([0, 0, 0, 1, 1, 1]))

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
struct RationalBSplineManifold{Dim,Deg,C,W,T,S<:NTuple{Dim, BSplineSpace}} <: AbstractManifold{Dim}
    controlpoints::Array{C,Dim}
    weights::Array{W,Dim}
    bsplinespaces::S
    function RationalBSplineManifold{Dim,Deg,C,W,T,S}(a::Array{C,Dim},w::Array{W,Dim},P::S) where {S<:NTuple{Dim, BSplineSpace{p,T} where p},C,W} where {Dim, Deg, T}
        new{Dim,Deg,C,W,T,S}(a,w,P)
    end
end

function RationalBSplineManifold(a::Array{C,Dim},w::Array{W,Dim},P::S) where {S<:Tuple{BSplineSpace{p,T} where p, Vararg{BSplineSpace{p,T} where p}}} where {Dim, T, C, W}
    if size(a) != dim.(P)
        msg = "The size of control points array $(size(a)) and dimensions of B-spline spaces $(dim.(P)) must be equal."
        throw(DimensionMismatch(msg))
    end
    if size(w) != dim.(P)
        msg = "The size of weights array $(size(w)) and dimensions of B-spline spaces $(dim.(P)) must be equal."
        throw(DimensionMismatch(msg))
    end
    Deg = degree.(P)
    return RationalBSplineManifold{Dim,Deg,C,W,T,S}(a, w, P)
end

function RationalBSplineManifold(a::Array{C,Dim},w::Array{W,Dim},P::S) where {S<:NTuple{Dim, BSplineSpace{p,T} where {p,T}},C,W} where {Dim}
    P′ = _promote_knottype(P)
    return RationalBSplineManifold(a, w, P′)
end

RationalBSplineManifold(a::Array{C,Dim},w::Array{T,Dim},Ps::Vararg{BSplineSpace, Dim}) where {C,Dim,T<:Real} = RationalBSplineManifold(a,w,Ps)

Base.:(==)(M1::RationalBSplineManifold, M2::RationalBSplineManifold) = (bsplinespaces(M1)==bsplinespaces(M2)) & (controlpoints(M1)==controlpoints(M2)) & (weights(M1)==weights(M2))

function Base.hash(M::RationalBSplineManifold{0}, h::UInt)
    hash(RationalBSplineManifold{0}, hash(weights(M), hash(controlpoints(M), h)))
end

function Base.hash(M::RationalBSplineManifold, h::UInt)
    hash(xor(hash.(bsplinespaces(M), h)...), hash(weights(M), hash(controlpoints(M), h)))
end

controlpoints(M::RationalBSplineManifold) = M.controlpoints
weights(M::RationalBSplineManifold) = M.weights
bsplinespaces(M::RationalBSplineManifold) = M.bsplinespaces

@generated function unbounded_mapping(M::RationalBSplineManifold{Dim,Deg}, t::Vararg{Real,Dim}) where {Dim,Deg}
    # Use `UnitRange` to support Julia v1.6 (LTS)
    # This can be replaced with `range` if we drop support for v1.6
    iter = CartesianIndices(UnitRange.(1, Deg .+ 1))
    exs = Expr[]
    for ci in iter
        ex_w = Expr(:call, [:getindex, :w, [:($(Symbol(:i,d))+$(ci[d]-1)) for d in 1:Dim]...]...)
        ex_a = Expr(:call, [:getindex, :a, [:($(Symbol(:i,d))+$(ci[d]-1)) for d in 1:Dim]...]...)
        ex = Expr(:call, [:*, [:($(Symbol(:b,d))[$(ci[d])]) for d in 1:Dim]..., ex_w]...)
        ex = Expr(:+=, :u, ex)
        push!(exs, ex)
        ex = Expr(:call, [:getindex, :a, [:($(Symbol(:i,d))+$(ci[d]-1)) for d in 1:Dim]...]...)
        ex = Expr(:call, [:*, [:($(Symbol(:b,d))[$(ci[d])]) for d in 1:Dim]..., ex_w, ex_a]...)
        ex = Expr(:+=, :v, ex)
        push!(exs, ex)
    end
    exs[1].head = :(=)
    exs[2].head = :(=)
    Expr(
        :block,
        Expr(:(=), Expr(:tuple, [Symbol(:P, i) for i in 1:Dim]...), :(bsplinespaces(M))),
        Expr(:(=), Expr(:tuple, [Symbol(:t, i) for i in 1:Dim]...), :t),
        :(a = controlpoints(M)),
        :(w = weights(M)),
        Expr(:(=), Expr(:tuple, [Symbol(:i, i) for i in 1:Dim]...), Expr(:tuple, [:(intervalindex($(Symbol(:P, i)), $(Symbol(:t, i)))) for i in 1:Dim]...)),
        Expr(:(=), Expr(:tuple, [Symbol(:b, i) for i in 1:Dim]...), Expr(:tuple, [:(bsplinebasisall($(Symbol(:P, i)), $(Symbol(:i, i)), $(Symbol(:t, i)))) for i in 1:Dim]...)),
        exs...,
        :(return v/u)
    )
end

## currying
# 1dim
@inline function (M::RationalBSplineManifold{1})(::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end

# 2dim
@inline function (M::RationalBSplineManifold{2})(::Colon,::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end
@inline function (M::RationalBSplineManifold{2,p})(t1::Real,::Colon) where p
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
@inline function (M::RationalBSplineManifold{2,p})(::Colon,t2::Real) where p
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
@inline function (M::RationalBSplineManifold{3})(::Colon,::Colon,::Colon)
    a = copy(controlpoints(M))
    w = copy(weights(M))
    Ps = bsplinespaces(M)
    return RationalBSplineManifold(a,w,Ps)
end
@inline function (M::RationalBSplineManifold{3,p})(t1::Real,::Colon,::Colon) where p
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
@inline function (M::RationalBSplineManifold{3,p})(::Colon,t2::Real,::Colon) where p
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
@inline function (M::RationalBSplineManifold{3,p})(::Colon,::Colon,t3::Real) where p
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
@inline function (M::RationalBSplineManifold{3,p})(t1::Real,t2::Real,::Colon) where p
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
@inline function (M::RationalBSplineManifold{3,p})(t1::Real,::Colon,t3::Real) where p
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
@inline function (M::RationalBSplineManifold{3,p})(::Colon,t2::Real,t3::Real) where p
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
