# B-spline manifold

abstract type AbstractManifold{Dim} end

# Broadcast like a scalar
Base.Broadcast.broadcastable(M::AbstractManifold) = Ref(M)

dim(::AbstractManifold{Dim}) where Dim = Dim

@doc raw"""
Construct B-spline manifold from given control points and B-spline spaces.

# Examples
```jldoctest
julia> using StaticArrays

julia> P = BSplineSpace{2}(KnotVector([0,0,0,1,1,1]))
BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([0, 0, 0, 1, 1, 1]))

julia> a = [SVector(1,0), SVector(1,1), SVector(0,1)]
3-element Vector{SVector{2, Int64}}:
 [1, 0]
 [1, 1]
 [0, 1]

julia> M = BSplineManifold(a, P);


julia> M(0.4)
2-element SVector{2, Float64} with indices SOneTo(2):
 0.84
 0.64

julia> M(1.2)
ERROR: DomainError with 1.2:
The input 1.2 is out of range.
[...]
```
"""
struct BSplineManifold{Dim,Deg,C,S<:NTuple{Dim, BSplineSpace}} <: AbstractManifold{Dim}
    bsplinespaces::S
    controlpoints::Array{C,Dim}
    function BSplineManifold(a::Array{C,Dim},Ps::S) where {S<:NTuple{Dim, BSplineSpace},C} where Dim
        if size(a) != dim.(Ps)
            msg = "The size of control points array $(size(a)) and dimensions of B-spline spaces $(dim.(Ps)) must be equal."
            throw(DimensionMismatch(msg))
        end
        Deg = degree.(Ps)
        new{Dim,Deg,C,S}(Ps,a)
    end
end

BSplineManifold(a::Array{C,Dim},Ps::Vararg{BSplineSpace, Dim}) where {C,Dim} = BSplineManifold(a,Ps)

Base.:(==)(M1::AbstractManifold, M2::AbstractManifold) = (bsplinespaces(M1)==bsplinespaces(M2)) & (controlpoints(M1)==controlpoints(M2))

bsplinespaces(M::BSplineManifold) = M.bsplinespaces
controlpoints(M::BSplineManifold) = M.controlpoints

function Base.hash(M::BSplineManifold{0}, h::UInt)
    hash(BSplineManifold{0}, hash(controlpoints(M), h))
end

function Base.hash(M::BSplineManifold, h::UInt)
    hash(xor(hash.(bsplinespaces(M), h)...), hash(controlpoints(M), h))
end

@doc raw"""

    unbounded_mapping(M::BSplineManifold{Dim}, t::Vararg{Real,Dim})

# Examples
```jldoctest
julia> P = BSplineSpace{1}(KnotVector([0,0,1,1]))
BSplineSpace{1, Int64, KnotVector{Int64}}(KnotVector([0, 0, 1, 1]))

julia> domain(P)
0 .. 1

julia> M = BSplineManifold([0,1], P);


julia> unbounded_mapping(M, 0.1)
0.1

julia> M(0.1)
0.1

julia> unbounded_mapping(M, 1.2)
1.2

julia> M(1.2)
ERROR: DomainError with 1.2:
The input 1.2 is out of range.
[...]
```
"""
function unbounded_mapping end

@generated function unbounded_mapping(M::BSplineManifold{Dim,Deg}, t::Vararg{Real,Dim}) where {Dim,Deg}
    iter = CartesianIndices(range.(1, Deg .+ 1))
    exs = Expr[]
    for i in iter
        a = [:getindex, :a, [:($(Symbol(:i,d))+$(i[d]-1)) for d in 1:Dim]...]
        b = Expr(:call, a...)
        c = [:*, [:($(Symbol(:b,d))[$(i[d])]) for d in 1:Dim]..., b]
        d = Expr(:call, c...)
        e = Expr(:+=, :v, d)
        push!(exs, e)
    end
    exs[1].head = :(=)
    Expr(
        :block,
        Expr(:(=), Expr(:tuple, [Symbol(:P, i) for i in 1:Dim]...), :(bsplinespaces(M))),
        Expr(:(=), Expr(:tuple, [Symbol(:t, i) for i in 1:Dim]...), :t),
        :(a = controlpoints(M)),
        Expr(:(=), Expr(:tuple, [Symbol(:i, i) for i in 1:Dim]...), Expr(:tuple, [:(intervalindex($(Symbol(:P, i)), $(Symbol(:t, i)))) for i in 1:Dim]...)),
        Expr(:(=), Expr(:tuple, [Symbol(:b, i) for i in 1:Dim]...), Expr(:tuple, [:(bsplinebasisall($(Symbol(:P, i)), $(Symbol(:i, i)), $(Symbol(:t, i)))) for i in 1:Dim]...)),
        exs...,
        :(return v)
    )
end

@generated function (M::AbstractManifold{Dim})(t::Vararg{Real, Dim}) where Dim
    Ps = [Symbol(:P,i) for i in 1:Dim]
    P = Expr(:tuple, Ps...)
    ts = [Symbol(:t,i) for i in 1:Dim]
    T = Expr(:tuple, ts...)
    exs = [:($(Symbol(:t,i)) in domain($(Symbol(:P,i))) || throw(DomainError($(Symbol(:t,i)), "The input "*string($(Symbol(:t,i)))*" is out of range."))) for i in 1:Dim]
    ret = Expr(:call,:unbounded_mapping,:M,[Symbol(:t,i) for i in 1:Dim]...)
    Expr(
        :block,
        :($(Expr(:meta, :inline))),
        :($T = t),
        :($P = bsplinespaces(M)),
        exs...,
        :(return $(ret))
    )
end


## currying
# 1dim
@inline function (M::BSplineManifold{1})(::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end

# 2dim
@inline function (M::BSplineManifold{2})(::Colon,::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end
@inline function (M::BSplineManifold{2,p})(t1::Real,::Colon) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    b = sum(a[j1+i1,:]*B1[1+i1] for i1 in 0:p1)
    return BSplineManifold(b,(P2,))
end
@inline function (M::BSplineManifold{2,p})(::Colon,t2::Real) where p
    p1, p2 = p
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[:,j2+i2]*B2[1+i2] for i2 in 0:p2)
    return BSplineManifold(b,(P1,))
end

# 3dim
@inline function (M::BSplineManifold{3})(::Colon,::Colon,::Colon)
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end
@inline function (M::BSplineManifold{3,p})(t1::Real,::Colon,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    B1 = bsplinebasisall(P1,j1,t1)
    b = sum(a[j1+i1,:,:]*B1[1+i1] for i1 in 0:p1)
    return BSplineManifold(b,(P2,P3))
end
@inline function (M::BSplineManifold{3,p})(::Colon,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[:,j2+i2,:]*B2[1+i2] for i2 in 0:p2)
    return BSplineManifold(b,(P1,P3))
end
@inline function (M::BSplineManifold{3,p})(::Colon,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j3 = intervalindex(P3,t3)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[:,:,j3+i3]*B3[1+i3] for i3 in 0:p3)
    return BSplineManifold(b,(P1,P2))
end
@inline function (M::BSplineManifold{3,p})(t1::Real,t2::Real,::Colon) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    j2 = intervalindex(P2,t2)
    B1 = bsplinebasisall(P1,j1,t1)
    B2 = bsplinebasisall(P2,j2,t2)
    b = sum(a[j1+i1,j2+i2,:]*B1[1+i1]*B2[1+i2] for i1 in 0:p1, i2 in 0:p2)
    return BSplineManifold(b,(P3,))
end
@inline function (M::BSplineManifold{3,p})(t1::Real,::Colon,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j1 = intervalindex(P1,t1)
    j3 = intervalindex(P3,t3)
    B1 = bsplinebasisall(P1,j1,t1)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[j1+i1,:,j3+i3]*B1[1+i1]*B3[1+i3] for i1 in 0:p1, i3 in 0:p3)
    return BSplineManifold(b,(P2,))
end
@inline function (M::BSplineManifold{3,p})(::Colon,t2::Real,t3::Real) where p
    p1, p2, p3 = p
    P1, P2, P3 = bsplinespaces(M)
    a = controlpoints(M)
    j2 = intervalindex(P2,t2)
    j3 = intervalindex(P3,t3)
    B2 = bsplinebasisall(P2,j2,t2)
    B3 = bsplinebasisall(P3,j3,t3)
    b = sum(a[:,j2+i2,j3+i3]*B2[1+i2]*B3[1+i3] for i2 in 0:p2, i3 in 0:p3)
    return BSplineManifold(b,(P1,))
end

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro
