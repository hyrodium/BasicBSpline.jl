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
struct BSplineManifold{Dim,Deg,C,T,S<:NTuple{Dim, BSplineSpace{p,T} where p}} <: AbstractManifold{Dim}
    controlpoints::Array{C,Dim}
    bsplinespaces::S
    function BSplineManifold{Dim,Deg,C,T,S}(a::Array{C,Dim},P::S) where {S<:NTuple{Dim, BSplineSpace{p,T} where p},C} where {Dim, Deg, T}
        new{Dim,Deg,C,T,S}(a,P)
    end
end

function BSplineManifold(a::Array{C,Dim},P::S) where {S<:Tuple{BSplineSpace{p,T} where p, Vararg{BSplineSpace{p,T} where p}}} where {Dim, T, C}
    if size(a) != dim.(P)
        msg = "The size of control points array $(size(a)) and dimensions of B-spline spaces $(dim.(P)) must be equal."
        throw(DimensionMismatch(msg))
    end
    Deg = degree.(P)
    return BSplineManifold{Dim,Deg,C,T,S}(a, P)
end

function BSplineManifold(a::Array{C,Dim},P::S) where {S<:NTuple{Dim, BSplineSpace{p,T} where {p,T}},C} where {Dim}
    P′ = _promote_knottype(P)
    return BSplineManifold(a, P′)
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
    # Use `UnitRange` to support Julia v1.6 (LTS)
    # This can be replaced with `range` if we drop support for v1.6
    iter = CartesianIndices(UnitRange.(1, Deg .+ 1))
    exs = Expr[]
    for ci in iter
        ex = Expr(:call, [:getindex, :a, [:($(Symbol(:i,d))+$(ci[d]-1)) for d in 1:Dim]...]...)
        ex = Expr(:call, [:*, [:($(Symbol(:b,d))[$(ci[d])]) for d in 1:Dim]..., ex]...)
        ex = Expr(:+=, :v, ex)
        push!(exs, ex)
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

@generated function (M::BSplineManifold{Dim})(t::Vararg{Real, Dim}) where Dim
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

# 2dim
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

# TODO: The performance of this method can be improved.
function (M::BSplineManifold{Dim,Deg})(t::Union{Real, Colon}...) where {Dim, Deg}
    P = bsplinespaces(M)
    a = controlpoints(M)
    t_real = _remove_colon(t...)
    P_real = _remove_colon(_get_on_real.(P, t)...)
    P_colon = _remove_colon(_get_on_colon.(P, t)...)
    j_real = intervalindex.(P_real, t_real)
    B = bsplinebasisall.(P_real, j_real, t_real)
    Deg_real = _remove_colon(_get_on_real.(Deg, t)...)
    ci = CartesianIndices(UnitRange.(0, Deg_real))
    next = _replace_noncolon((), j_real, t...)
    a′ = view(a, next...) .* *(getindex.(B, 1)...)
    l = length(ci)
    for i in view(ci,2:l)
        next = _replace_noncolon((), j_real .+ i.I, t...)
        b = *(getindex.(B, i.I .+ 1)...)
        a′ .+= view(a, next...) .* b
    end
    return BSplineManifold(a′, P_colon)
end

@inline function (M::BSplineManifold{Dim})(::Vararg{Colon, Dim}) where Dim
    a = copy(controlpoints(M))
    Ps = bsplinespaces(M)
    return BSplineManifold(a,Ps)
end

@inline function (M::BSplineManifold{0})()
    a = controlpoints(M)
    return a[]
end
