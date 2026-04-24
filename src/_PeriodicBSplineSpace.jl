# Periodic B-spline space

@doc raw"""
Construct periodic B-spline space from given polynominal degree and knot vector.

The knot vector ``k = (k_1, \dots, k_{n+1})`` represents one period of the periodic knot
sequence. The period length is ``L = k_{n+1} - k_1``, and the space has ``n`` cyclically
indexed basis functions (so ``\dim = \#k - 1``, independent of the degree ``p``).

# Examples
```jldoctest
julia> p = 2
2

julia> k = KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])
KnotVector([0.0, 1.0, 2.0, 3.0, 4.0])

julia> PeriodicBSplineSpace{p}(k)
PeriodicBSplineSpace{2, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.0, 2.0, 3.0, 4.0]))
```
"""
struct PeriodicBSplineSpace{p, T<:Real, K<:AbstractKnotVector{T}} <: AbstractBSplineSpace{p, T}
    knotvector::K
    global unsafe_periodicbsplinespace(::Val{p}, k::K) where {p, K<:AbstractKnotVector{T}} where {T<:Real} = new{p,T,K}(k)
end

function PeriodicBSplineSpace{p}(k::AbstractKnotVector) where p
    if p < 0
        throw(DomainError(p, "degree of polynominal must be non-negative"))
    end
    if length(k) < p + 2
        throw(DomainError(k, "knot vector has too few knots for the given degree"))
    end
    return unsafe_periodicbsplinespace(Val{p}(), k)
end
function PeriodicBSplineSpace{p,T}(k::AbstractKnotVector{T}) where {p, T}
    return PeriodicBSplineSpace{p}(k)
end
function PeriodicBSplineSpace{p,T1}(k::AbstractKnotVector{T2}) where {p, T1, T2}
    return PeriodicBSplineSpace{p,T1}(AbstractKnotVector{T1}(k))
end

"""
Convert PeriodicBSplineSpace to PeriodicBSplineSpace
"""
function PeriodicBSplineSpace(P::PeriodicBSplineSpace{p}) where p
    return P
end
function PeriodicBSplineSpace{p}(P::PeriodicBSplineSpace{p}) where p
    return P
end
function PeriodicBSplineSpace{p,T}(P::PeriodicBSplineSpace{p}) where {p,T}
    return PeriodicBSplineSpace{p}(AbstractKnotVector{T}(knotvector(P)))
end
function PeriodicBSplineSpace{p,T,K}(P::PeriodicBSplineSpace{p}) where {p,T,K}
    return PeriodicBSplineSpace{p,T}(K(knotvector(P)))
end
function PeriodicBSplineSpace{p,T,K}(k::AbstractKnotVector) where {p,T,K}
    return PeriodicBSplineSpace{p,T}(K(k))
end

# Broadcast like a scalar
Base.Broadcast.broadcastable(P::PeriodicBSplineSpace) = Ref(P)

# Equality
@inline Base.:(==)(P1::PeriodicBSplineSpace{p}, P2::PeriodicBSplineSpace{p}) where p = knotvector(P1) == knotvector(P2)
@inline Base.:(==)(P1::PeriodicBSplineSpace{p1}, P2::PeriodicBSplineSpace{p2}) where {p1, p2} = false

Base.copy(P::PeriodicBSplineSpace{p}) where p = PeriodicBSplineSpace{p}(copy(P.knotvector))

function Base.hash(P::PeriodicBSplineSpace{p}, h::UInt) where p
    k = knotvector(P)
    return hash(PeriodicBSplineSpace{p}, hash(_vec(k), h))
end

bsplinespace(P::PeriodicBSplineSpace) = P

@inline function knotvector(P::PeriodicBSplineSpace)
    return P.knotvector
end

@doc raw"""
Return the period of a periodic B-spline space, ``L = k_{n+1} - k_1``.
"""
function period(P::PeriodicBSplineSpace)
    k = knotvector(P)
    return k[end] - k[1]
end

@doc raw"""
Return the fundamental domain ``[k_1, k_{n+1}]`` of a periodic B-spline space.
"""
function domain(P::PeriodicBSplineSpace)
    k = knotvector(P)
    return k[1]..k[end]
end

@doc raw"""
Return dimension of a periodic B-spline space.
```math
\dim(\mathcal{P}^{\mathrm{per}}[p,k]) = \# k - 1
```
(Independent of ``p`` — the last knot is interpreted as the period marker.)

# Examples
```jldoctest
julia> dim(PeriodicBSplineSpace{1}(KnotVector([1,2,3,4,5,6,7])))
6

julia> dim(PeriodicBSplineSpace{2}(KnotVector([1,2,3,4,5,6,7])))
6
```
"""
function dim(P::PeriodicBSplineSpace)
    k = knotvector(P)
    return length(k) - 1
end

##########################
#= Cyclic knot indexing =#
##########################

@doc raw"""
Reduce `t` to the fundamental period `[k_1, k_{n+1})`.
"""
@inline function _periodic_reduce(P::PeriodicBSplineSpace, t::Real)
    k = knotvector(P)
    L = period(P)
    t0 = k[1]
    return mod(t - t0, L) + t0
end

@doc raw"""
Cyclically-extended knot value: for any integer `j`, return `k[mod1(j, n)] + L * fld(j-1, n)`,
where `n = length(k) - 1` is the number of distinct knots in one period.
"""
@inline function _periodic_knot(P::PeriodicBSplineSpace, j::Integer)
    k = knotvector(P)
    L = period(P)
    n = length(k) - 1
    q = fld(j - 1, n)
    r = mod(j - 1, n) + 1
    return k[r] + L * q
end
