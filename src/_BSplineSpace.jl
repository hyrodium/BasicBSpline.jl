# B-spline space

abstract type AbstractFunctionSpace{T} end

@doc raw"""
Construct B-spline space from given polynominal degree and knot vector.
```math
\mathcal{P}[p,k]
```

# Examples
```jldoctest
julia> p = 2
2

julia> k = KnotVector([1,3,5,6,8,9])
KnotVector([1, 3, 5, 6, 8, 9])

julia> BSplineSpace{p}(k)
BSplineSpace{2, Int64, KnotVector{Int64}}(KnotVector([1, 3, 5, 6, 8, 9]))
```
"""
struct BSplineSpace{p, T<:Real, K<:AbstractKnotVector{T}} <: AbstractFunctionSpace{T}
    knotvector::K
    global unsafe_bsplinespace(::Val{p}, k::K) where {p, K<:AbstractKnotVector{T}} where {T<:Real} = new{p,T,K}(k)
end
function BSplineSpace{p}(k::AbstractKnotVector) where p
    if p < 0
        throw(DomainError(p, "degree of polynominal must be non-negative"))
    end
    return unsafe_bsplinespace(Val{p}(), k)
end
function BSplineSpace{p,T}(k::AbstractKnotVector{T}) where {p, T}
    return BSplineSpace{p}(k)
end
function BSplineSpace{p,T1}(k::AbstractKnotVector{T2}) where {p, T1, T2}
    return BSplineSpace{p,T1}(AbstractKnotVector{T1}(k))
end

const UniformBSplineSpace{p, T, R} = BSplineSpace{p, T, UniformKnotVector{T, R}}

"""
Convert BSplineSpace to BSplineSpace
"""
function BSplineSpace(P::BSplineSpace{p}) where p
    return P
end
function BSplineSpace{p}(P::BSplineSpace{p}) where p
    return P
end
function BSplineSpace{p,T}(P::BSplineSpace{p}) where {p,T}
    return BSplineSpace{p}(AbstractKnotVector{T}(knotvector(P)))
end
function BSplineSpace{p,T,K}(P::BSplineSpace{p}) where {p,T,K}
    return BSplineSpace{p,T}(K(knotvector(P)))
end
function BSplineSpace{p,T,K}(k::AbstractKnotVector) where {p,T,K}
    return BSplineSpace{p,T}(K(k))
end

# Broadcast like a scalar
Base.Broadcast.broadcastable(P::BSplineSpace) = Ref(P)
Base.iterate(P::AbstractFunctionSpace) = (P, nothing)
Base.iterate(::AbstractFunctionSpace, ::Any) = nothing

# Equality
@inline Base.:(==)(P1::BSplineSpace{p}, P2::BSplineSpace{p}) where p = knotvector(P1) == knotvector(P2)
@inline Base.:(==)(P1::BSplineSpace{p1}, P2::BSplineSpace{p2}) where {p1, p2} = false

Base.copy(P::BSplineSpace{p}) where p = BSplineSpace{p}(copy(P.knotvector))

bsplinespace(P::BSplineSpace) = P

@inline function degree(::BSplineSpace{p}) where p
    return p
end

@inline function knotvector(P::BSplineSpace)
    return P.knotvector
end

@generated function _promote_knottype(P::NTuple{Dim,BSplineSpace}) where Dim
    return Expr(
        :block,
        Expr(:(=), Expr(:tuple, [Symbol(:P, i) for i in 1:Dim]...), :P),
        Expr(:(=), Expr(:tuple, [Symbol(:k, i) for i in 1:Dim]...), Expr(:tuple, [:(knotvector($(Symbol(:P, i)))) for i in 1:Dim]...)),
        Expr(:(=), :T, Expr(:call, :promote_type, [:(eltype($(Symbol(:k, i)))) for i in 1:Dim]...)),
        Expr(:(=), Expr(:tuple, [Symbol(:k, i, :(var"′")) for i in 1:Dim]...), Expr(:tuple, [:(AbstractKnotVector{T}($(Symbol(:k, i)))) for i in 1:Dim]...)),
        Expr(:(=), :P′, Expr(:tuple, [:(BSplineSpace{degree($(Symbol(:P, i)))}($(Symbol(:k, i, :(var"′"))))) for i in 1:Dim]...)),
        :(return P′)
    )
end

function domain(P::BSplineSpace)
    p = degree(P)
    k = knotvector(P)
    return k[1+p]..k[end-p]
end

@doc raw"""
Return dimention of a B-spline space.
```math
\dim(\mathcal{P}[p,k])
=\# k - p -1
```

# Examples
```jldoctest
julia> dim(BSplineSpace{1}(KnotVector([1,2,3,4,5,6,7])))
5

julia> dim(BSplineSpace{1}(KnotVector([1,2,4,4,4,6,7])))
5

julia> dim(BSplineSpace{1}(KnotVector([1,2,3,5,5,5,7])))
5
```
"""
function dim(bsplinespace::BSplineSpace{p}) where p
    k = knotvector(bsplinespace)
    return length(k) - p - 1
end

@doc raw"""
Check inclusive relationship between B-spline spaces.
```math
\mathcal{P}[p,k]
\subseteq\mathcal{P}[p',k']
```

# Examples
```jldoctest
julia> P1 = BSplineSpace{1}(KnotVector([1,3,5,8]));


julia> P2 = BSplineSpace{1}(KnotVector([1,3,5,6,8,9]));


julia> P3 = BSplineSpace{2}(KnotVector([1,1,3,3,5,5,8,8]));


julia> P1 ⊆ P2
true

julia> P1 ⊆ P3
true

julia> P2 ⊆ P3
false

julia> P2 ⊈ P3
true
```
"""
function Base.issubset(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
    p₊ = p′ - p
    p₊ < 0 && return false
    k = knotvector(P)
    k′ = knotvector(P′)
    l = length(k)
    i = 1
    c = 0
    for j in 1:length(k′)
        if l < i
            return true
        end
        if k[i] < k′[j]
            return false
        end
        if k[i] == k′[j]
            c += 1
            if c == p₊+1
                i += 1
                c = 0
            elseif (i ≠ l) && (k[i] == k[i+1])
                i += 1
                c = 0
            end
        end
    end
    return l < i
    # The above implementation is the same as `k + p₊ * unique(k) ⊆ k′`, but that's much faster.
end

@doc raw"""
Check inclusive relationship between B-spline spaces.
```math
\mathcal{P}[p,k]
\sqsubseteq\mathcal{P}[p',k']
\Leftrightarrow
\mathcal{P}[p,k]|_{[k_{p+1},k_{l-p}]}
\subseteq\mathcal{P}[p',k']|_{[k'_{p'+1},k'_{l'-p'}]}
```
"""
function issqsubset(P::BSplineSpace{p}, P′::BSplineSpace{p′}) where {p, p′}
    p₊ = p′ - p
    p₊ < 0 && return false

    k = knotvector(P)
    k′ = knotvector(P′)
    !(k[1+p] == k′[1+p′] < k′[end-p′] == k[end-p]) && return false

    l = length(k)
    l′ = length(k′)
    inner_knotvector = view(k, p+2:l-p-1)
    inner_knotvector′ = view(k′, p′+2:l′-p′-1)

    _P = BSplineSpace{p}(inner_knotvector)
    _P′ = BSplineSpace{p′}(inner_knotvector′)
    return _P ⊆ _P′
end

const ⊑ = issqsubset
⊒(l, r) = r ⊑ l
⋢(l, r) = !⊑(l, r)
⋣(l, r) = r ⋢ l

≃(P1::BSplineSpace, P2::BSplineSpace) = (P1 ⊑ P2) & (P2 ⊑ P1)

Base.:⊊(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊆ B)
Base.:⊋(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊇ B)
⋤(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊑ B)
⋥(A::AbstractFunctionSpace, B::AbstractFunctionSpace) = (A ≠ B) & (A ⊒ B)

function isdegenerate_R(P::BSplineSpace{p}, i::Integer) where p
    k = knotvector(P)
    return k[i] == k[i+p+1]
end

function isdegenerate_I(P::BSplineSpace{p}, i::Integer) where p
    return iszero(width(bsplinesupport(P,i) ∩ domain(P)))
end

function _iszeros_R(P::BSplineSpace{p}) where p
    return [isdegenerate_R(P,i) for i in 1:dim(P)]
end

function _iszeros_I(P::BSplineSpace{p}) where p
    return [isdegenerate_I(P,i) for i in 1:dim(P)]
end

@doc raw"""
Check if given B-spline space is non-degenerate.

# Examples
```jldoctest
julia> isnondegenerate(BSplineSpace{2}(KnotVector([1,3,5,6,8,9])))
true

julia> isnondegenerate(BSplineSpace{1}(KnotVector([1,3,3,3,8,9])))
false
```
"""
isnondegenerate(P::BSplineSpace)

@doc raw"""
Check if given B-spline space is degenerate.

# Examples
```jldoctest
julia> isdegenerate(BSplineSpace{2}(KnotVector([1,3,5,6,8,9])))
false

julia> isdegenerate(BSplineSpace{1}(KnotVector([1,3,3,3,8,9])))
true
```
"""
isdegenerate(P::BSplineSpace)

for (f, fnon) in ((:isdegenerate_R, :isnondegenerate_R), (:isdegenerate_I, :isnondegenerate_I))
    @eval function $f(P::AbstractFunctionSpace)
        for i in 1:dim(P)
            $f(P,i) && return true
        end
        return false
    end
    @eval $fnon(P, i::Integer) = !$f(P, i)
    @eval $fnon(P) = !$f(P)
end
const isdegenerate = isdegenerate_R
const isnondegenerate = isnondegenerate_R

function _nondegeneratize_R(P::BSplineSpace{p}) where p
    k = copy(KnotVector(knotvector(P)))
    I = Int[]
    for i in 1:dim(P)
        isdegenerate_R(P,i) && push!(I, i)
    end
    deleteat!(BasicBSpline._vec(k), I)
    return BSplineSpace{p}(k)
end
function _nondegeneratize_I(P::BSplineSpace{p}) where p
    k = copy(KnotVector(knotvector(P)))
    l = length(k)
    I = Int[]
    for i in 1:dim(P)
        if isdegenerate_I(P,i)
            if k[l-p] < k[i+p+1]
                push!(I, i+p+1)
            else
                push!(I, i)
            end
        end
    end
    deleteat!(BasicBSpline._vec(k), I)
    return BSplineSpace{p}(k)
end

"""
Exact dimension of a B-spline space.

# Examples
```jldoctest
julia> exactdim_R(BSplineSpace{1}(KnotVector([1,2,3,4,5,6,7])))
5

julia> exactdim_R(BSplineSpace{1}(KnotVector([1,2,4,4,4,6,7])))
4

julia> exactdim_R(BSplineSpace{1}(KnotVector([1,2,3,5,5,5,7])))
4
```
"""
function exactdim_R(P::BSplineSpace)
    n = dim(P)
    for i in 1:dim(P)
        n -= isdegenerate_R(P,i)
    end
    return n
    # The above implementation is the same as `dim(P)-sum(_iszeros_R(P))`, but that's much faster.
end

"""
Exact dimension of a B-spline space.

# Examples
```jldoctest
julia> exactdim_I(BSplineSpace{1}(KnotVector([1,2,3,4,5,6,7])))
5

julia> exactdim_I(BSplineSpace{1}(KnotVector([1,2,4,4,4,6,7])))
4

julia> exactdim_I(BSplineSpace{1}(KnotVector([1,2,3,5,5,5,7])))
3
```
"""
function exactdim_I(P::BSplineSpace)
    n = dim(P)
    for i in 1:dim(P)
        n -= isdegenerate_I(P,i)
    end
    return n
    # The above implementation is the same as `dim(P)-sum(_iszeros_I(P))`, but that's much faster.
end

const exactdim = exactdim_R

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

julia> P = BSplineSpace{2}(k)
BSplineSpace{2, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> bsplinesupport_R(P,1)
0.0 .. 5.5

julia> bsplinesupport_R(P,2)
1.5 .. 8.0
```
"""
function bsplinesupport_R(P::BSplineSpace{p}, i::Integer) where p
    k = knotvector(P)
    return k[i]..k[i+p+1]
end

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)}|_I)=[k_{i},k_{i+p+1}] \cap I \qquad (I = [k_{1+p}, k_{l-p}])
```

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

julia> P = BSplineSpace{2}(k)
BSplineSpace{2, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> bsplinesupport_I(P,1)
2.5 .. 5.5

julia> bsplinesupport_I(P,2)
2.5 .. 8.0
```
"""
function bsplinesupport_I(P::BSplineSpace{p}, i::Integer) where p
    return bsplinesupport_R(P,i) ∩ domain(P)
end

@doc raw"""
Return the support of ``i``-th B-spline basis function.
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])
KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0])

julia> P = BSplineSpace{2}(k)
BSplineSpace{2, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]))

julia> bsplinesupport(P,1)
0.0 .. 5.5

julia> bsplinesupport(P,2)
1.5 .. 8.0
```
"""
const bsplinesupport = bsplinesupport_R

@doc raw"""
Internal methods for obtaining a B-spline space with one degree lower.
```math
\begin{aligned}
\mathcal{P}[p,k] &\mapsto \mathcal{P}[p-1,k] \\
D^r\mathcal{P}[p,k] &\mapsto D^{r-1}\mathcal{P}[p-1,k]
\end{aligned}
```
"""
_lower_R

_lower_R(P::BSplineSpace{p,T}) where {p,T} = BSplineSpace{p-1}(knotvector(P))

function _lower_I(P::BSplineSpace{p,T}) where {p,T}
    k = knotvector(P)
    l = length(k)
    return BSplineSpace{p-1}(view(k,2:l-1))
end

"""
Return an index of a interval in the domain of B-spline space

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]);


julia> P = BSplineSpace{2}(k);


julia> domain(P)
2.5 .. 9.0

julia> intervalindex(P,2.6)
1

julia> intervalindex(P,5.6)
2

julia> intervalindex(P,8.5)
3

julia> intervalindex(P,9.5)
3
```
"""
function intervalindex(P::BSplineSpace{p},t::Real) where p
    k = knotvector(P)
    l = length(k)
    v = view(_vec(k),2+p:l-p-1)
    return searchsortedlast(v,t)+1
end

"""
Expand B-spline space with given additional degree and knotvector.
This function is compatible with `issqsubset` (`⊑`)

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]);


julia> P = BSplineSpace{2}(k);


julia> P′ = expandspace_I(P, Val(1), KnotVector([6.0]))
BSplineSpace{3, Float64, KnotVector{Float64}}(KnotVector([0.0, 1.5, 2.5, 2.5, 5.5, 5.5, 6.0, 8.0, 8.0, 9.0, 9.0, 9.5, 10.0]))

julia> P ⊆ P′
false

julia> P ⊑ P′
true

julia> domain(P)
2.5 .. 9.0

julia> domain(P′)
2.5 .. 9.0
```
"""
function expandspace_I end

function expandspace_I(P::BSplineSpace{p,T}, ::Val{p₊}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,p₊,T}
    k = knotvector(P)
    I = domain(P)
    l = length(k)
    l₊ = length(k₊)
    if l₊ ≥ 1
        k₊[1] ∉ I && throw(DomainError(k₊, "input knot vector is out of domain."))
    end
    if l₊ ≥ 2
        k₊[end] ∉ I && throw(DomainError(k₊, "input knot vector is out of domain."))
    end
    k̂ = unique(view(k, 1+p:l-p))
    p′ = p + p₊
    k′ = k + p₊*k̂ + k₊
    P′ = BSplineSpace{p′}(k′)
    return P′
end

function expandspace_I(P::BSplineSpace{p,T}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,T}
    k = knotvector(P)
    I = domain(P)
    l₊ = length(k₊)
    if l₊ ≥ 1
        k₊[1] ∉ I && throw(DomainError(k₊, "input knot vector is out of domain."))
    end
    if l₊ ≥ 2
        k₊[end] ∉ I && throw(DomainError(k₊, "input knot vector is out of domain."))
    end
    k′ = k + k₊
    P′ = BSplineSpace{p}(k′)
    return P′
end

"""
Expand B-spline space with given additional degree and knotvector.
This function is compatible with `issubset` (`⊆`).

# Examples
```jldoctest
julia> k = KnotVector([0.0, 1.5, 2.5, 5.5, 8.0, 9.0, 9.5, 10.0]);


julia> P = BSplineSpace{2}(k);


julia> P′ = expandspace_R(P, Val(1), KnotVector([6.0]))
BSplineSpace{3, Float64, KnotVector{Float64}}(KnotVector([0.0, 0.0, 1.5, 1.5, 2.5, 2.5, 5.5, 5.5, 6.0, 8.0, 8.0, 9.0, 9.0, 9.5, 9.5, 10.0, 10.0]))

julia> P ⊆ P′
true

julia> P ⊑ P′
false

julia> domain(P)
2.5 .. 9.0

julia> domain(P′)
1.5 .. 9.5
```
"""
function expandspace_R end

function expandspace_R(P::BSplineSpace{p,T}, ::Val{p₊}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,p₊,T}
    k = knotvector(P)
    p′ = p + p₊
    k′ = k + p₊*k + k₊
    P′ = BSplineSpace{p′}(k′)
    return P′
end

function expandspace_R(P::BSplineSpace{p,T}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,T}
    k = knotvector(P)
    k′ = k + k₊
    P′ = BSplineSpace{p}(k′)
    return P′
end

"""
Expand B-spline space with given additional degree and knotvector.
The behavior of `expandspace` is same as `expandspace_I`.
"""
function expandspace(P::BSplineSpace{p,T}, _p₊::Val{p₊}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,p₊,T}
    return expandspace_I(P,_p₊,k₊)
end

function expandspace(P::BSplineSpace{p,T}, k₊::AbstractKnotVector=EmptyKnotVector{T}()) where {p,T}
    return expandspace_I(P,k₊)
end

function Base.hash(P::BSplineSpace{p}, h::UInt) where p
    k = knotvector(P)
    return hash(BSplineSpace{p}, hash(_vec(k), h))
end

function clamp!(P::BSplineSpace{p, T, <:KnotVector}) where {p, T}
    v = _vec(knotvector(P))
    v[1:p] .= v[p+1]
    v[end-p+1:end] .= v[end-p]
    return P
end
function clamp(P::BSplineSpace{p, T, <:KnotVector}) where {p, T}
    return clamp!(copy(P))
end
function clamp(P::BSplineSpace{p, T, <:UniformKnotVector}) where {p, T}
    k = knotvector(P)
    return clamp(BSplineSpace{p}(KnotVector(k)))
end

function isclamped(P::BSplineSpace{p}) where p
    k = knotvector(P)
    return k[1] == k[p+1] && k[end-p] == k[end]
end

function Base.:+(P1::BSplineSpace{p}, P2::BSplineSpace{p}) where {p}
    return BSplineSpace{p}(knotvector(P1) ∪ knotvector(P2))
end

function expand_domain(P::BSplineSpace{p,T1}, Δt::T2) where {p,T1,T2}
    U = promote_type(T1, T2)
    v = Vector{U}(_vec(knotvector(P)))
    v[1:p+1] .= v[p+1]-Δt
    v[end-p:end] .= v[end-p]+Δt
    return BSplineSpace{p}(KnotVector(v))
end
