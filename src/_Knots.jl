# Knots

@doc raw"""
Construct knot vector from given array.
```math
k=(k_1,\dots,k_l)
```

# Examples
```jldoctest; setup = :(using BasicBSpline)
julia> k = Knots(1:3)
Knots([1.0, 2.0, 3.0])
```
"""
struct Knots
    vector::Vector{Float64}
    function Knots(vector::AbstractVector{<:Real})
        return new(sort(convert(Vector{Float64}, vector)))
    end
end
function Knots(vector::AbstractVector)
    if isempty(vector)
        return Knots(Float64[])
    else
        return Knots(convert(Vector{Float64}, vector))
    end
end
function Knots(knot::Real...)
    return Knots(collect(knot))
end

Base.zero(::Type{Knots}) = Knots(Float64[])
Base.:(==)(k₁::Knots, k₂::Knots) = (k₁.vector == k₂.vector)
Base.:+(k₁::Knots, k₂::Knots) = Knots(sort([k₁.vector..., k₂.vector...]))
function Base.:*(p₊::Integer, k::Knots)
    if p₊ == 0
        return zero(Knots)
    elseif p₊ > 0
        return sum(k for _ in 1:p₊)
    else
        throw(DomainError(p₊, "additional degree of polynominal must be non-negative"))
    end
end

Base.in(r::Real, k::Knots) = in(r, k.vector)
Base.getindex(k::Knots, i::Integer) = k.vector[i]
Base.getindex(k::Knots, v::AbstractVector{<:Integer}) = Knots(k.vector[v])
Base.length(k::Knots) = length(k.vector)
♯(k::Knots) = length(k::Knots)
Base.firstindex(k::Knots) = 1
Base.lastindex(k::Knots) = length(k)
Base.unique(k::Knots) = Knots(unique(k.vector))
Base.iterate(k::Knots) = iterate(k.vector)
Base.iterate(k::Knots, i::Integer) = iterate(k.vector, i)

@doc raw"""
Check a inclusive relation ship ``k\subset k'``.
```math
(1,2,3) \subseteq (1,2,3,4)
(1,2,5) \nsubseteq (1,2,3,4)
```

# Examples
```jldoctest; setup = :(using BasicBSpline)
julia> k1, k2 = Knots(1:3), Knots(-2,1,3,2,4);

julia> k1 ⊆ k2
true
```
"""
function Base.issubset(k::Knots, k′::Knots)
    K′ = copy(k′.vector)
    for kᵢ in k.vector
        i = findfirst(==(kᵢ), K′)
        if isnothing(i)
            return false
        end
        deleteat!(K′, i)
    end
    return true
end

𝔫(k::Knots, t::Real) = count(==(t), k.vector)

"""
Find an index ``i`` such that ``k_{i} ≤ t < k_{i+1}``.
"""
function _knotindex₊₀(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> k[i]≤t<k[i+1], 1:length(k)-1)
end

"""
Find an index ``i`` such that ``k_{i} < t ≤ k_{i+1}``.
"""
function _knotindex₋₀(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> k[i]<t≤k[i+1], 1:length(k)-1)
end

"""
Find an index ``i`` such that ``k_{i} ≤ t < k_{i+1}`` or ``k_{i} < t = k_{i+1} = k_{\text{end}}``.
"""
function _knotindex(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> (k[i]≤t<k[i+1])|(k[i]<t==k[i+1]==k[end]), 1:length(k)-1)
end
