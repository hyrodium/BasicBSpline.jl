# Knots

@doc raw"""
Construct knot vector from given array.
```math
k=(k_1,\dots,k_l)
```
"""
struct Knots
    vector::Array{Float64,1}
    function Knots(vector::AbstractVector{<:Real})
        return new(sort(convert(Vector{Float64}, vector)))
    end
end
function Knots(vector::Array{Any,1})
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
        error("Polynominal degree p₊ must be non-negative.")
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

function Base.:⊆(k::Knots, k′::Knots)
    K′ = copy(k′.vector)
    for kᵢ in k.vector
        i = findfirst(x -> x == kᵢ, K′)
        if i isa Nothing
            return false
        end
        deleteat!(K′, i)
    end
    return true
end

𝔫(k::Knots, t::Real) = count(==(t), k.vector)

function _knotindex₊₀(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> k[i]≤t<k[i+1], 1:length(k)-1)
end

function _knotindex₋₀(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> k[i]<t≤k[i+1], 1:length(k)-1)
end

function _knotindex(k::Union{Knots, AbstractVector{<:Real}}, t::Real)
    return findfirst(i -> (k[i]≤t<k[i+1])|(k[i]<t==k[i+1]==k[end]), 1:length(k)-1)
end
