# Knots
@doc raw"""
Construct knot vector from given array.
```math
k=(k_1,\dots,k_l)
```
"""
struct Knots
    vector :: Array{Float64,1}
    function Knots(vector::AbstractArray{T,1} where T<:Real)
        return new(sort(convert(Array{Float64,1},vector)))
    end
    function Knots(vector::Array{Any,1})
        if isempty(vector)
            return new(Float64[])
        else
            error("The elements of given vector must be real number.")
        end
    end
    function Knots(knot::Real...)
        return Knots(collect(knot))
    end
end

Base.zero(::Type{Knots}) = Knots([])
Base. ==(k₁::Knots, k₂::Knots) = (k₁.vector==k₂.vector)
Base.:+(k₁::Knots, k₂::Knots) = Knots(sort([k₁.vector...,k₂.vector...]))
Base.:*(p₊::Int, k::Knots) = (
        if p₊ == 0
            zero(Knots)
        elseif p₊ > 0
            sum(k for _ ∈ 1:p₊)
        else
            error("Polynominal degree p₊ must be non-negative.")
        end
    )

Base.in(r::Real, k::Knots) = in(r,k.vector)
Base.getindex(k::Knots, i::Int) = k.vector[i]
Base.getindex(k::Knots, v::AbstractArray{Int64,1}) = Knots(k.vector[v])
Base.length(k::Knots) = length(k.vector)
♯(k::Knots) = length(k::Knots)
Base.firstindex(k) = 1
Base.lastindex(k) = length(k)
Base.unique(k::Knots) = Knots(unique(k.vector))

function Base.:⊆(k::Knots, k′::Knots)
    K′ = copy(k′.vector)
    for kᵢ ∈ k.vector
        i = findfirst(x -> x == kᵢ,K′)
        if i isa Nothing
            return false
        end
        deleteat!(K′,i)
    end
    return true
end

𝔫(k::Knots, t::Real) = count(s -> (t == s), k.vector)
