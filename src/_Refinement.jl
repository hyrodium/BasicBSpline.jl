# Refinement

@doc raw"""
Refinement of B-spline manifold with given B-spline spaces.
"""
function refinement(M::AbstractBSplineManifold, Ps′::Vector{<:AbstractBSplineSpace})
    Ps = collect(bsplinespaces(M))
    a = controlpoints(M)
    d = length(Ps)
    n = dim.(Ps)
    n′ = dim.(Ps′)

    A = changebasis.(Ps, Ps′)
    # TODO: general dimension
    if d == 1
        a′ = [sum(A[1][I₁, J₁] * a[I₁] for I₁ in 1:n[1]) for J₁ in 1:n′[1]]
    elseif d == 2
        a′ = [sum(A[1][I₁, J₁] * A[2][I₂, J₂] * a[I₁, I₂] for I₁ in 1:n[1], I₂ in 1:n[2]) for J₁ in 1:n′[1], J₂ in 1:n′[2]]
    elseif d == 3
        a′ = [
            sum(A[1][I₁, J₁] * A[2][I₂, J₂] * A[3][I₃, J₃] * a[I₁, I₂, I₃] for I₁ in 1:n[1], I₂ in 1:n[2], I₃ in 1:n[3])
            for J₁ in 1:n′[1], J₂ in 1:n′[2], J₃ in 1:n′[3]
        ]
    end
    return (typeof(M).name.wrapper)(Ps′, a′)
end


@doc raw"""
Refinement of B-spline manifold with additional degree and knots.
"""
function refinement(M::AbstractBSplineManifold; p₊::Union{Nothing,AbstractVector{<:Integer}} = nothing, k₊::Union{Nothing,Vector{<:Knots}} = nothing)
    Ps = collect(bsplinespaces(M))
    𝒂 = controlpoints(M)
    d = length(Ps)
    n = dim.(Ps)
    if isnothing(p₊)
        p₊ = zeros(Int, d)
    elseif length(Ps) ≠ length(p₊)
        throw(DimensionMismatch())
    end
    if isnothing(k₊)
        k₊ = zeros(Knots, d)
    elseif length(Ps) ≠ length(k₊)
        throw(DimensionMismatch())
    end

    Ps′ = [(P = Ps[i];
    p = degree(P);
    k = knots(P);
    k_unique = unique(k[1+p:end-p]);
    typeof(P).name.wrapper{p + p₊[i]}(k + p₊[i] * k_unique + k₊[i])) for i in eachindex(Ps)]

    return refinement(M, Ps′)
end
