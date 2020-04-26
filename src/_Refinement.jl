# Refinement
function BSplineCoefficient(P::BSplineSpace, P′::BSplineSpace)::Array{Float64,2}
    p = P.degree
    k = P.knots
    p′ = P′.degree
    k′ = P′.knots
    p₊ = p′-p
    if P ⊈ P′
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end

    if p == 0
        n=length(k)-1
        n′=length(k′)-p₊-1
        A⁰=Float64[BSplineSupport(j,𝒫(p₊,k′)) ⊆ BSplineSupport(i,𝒫(0,k)) for i ∈ 1:n, j ∈ 1:n′]
        A⁰[:,findall(iszeros(P′))].=NaN
        return A⁰
    end

    Aᵖ⁻¹=BSplineCoefficient(𝒫(p-1, k), 𝒫(p′-1, k′))
    n = dim(P)
    n′=dim(P′)
    Z = iszeros(𝒫(p′-1,k′))
    W = findall(Z)
    K′ = [k′[i+p′]-k′[i] for i ∈ 1:n′+1]
    K = [ifelse(k[i+p]≠k[i], 1/(k[i+p]-k[i]), 0.0) for i ∈ 1:n+1]
    Δ = (p/p′)*[K′[j]*(K[i]*Aᵖ⁻¹[i,j]-K[i+1]*Aᵖ⁻¹[i+1,j]) for i ∈ 1:n, j ∈ 1:n′+1]
    Aᵖ = zeros(n,n′)
    Aᵖ[:,1] = Δ[:,1]
    Aᵖ[:,n′] = -Δ[:,n′+1]

    if length(W) == 0
        Q = [1:n′]
    else
        Q = [1:W[1]-1,[W[i]:W[i+1]-1 for i ∈ 1:length(W)-1]...,W[end]:n′]
    end
    l = length(Q)
    L = length.(Q)
    Ãᵖ = [Aᵖ[:,q] for q ∈ Q]

    for ȷ ∈ 2:l-1
        if L[ȷ] == 1
            Ãᵖ[ȷ] .= NaN
        end
    end
    for ȷ ∈ 1:l-1
        if L[ȷ] ≥ 2
            t = k′[W[ȷ]]
            Ãᵖ[ȷ][:,end] = BSplineBasis₋₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 2:l
        if L[ȷ] ≥ 2
            t = k′[W[ȷ-1]+p]
            Ãᵖ[ȷ][:,1] = BSplineBasis₊₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 1:l
        if L[ȷ] ≥ 3
            r = Q[ȷ]
            A₊ = copy(Ãᵖ[ȷ])
            A₋ = copy(Ãᵖ[ȷ])
            for j ∈ 1:L[ȷ]-2
                A₊[:,j+1] = A₊[:,j]+Δ[:,j+r[1]]
                A₋[:,L[ȷ]-j] = A₋[:,L[ȷ]-j+1]-Δ[:,L[ȷ]-j+r[1]]
            end
            Ãᵖ[ȷ] = (A₊+A₋)/2
        end
    end
    Aᵖ = hcat(Ãᵖ...)
    return Aᵖ .* Float64[BSplineSupport(j,P′) ⊆ BSplineSupport(i,P) for i ∈ 1:n, j ∈ 1:n′]
end


@doc raw"""
Refinement of B-spline manifold.
"""
function Refinement(M::BSplineManifold, Ps′::Array{BSplineSpace,1})
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d̂ = size(𝒂)[end]
    n = dim.(Ps)
    n′ = dim.(Ps′)
    if prod(Ps .⊆ Ps′)
        A = BSplineCoefficient.(Ps,Ps′)
        𝒂′ = [sum(A[1][I₁,J₁]*A[2][I₂,J₂]*𝒂[I₁,I₂,i] for I₁ ∈ 1:n[1], I₂ ∈ 1:n[2]) for J₁ ∈ 1:n′[1], J₂ ∈ 1:n′[2], i ∈ 1:d̂]
        return BSplineManifold(Ps′, 𝒂′)
    else
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end
end

@doc raw"""
Refinement of B-spline manifold.
"""
function Refinement(M::BSplineManifold; p₊::Union{Nothing,Array{Int,1}}=nothing, k₊::Union{Nothing,Array{Knots,1}}=nothing)
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    d = length(Ps)
    d̂ = size(𝒂)[end]
    n = dim.(Ps)
    if p₊ == nothing
        p₊=zeros(Int,d)
    elseif length(Ps) ≠ length(p₊)
        error("dimension does not match")
    end
    if k₊ == nothing
        k₊=zeros(Knots,d)
    elseif length(Ps) ≠ length(k₊)
        error("dimension does not match")
    end

    Ps′ = BSplineSpace[]
    for i ∈ 1:length(Ps)
        P = Ps[i]
        p = P.degree
        k = P.knots
        push!(Ps′,𝒫(p+p₊[i], k+p₊[i]*unique(k)+k₊[i]))
    end

    return Refinement(M, Ps′)
end
