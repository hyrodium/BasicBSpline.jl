# Change Basis between B-Spline Spaces

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊆ P′``
"""
function _changebasis_R(P::BSplineSpace{p,T}, P′::BSplineSpace{p′,T})::Matrix{T} where {p,p′,T}
    k = knots(P)
    k′ = knots(P′)
    p₊ = p′ - p

    if p == 0
        n = length(k) - 1
        n′ = length(k′) - p₊ - 1
        A⁰ = Float64[bsplinesupport(BSplineSpace{p₊}(k′), j) ⊆ bsplinesupport(BSplineSpace{0}(k), i) for i in 1:n, j in 1:n′]
        A⁰[:, findall(iszeros(P′))] .= NaN
        return A⁰
    end

    Aᵖ⁻¹ = _changebasis_R(lower(P), lower(P′)) # (n+1) × (n′+1) matrix
    n = dim(P)
    n′ = dim(P′)
    Z = iszeros(lower(P′))
    W = findall(Z)
    K′ = [k′[i+p′] - k′[i] for i in 1:n′+1]
    K = [ifelse(k[i+p] ≠ k[i], 1 / (k[i+p] - k[i]), 0.0) for i in 1:n+1]
    Δ = (p / p′) * [K′[j] * (K[i] * Aᵖ⁻¹[i, j] - K[i+1] * Aᵖ⁻¹[i+1, j]) for i in 1:n, j in 1:n′+1]
    Aᵖ = zeros(n, n′)
    Aᵖ[:, 1] = Δ[:, 1]
    Aᵖ[:, n′] = -Δ[:, n′+1]

    # split Aᵖ for sub-block
    if length(W) == 0
        Q = [1:n′]
    else
        Q = [1:W[1]-1, [W[i]:W[i+1]-1 for i in 1:length(W)-1]..., W[end]:n′]
    end
    λ = length(Q)
    Λ = length.(Q)
    Ãᵖ = [Aᵖ[:, q] for q in Q]

    for ȷ in 2:λ-1
        if Λ[ȷ] == 1
            # if B(i,p′,k′) = 0
            Ãᵖ[ȷ] .= NaN
        end
    end
    for ȷ in 1:λ-1
        if Λ[ȷ] ≥ 2
            t = k′[W[ȷ]]
            for i in 1:n
                Ãᵖ[ȷ][i, end] = bsplinebasis₋₀(P, i, t)
            end
        end
    end
    for ȷ in 2:λ
        if Λ[ȷ] ≥ 2
            t = k′[W[ȷ-1]+p]
            for i in 1:n
                Ãᵖ[ȷ][i, 1] = bsplinebasis₊₀(P, i, t)
            end
        end
    end
    for ȷ in 1:λ
        if Λ[ȷ] ≥ 3
            r = Q[ȷ]
            A₊ = copy(Ãᵖ[ȷ])
            A₋ = copy(Ãᵖ[ȷ])
            for j in 1:Λ[ȷ]-2
                A₊[:, j+1] = A₊[:, j] + Δ[:, j+r[1]]
                A₋[:, Λ[ȷ]-j] = A₋[:, Λ[ȷ]-j+1] - Δ[:, Λ[ȷ]-j+r[1]]
            end
            Ãᵖ[ȷ] = (A₊ + A₋) / 2
        end
    end
    Aᵖ = hcat(Ãᵖ...) # n × n′ matrix
    return Aᵖ .* Float64[bsplinesupport(P′,j) ⊆ bsplinesupport(P,i) for i in 1:n, j in 1:n′]
end

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p_1,k_1)} = \sum_{j}A_{i,j}B_{(j,p_2,k_2)}
```

Assumption:
* ``P1 ⊑ P2``
* ``P2 ⊑ P1``
"""
function _changebasis_sim(P1::AbstractBSplineSpace, P2::AbstractBSplineSpace)
    # if P1 ⋢ P2
    #     error("P1 ⋢ P2")
    # end
    # if P2 ⋢ P1
    #     error("P2 ⋢ P1")
    # end
    n = dim(P1)
    p = degree(P1)
    v = (knots(P1).vector)[p+1:end-p]

    if length(v) > p
        A = Matrix{Float64}(I, n, n)
        # TODO: Fix below
        vvv = [v[1] * (p - i + 1) / (p + 1) + v[i+1] * (i) / (p + 1) for i in 1:p]
        A1 = [bsplinebasis₊₀(P1, i, t) for i in 1:p, t in vvv]
        A2 = [bsplinebasis₊₀(P2, i, t) for i in 1:p, t in vvv]
        A[1:p, 1:p] = A1 * inv(A2)
        vvv = [v[end-p+i-1] * (p - i + 1) / (p + 1) + v[end] * (i) / (p + 1) for i in 1:p]
        A1 = [bsplinebasis₋₀(P1, i, t) for i in n-p+1:n, t in vvv]
        A2 = [bsplinebasis₋₀(P2, i, t) for i in n-p+1:n, t in vvv]
        A[n-p+1:n, n-p+1:n] = A1 * inv(A2)
        # TODO: Fix above
    else
        # TODO: Fix below
        vvv = [v[1] * (n - i + 1) / (n + 1) + v[end] * (i) / (n + 1) for i in 1:n]
        A1 = [bsplinebasis₋₀(P1, i, t) for i in 1:n, t in vvv]
        A2 = [bsplinebasis₋₀(P2, i, t) for i in 1:n, t in vvv]
        A = A1 * inv(A2)
        # TODO: Fix above
    end
    return A
end

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊑ P′``
"""
function _changebasis_I(P::BSplineSpace{p,T}, P′::BSplineSpace{p′,T})::Matrix{T} where {p, p′, T}
    I = bsplineunity(P)
    k = knots(P)
    k′ = knots(P′)
    p₊ = p′ - p

    _P = BSplineSpace{p}(k[1+p:end-p] + p * Knots(k[1+p], k[end-p]))
    # if dim(_P) ≠ dim(P)
    #     error("dim(_P)≠dim(P)")
    # end
    _P′ = BSplineSpace{p′}(k′[1+p′:end-p′] + p′ * Knots(k′[1+p′], k′[end-p′]))
    # if dim(_P′) ≠ dim(P′)
    #     error("dim(_P′)≠dim(P′)")
    # end
    _A = _changebasis_R(_P, _P′)
    Asim = _changebasis_sim(P, _P)
    Asim′ = _changebasis_sim(_P′, P′)
    A = Asim * _A * Asim′

    return A
end

function changebasis(P::BSplineSpace, P′::BSplineSpace)
    if P ⊆ P′
        return _changebasis_R(P, P′)
    elseif P ⊑ P′
        return _changebasis_I(P, P′)
    else
        throw(DomainError((P, P′),"𝒫[p,k] ⊆ 𝒫[p′,k′] or 𝒫[p,k] ⊑ 𝒫[p′,k′] must hold."))
    end
end

function changebasis(P::AbstractBSplineSpace, P′::AbstractBSplineSpace)
    changebasis(BSplineSpace(P),BSplineSpace(P′))
end
