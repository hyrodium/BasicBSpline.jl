# Change Basis between B-Spline Spaces
# See https://hackmd.io/lpA0D0ySTQ6Hq1CdaaONxQ for more information.

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊆ P^{\prime}``
"""
function _changebasis_R end

function _changebasis_R(P::AbstractBSplineSpace{0,T}, P′::AbstractBSplineSpace{p′,T′}) where {p′,T,T′}
    P ⊆ P′ || throw(DomainError((P,P′),"P ⊆ P′ should be hold."))
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    n = dim(P)
    n′ = dim(P′)
    A⁰ = U[bsplinesupport(P′, j) ⊆ bsplinesupport(P, i) for i in 1:n, j in 1:n′]
    A⁰[:, findall(_iszeros(P′))] .= zero(U)  # Strictly this should be NaN, but we use zero here to support Rational.
    return A⁰
end

function _changebasis_R(P::AbstractBSplineSpace{p,T}, P′::AbstractBSplineSpace{p′,T′}) where {p,p′,T,T′}
    P ⊆ P′ || throw(DomainError((P,P′),"P ⊆ P′ should be hold."))
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    n = dim(P)
    n′ = dim(P′)

    Aᵖ⁻¹ = _changebasis_R(_lower(P), _lower(P′))  # (n+1) × (n′+1) matrix
    Aᵖ = zeros(U, n, n′)  # n × n′ matrix

    Z = _iszeros(_lower(P′))
    W = findall(Z)
    K′ = [k′[i+p′] - k′[i] for i in 1:n′+1]
    K = [ifelse(k[i+p] ≠ k[i], 1 / (k[i+p] - k[i]), 0.0) for i in 1:n+1]
    Δ = (p / p′) * [K′[j] * (K[i] * Aᵖ⁻¹[i, j] - K[i+1] * Aᵖ⁻¹[i+1, j]) for i in 1:n, j in 1:n′+1]
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
            Ãᵖ[ȷ] .= zero(U)  # Strictly this should be NaN, but we use zero here to support Rational.
        end
    end
    for ȷ in 1:λ-1
        if Λ[ȷ] ≥ 2
            t = k′[W[ȷ]]
            for i in 1:n
                Ãᵖ[ȷ][i, end] = bsplinebasis₋₀(P,i,t)
            end
        end
    end
    for ȷ in 2:λ
        if Λ[ȷ] ≥ 2
            t = k′[W[ȷ-1]+p]
            for i in 1:n
                Ãᵖ[ȷ][i, 1] = bsplinebasis₊₀(P,i,t)
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
    return Aᵖ .* T[bsplinesupport(P′,j) ⊆ bsplinesupport(P,i) for i in 1:n, j in 1:n′]
end

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p_1,k_1)} = \sum_{j}A_{i,j}B_{(j,p_2,k_2)}
```

Assumption:
* ``P_1 ≃ P_2``
"""
function _changebasis_sim(P1::AbstractBSplineSpace{p,T1}, P2::AbstractBSplineSpace{p,T2}) where {p,T1,T2}
    P1 ≃ P2 || throw(DomainError((P1,P2),"P1 ≃ P2 should be hold."))
    U = StaticArrays.arithmetic_closure(promote_type(T1,T2))
    n = dim(P1)
    a, b = extrema(domain(P1))

    A = Matrix{U}(I, n, n)
    A1 = @MMatrix zeros(U,p,p)
    A2 = @MMatrix zeros(U,p,p)
    for r in 1:p
        A1[:,r] = bsplinebasis₊₀.(BSplineDerivativeSpace{r-1}(P1),1:p,a)
        A2[:,r] = bsplinebasis₊₀.(BSplineDerivativeSpace{r-1}(P2),1:p,a)
    end
    A[1:p, 1:p] = A1/A2
    for r in 1:p
        A1[:,r] = bsplinebasis₋₀.(BSplineDerivativeSpace{r-1}(P1),n-p+1:n,b)
        A2[:,r] = bsplinebasis₋₀.(BSplineDerivativeSpace{r-1}(P2),n-p+1:n,b)
    end
    A[n-p+1:n, n-p+1:n] = A1/A2
    return A
end

@doc raw"""
Return a coefficient matrix A which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊑ P^{\prime}``
"""
function _changebasis_I(P::AbstractBSplineSpace{p,T}, P′::AbstractBSplineSpace{p′,T′}) where {p,p′,T,T′}
    P ⊑ P′ || throw(DomainError((P,P′),"P ⊑ P′ should be hold."))
    k = knotvector(P)
    k′ = knotvector(P′)

    _P = BSplineSpace{p}(k[1+p:end-p] + p * KnotVector{T}(k[1+p], k[end-p]))
    _P′ = BSplineSpace{p′}(k′[1+p′:end-p′] + p′ * KnotVector{T′}(k′[1+p′], k′[end-p′]))
    _A = _changebasis_R(_P, _P′)
    Asim = _changebasis_sim(P, _P)
    Asim′ = _changebasis_sim(_P′, P′)
    A = Asim * _A * Asim′

    return A
end

## UniformBSplineSpace
function _changebasis_R(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    P ⊆ P′ || throw(DomainError((P,P′),"P ⊆ P′ should be hold."))
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    j = findfirst(==(k[1]), _vec(k′))
    A = zeros(StaticArrays.arithmetic_closure(T), n, n′)
    for i in 1:n
        A[i, j+r*(i-1):j+(r-1)*(p+1)+r*(i-1)] = block
    end
    return A
end
function _changebasis_I(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    P ⊑ P′ || throw(DomainError((P,P′),"P ⊑ P′ should be hold."))
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    A = zeros(StaticArrays.arithmetic_closure(T), n, n′)
    for i in 1:n
        a = r*i-(r-1)*(p+1)
        b = r*i
        rr = a:b
        if rr ⊆ 1:n′
            A[i,rr] = block
        else
            A[i,max(a,1):min(b,n′)] = block[max(a,1)-a+1:min(b,n′)-b+(r-1)*(p+1)+1]
        end
    end
    return A
end

## BSplineDerivativeSpace
function _changebasis_R(dP::BSplineDerivativeSpace{r,<:AbstractBSplineSpace{p}}, P′::AbstractBSplineSpace) where {r,p}
    dP ⊆ P′ || throw(DomainError((P,P′),"dP ⊆ P′ should be hold."))
    k = knotvector(dP)
    n = dim(dP)
    A = Matrix(I(n))
    for _r in 0:r-1
        _p = p - _r
        _A = zeros(n+_r,n+1+_r)
        for i in 1:n+_r
            _A[i,i] = _p/(k[i+_p]-k[i])
            _A[i,i+1] = -_p/(k[i+_p+1]-k[i+1])
        end
        A = A*_A
    end
    _P = BSplineSpace{degree(dP)}(k)
    A = A*_changebasis_R(_P, P′)
    return A
end
function _changebasis_R(dP::BSplineDerivativeSpace, dP′::BSplineDerivativeSpace{0})
    dP ⊆ dP′ || throw(DomainError((dP,dP′),"dP ⊆ dP′ should be hold."))
    P′ = bsplinespace(dP′)
    return _changebasis_R(dP, P′)
end
function _changebasis_R(dP::BSplineDerivativeSpace{r}, dP′::BSplineDerivativeSpace{r′}) where {r,r′}
    dP ⊆ dP′ || throw(DomainError((dP,dP′),"dP ⊆ dP′ should be hold."))
    if r > r′
        P = bsplinespace(dP)
        P′ = bsplinespace(dP′)
        _dP = BSplineDerivativeSpace{r-r′}(P)
        return _changebasis_R(_dP, P′)
    elseif r == r′
        P = bsplinespace(dP)
        P′ = bsplinespace(dP′)
        return _changebasis_R(P, P′)
    end
end
function _changebasis_R(P::AbstractBSplineSpace, dP′::BSplineDerivativeSpace{0})
    P′ = bsplinespace(dP′)
    return _changebasis_R(P, P′)
end

function changebasis(P::AbstractFunctionSpace, P′::AbstractFunctionSpace)
    P ⊆ P′ && return _changebasis_R(P, P′)
    P ⊑ P′ && return _changebasis_I(P, P′)
    throw(DomainError((P, P′),"𝒫[p,k] ⊆ 𝒫[p′,k′] or 𝒫[p,k] ⊑ 𝒫[p′,k′] must hold."))
end
