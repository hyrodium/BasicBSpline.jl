# Change Basis between B-Spline Spaces
# See https://hackmd.io/lpA0D0ySTQ6Hq1CdaaONxQ for more information.

@doc raw"""
Return a coefficient matrix ``A`` which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊆ P^{\prime}``
"""
function changebasis_R(P::AbstractFunctionSpace, P′::AbstractFunctionSpace)
    P ⊆ P′ || throw(DomainError((P,P′),"P ⊆ P′ should be hold."))
    return _changebasis_R(P, P′)
end

function _changebasis_R(P::BSplineSpace{p,T}, P′::BSplineSpace{p′,T′}) where {p,T,p′,T′}
    _P = BSplineSpace{p,T,KnotVector{T}}(P)
    _P′ = BSplineSpace{p′,T′,KnotVector{T′}}(P′)
    return _changebasis_R(_P,_P′)
end

function _changebasis_R(P::BSplineSpace{0,T,KnotVector{T}}, P′::BSplineSpace{p′,T′,KnotVector{T′}}) where {p′,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T, T′))
    n = dim(P)
    n′ = dim(P′)
    n′_exact = exactdim_R(P′)
    k = knotvector(P)
    k′ = knotvector(P′)
    I = Vector{Int32}(undef, n′_exact)
    J = Vector{Int32}(undef, n′_exact)
    s = 1
    local j_begin
    j_end = 0
    for i in 1:n
        isdegenerate(P, i) && continue
        for j in (j_end+1):n′
            k′[j] == k[i] && (j_begin = j; break)
        end
        for j in j_begin:n′
            k′[j+p′+1] == k[i+1] && (j_end = j + p′; break)
        end
        for j in j_begin:j_end
            isdegenerate(P′, j) && continue
            I[s] = i
            J[s] = j
            s += 1
        end
    end
    A⁰ = sparse(view(I,1:s-1), view(J,1:s-1), fill(one(U), s-1), n, n′)
    return A⁰
end

function _find_j_begin_end_R(P::BSplineSpace{p}, P′::BSplineSpace{p′}, i, j_begin, j_end):: Tuple{Int, Int} where {p, p′}
    k = knotvector(P)
    k′ = knotvector(P′)
    n′ = dim(P′)
    Pi = BSplineSpace{p}(view(k, i:i+p+1))
    j_end::Int = findnext(j->Pi ⊆ BSplineSpace{p′}(view(k′, j_begin:j+p′+1)), 1:n′, j_end)
    j_begin::Int = findprev(j->Pi ⊆ BSplineSpace{p′}(view(k′, j:j_end+p′+1)), 1:n′, j_end)
    return j_begin, j_end
end

function _ΔAᵖ_R(Aᵖ⁻¹::AbstractMatrix, K::AbstractVector, K′::AbstractVector, i::Integer, j::Integer)
    return K′[j] * (K[i] * Aᵖ⁻¹[i, j] - K[i+1] * Aᵖ⁻¹[i+1, j])
end

function _changebasis_R(P::BSplineSpace{p,T,KnotVector{T}}, P′::BSplineSpace{p′,T′,KnotVector{T′}}) where {p,p′,T,T′}
    #=
    Example matrix: n=3, n′=4

    Aᵖ⁻¹₁₁      Aᵖ⁻¹₁₂      Aᵖ⁻¹₁₃      Aᵖ⁻¹₁₄      Aᵖ⁻¹₁₅
      ├    Aᵖ₁₁   ┼    Aᵖ₁₂   ┼    Aᵖ₁₃   ┼    Aᵖ₁₄   ┤
    Aᵖ⁻¹₂₁      Aᵖ⁻¹₂₂      Aᵖ⁻¹₂₃      Aᵖ⁻¹₂₄      Aᵖ⁻¹₂₅
      ├    Aᵖ₂₁   ┼    Aᵖ₂₂   ┼    Aᵖ₂₃   ┼    Aᵖ₂₄   ┤
    Aᵖ⁻¹₃₁      Aᵖ⁻¹₃₂      Aᵖ⁻¹₃₃      Aᵖ⁻¹₃₄      Aᵖ⁻¹₃₅
      ├    Aᵖ₃₁   ┼    Aᵖ₃₂   ┼    Aᵖ₃₃   ┼    Aᵖ₃₄   ┤
    Aᵖ⁻¹₄₁      Aᵖ⁻¹₄₂      Aᵖ⁻¹₄₃      Aᵖ⁻¹₄₄      Aᵖ⁻¹₄₅

    === i-rule ===
    - isdegenerate_R(P, i) ⇒ Aᵖᵢⱼ = 0

    === j-rules for fixed i ===
    - Rule-0: B₍ᵢ,ₚ,ₖ₎ ∈ span₍ᵤ₌ₛ…ₜ₎(B₍ᵤ,ₚ′,ₖ′₎), s≤j≤t ⇒ Aᵖᵢⱼ ≠ 0
              supp(B₍ⱼ,ₚ′,ₖ′₎) ⊈ supp(B₍ᵢ,ₚ,ₖ₎) ⇒ Aᵖᵢⱼ = 0 is almost equivalent (if there are not duplicated knots).
    - Rule-1: isdegenerate_R(P′, j) ⇒ Aᵖᵢⱼ = 0
              Ideally this should be NaN, but we need to support `Rational` which does not include `NaN`.
    - Rule-2: k′ⱼ   = k′ⱼ₊ₚ′ ⇒ B₍ᵢ,ₚ,ₖ₎(t) = ∑ⱼAᵖᵢⱼB₍ⱼ,ₚ′,ₖ′₎(t)  (t → k′ⱼ   + 0)
    - Rule-3: k′ⱼ₊₁ = k′ⱼ₊ₚ′ ⇒ B₍ᵢ,ₚ,ₖ₎(t) = ∑ⱼAᵖᵢⱼB₍ⱼ,ₚ′,ₖ′₎(t)  (t → k′ⱼ₊₁ - 0)
    - Rule-6: Aᵖᵢⱼ = Aᵖᵢⱼ₋₁ + (recursive formula based on Aᵖ⁻¹)
    - Rule-7: Aᵖᵢⱼ = Aᵖᵢⱼ₊₁ - (recursive formula based on Aᵖ⁻¹)
    =#
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    n = dim(P)
    n′ = dim(P′)
    K′ = [k′[j+p′] - k′[j] for j in 1:n′+1]
    K = U[ifelse(k[i+p] ≠ k[i], U(1 / (k[i+p] - k[i])), zero(U)) for i in 1:n+1]
    Aᵖ⁻¹ = _changebasis_R(_lower_R(P), _lower_R(P′))  # (n+1) × (n′+1) matrix
    n_nonzero = exactdim_R(P′)*(p+1)  # This is a upper bound of the number of non-zero elements of Aᵖ (rough estimation).
    I = Vector{Int32}(undef, n_nonzero)
    J = Vector{Int32}(undef, n_nonzero)
    V = Vector{U}(undef, n_nonzero)
    s = 1
    j_begin = 1
    j_end = 1
    for i in 1:n
        # Skip for degenerated basis
        isdegenerate_R(P,i) && continue

        # The following indices `j*` have the following relationships.
        #=

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
            |----*----------------*------>|
                 j_prev  j_mid    j_next
                 |       |        |
                  |--j₊->||<-j₋--|
                 266666666777777773
                 366666666777777773

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
           *|---------------*------------>|
           j_prev  j_mid    j_next
           |       |        |
            |--j₊->||<-j₋--|
           066666666777777773

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
            |-------------*-------------->|*
                          j_prev  j_mid    j_next
                          |       |        |
                           |--j₊->||<-j₋--|
                          266666666777777770
                          366666666777777770

         1                                    n′
         *-----------------------------*----->|
         j_begin                       j_end
         *----------------*----------->|
        j_prev  j_mid    j_next
        |       |        |
         |--j₊->||<-j₋--|
        066666666777777773

         1                                    n′
         |------*---------------------------->*
                j_begin                       j_end
                |-------------*--------------->*
                              j_prev  j_mid    j_next
                              |       |        |
                               |--j₊->||<-j₋--|
                              266666666777777770
                              366666666777777770

        =#

        # Precalculate the range of j
        j_begin, j_end = _find_j_begin_end_R(P, P′, i, j_begin, j_end)
        j_range = j_begin:j_end

        # Rule-0: outside of j_range
        j_prev = j_begin-1
        Aᵖᵢⱼ_prev = zero(U)
        for j_next in j_range
            # Rule-1: zero
            isdegenerate_R(P′,j_next) && continue
            # Rule-2: right limit
            if k′[j_next] == k′[j_next+p′]
                I[s], J[s] = i, j_next
                V[s] = Aᵖᵢⱼ_prev = bsplinebasis₊₀(P,i,k′[j_next+1])
                s += 1
                j_prev = j_next
            # Rule-3: left limit (or both limit)
            elseif k′[j_next+1] == k′[j_next+p′]
                j_mid = (j_prev + j_next) ÷ 2
                # Rule-6: right recursion
                for j₊ in (j_prev+1):j_mid
                    I[s], J[s] = i, j₊
                    V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * _ΔAᵖ_R(Aᵖ⁻¹,K,K′,i,j₊) / p′
                    s += 1
                end
                I[s], J[s] = i, j_next
                V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_next = bsplinebasis₋₀(P,i,k′[j_next+1])
                s += 1
                # Rule-7: left recursion
                for j₋ in reverse((j_mid+1):(j_next-1))
                    I[s], J[s] = i, j₋
                    V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * _ΔAᵖ_R(Aᵖ⁻¹,K,K′,i,j₋+1) / p′
                    s += 1
                end
                j_prev = j_next
            end
        end
        # Rule-0: outside of j_range
        j_next = j_end + 1
        j_mid = (j_prev + j_next) ÷ 2
        Aᵖᵢⱼ_next = zero(U)
        # Rule-6: right recursion
        for j₊ in (j_prev+1):j_mid
            I[s], J[s] = i, j₊
            V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * _ΔAᵖ_R(Aᵖ⁻¹,K,K′,i,j₊) / p′
            s += 1
        end
        # Rule-7: left recursion
        for j₋ in reverse((j_mid+1):(j_next-1))
            I[s], J[s] = i, j₋
            V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * _ΔAᵖ_R(Aᵖ⁻¹,K,K′,i,j₋+1) / p′
            s += 1
        end
    end

    Aᵖ = sparse(view(I,1:s-1), view(J,1:s-1), view(V,1:s-1), n, n′)
    return Aᵖ
end

@doc raw"""
Return a coefficient matrix ``A`` which satisfy
```math
B_{(i,p_1,k_1)} = \sum_{j}A_{i,j}B_{(j,p_2,k_2)}
```

Assumption:
* ``P_1 ≃ P_2``
"""
function _changebasis_sim(P1::BSplineSpace{p,T1}, P2::BSplineSpace{p,T2}) where {p,T1,T2}
    P1 ≃ P2 || throw(DomainError((P1,P2),"P1 ≃ P2 should be hold."))
    U = StaticArrays.arithmetic_closure(promote_type(T1,T2))
    k1 = knotvector(P1)
    k2 = knotvector(P2)
    n = dim(P1)     # == dim(P2)
    l = length(k1)  # == length(k2)

    A = SparseMatrixCSC{U,Int32}(I,n,n)

    A1 = _derivatives_at_left(P1)
    A2 = _derivatives_at_left(P2)
    A12 = A1/A2
    Δ = p
    for i in reverse(1:p)
        k1[i] ≠ k2[i] && break
        Δ -= 1
    end
    # A12 must be lower-triangular
    for i in 1:Δ, j in 1:i
        A[i, j] = A12[i,j]
    end

    A1 = _derivatives_at_right(P1)
    A2 = _derivatives_at_right(P2)
    A12 = A1/A2
    Δ = p
    for i in l-p+1:l
        k1[i] ≠ k2[i] && break
        Δ -= 1
    end
    # A12 must be upper-triangular
    for i in 1:Δ, j in i:Δ
        A[n-Δ+i, n-Δ+j] = A12[i+p-Δ,j+p-Δ]
    end

    return A
end

@generated function _derivatives_at_left(P::BSplineSpace{p,T}) where {p, T}
    args = [:(pop(bsplinebasisall(BSplineDerivativeSpace{$(r)}(P),1,a))) for r in 0:p-1]
    quote
        a, _ = extrema(domain(P))
        $(Expr(:call, :hcat, args...))
    end
end
function _derivatives_at_left(::BSplineSpace{0,T}) where {T}
    U = StaticArrays.arithmetic_closure(T)
    SMatrix{0,0,U}()
end

@generated function _derivatives_at_right(P::BSplineSpace{p,T}) where {p, T}
    args = [:(popfirst(bsplinebasisall(BSplineDerivativeSpace{$(r)}(P),n-p,b))) for r in 0:p-1]
    quote
        n = dim(P)
        _, b = extrema(domain(P))
        $(Expr(:call, :hcat, args...))
    end
end
function _derivatives_at_right(::BSplineSpace{0,T}) where {T}
    U = StaticArrays.arithmetic_closure(T)
    SMatrix{0,0,U}()
end

@doc raw"""
Return a coefficient matrix ``A`` which satisfy
```math
B_{(i,p,k)} = \sum_{j}A_{i,j}B_{(j,p',k')}
```

Assumption:
* ``P ⊑ P^{\prime}``
"""
function changebasis_I(P::AbstractFunctionSpace, P′::AbstractFunctionSpace)
    P ⊑ P′ || throw(DomainError((P,P′),"P ⊑ P′ should be hold."))
    return _changebasis_I(P, P′)
end

function _changebasis_I(P::BSplineSpace{0,T,<:AbstractKnotVector{T}}, P′::BSplineSpace{p′,T′,<:AbstractKnotVector{T′}}) where {p′,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T, T′))
    n = dim(P)
    n′ = dim(P′)
    n′_exact = exactdim_R(P′)
    k = knotvector(P)
    k′ = knotvector(P′)
    I = Vector{Int32}(undef, n′_exact)
    J = Vector{Int32}(undef, n′_exact)
    s = 1
    local j_begin
    j_end = 0
    for i in 1:n
        isdegenerate(P, i) && continue
        for j in (j_end+1):n′
            k′[j] ≤ k[i] && (j_begin = j; break)
        end
        for j in j_begin:n′
            k′[j+p′+1] ≥ k[i+1] && (j_end = j + p′; break)
        end
        for j in j_begin:j_end
            isdegenerate(P′, j) && continue
            I[s] = i
            J[s] = j
            s += 1
        end
    end
    A⁰ = sparse(view(I,1:s-1), view(J,1:s-1), fill(one(U), s-1), n, n′)
    return A⁰
end

function _changebasis_I_old(P::BSplineSpace{p,T}, P′::BSplineSpace{p′,T′}) where {p,p′,T,T′}
    k = knotvector(P)
    k′ = knotvector(P′)

    _P = BSplineSpace{p}(k[1+p:end-p] + p * KnotVector{T}([k[1+p], k[end-p]]))
    _P′ = BSplineSpace{p′}(k′[1+p′:end-p′] + p′ * KnotVector{T′}([k′[1+p′], k′[end-p′]]))
    _A = _changebasis_R(_P, _P′)
    Asim = _changebasis_sim(P, _P)
    Asim′ = _changebasis_sim(_P′, P′)
    A = Asim * _A * Asim′

    return A
end

function _find_j_begin_end_I(P::BSplineSpace{p}, P′::BSplineSpace{p′}, i, j_begin, j_end):: Tuple{Int, Int} where {p, p′}
    # TODO: avoid `_changebasis_I_old`
    Aᵖ_old = _changebasis_I_old(P,P′)
    j_begin = findfirst(e->abs(e)>1e-14, Aᵖ_old[i, :])
    j_end = findlast(e->abs(e)>1e-14, Aᵖ_old[i, :])
    return j_begin, j_end
end

function _ΔAᵖ_I(Aᵖ⁻¹::AbstractMatrix, K::AbstractVector, K′::AbstractVector, i::Integer, j::Integer)
    n = length(K)-1
    if i == 1
        return - K′[j] * K[i+1] * Aᵖ⁻¹[i, j-1]
    elseif i == n
        return K′[j] * K[i] * Aᵖ⁻¹[i-1, j-1]
    else
        return K′[j] * (K[i] * Aᵖ⁻¹[i-1, j-1] - K[i+1] * Aᵖ⁻¹[i, j-1])
    end
end

function _changebasis_I(P::BSplineSpace{p,T,<:AbstractKnotVector{T}}, P′::BSplineSpace{p′,T′,<:AbstractKnotVector{T′}}) where {p,p′,T,T′}
    #=
    Example matrix: n=4, n′=5

    Aᵖ₁₁   ┬    Aᵖ₁₂   ┬    Aᵖ₁₃   ┬    Aᵖ₁₄   ┬    Aᵖ₁₅
         Aᵖ⁻¹₁₁      Aᵖ⁻¹₁₂      Aᵖ⁻¹₁₃      Aᵖ⁻¹₁₄
    Aᵖ₂₁   ┼    Aᵖ₂₂   ┼    Aᵖ₂₃   ┼    Aᵖ₂₄   ┼    Aᵖ₂₅
         Aᵖ⁻¹₂₁      Aᵖ⁻¹₂₂      Aᵖ⁻¹₂₃      Aᵖ⁻¹₂₄
    Aᵖ₃₁   ┼    Aᵖ₃₂   ┼    Aᵖ₃₃   ┼    Aᵖ₃₄   ┼    Aᵖ₃₅
         Aᵖ⁻¹₃₁      Aᵖ⁻¹₃₂      Aᵖ⁻¹₃₃      Aᵖ⁻¹₃₄
    Aᵖ₄₁   ┴    Aᵖ₄₂   ┴    Aᵖ₄₃   ┴    Aᵖ₄₄   ┴    Aᵖ₄₅

    === i-rule ===
    - isdegenerate_I(P, i) ⇒ Aᵖᵢⱼ = 0

    === j-rules for fixed i ===
    - Rule-0: B₍ᵢ,ₚ,ₖ₎ ∈ span₍ᵤ₌ₛ…ₜ₎(B₍ᵤ,ₚ′,ₖ′₎), s≤j≤t ⇒ Aᵖᵢⱼ ≠ 0 (on B-spline space domain)
              supp(B₍ⱼ,ₚ′,ₖ′₎) ⊈ supp(B₍ᵢ,ₚ,ₖ₎) ⇒ Aᵖᵢⱼ = 0 is almost equivalent (if there are not duplicated knots).
    - Rule-1: isdegenerate_I(P′, j) ⇒ Aᵖᵢⱼ = 0
              Ideally this should be NaN, but we need to support `Rational` which does not include `NaN`.
    - Rule-2: k′ⱼ   = k′ⱼ₊ₚ′ ⇒ B₍ᵢ,ₚ,ₖ₎(t) = ∑ⱼAᵖᵢⱼB₍ⱼ,ₚ′,ₖ′₎(t)  (t → k′ⱼ   + 0)
    - Rule-3: k′ⱼ₊₁ = k′ⱼ₊ₚ′ ⇒ B₍ᵢ,ₚ,ₖ₎(t) = ∑ⱼAᵖᵢⱼB₍ⱼ,ₚ′,ₖ′₎(t)  (t → k′ⱼ₊₁ - 0)
    - Rule-6: Aᵖᵢⱼ = Aᵖᵢⱼ₋₁ + (recursive formula based on Aᵖ⁻¹)
    - Rule-7: Aᵖᵢⱼ = Aᵖᵢⱼ₊₁ - (recursive formula based on Aᵖ⁻¹)
    - Rule-8: Aᵖᵢⱼ = Aᵖᵢⱼ₋₁ + (recursive formula based on Aᵖ⁻¹) with postprocess Δ
    =#
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    n = dim(P)
    n′ = dim(P′)
    K′ = [k′[j+p′] - k′[j] for j in 1:n′+1]
    K = U[ifelse(k[i+p] ≠ k[i], U(1 / (k[i+p] - k[i])), zero(U)) for i in 1:n+1]
    Aᵖ⁻¹ = _changebasis_I(_lower_I(P), _lower_I(P′))  # (n-1) × (n′-1) matrix
    n_nonzero = exactdim_I(P′)*(p+1)  # This is a upper bound of the number of non-zero elements of Aᵖ (rough estimation).
    I = Vector{Int32}(undef, n_nonzero)
    J = Vector{Int32}(undef, n_nonzero)
    V = Vector{U}(undef, n_nonzero)
    s = 1
    j_begin = 1
    j_end = 1
    for i in 1:n
        # Skip for degenerated basis
        isdegenerate_I(P,i) && continue

        # The following indices `j*` have the following relationships.
        #=

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
            |----*----------------*------>|
                 j_prev  j_mid    j_next
                 |       |        |
                  |--j₊->||<-j₋--|
                 266666666777777773
                 366666666777777773

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
           *|---------------*------------>|
           j_prev  j_mid    j_next
           |       |        |
            |--j₊->||<-j₋--|
           066666666777777773

         1                                    n′
         |--*-----------------------------*-->|
            j_begin                       j_end
            |-------------*-------------->|*
                          j_prev  j_mid    j_next
                          |       |        |
                           |--j₊->||<-j₋--|
                          266666666777777770
                          366666666777777770

         1                                    n′
         *-----------------------------*----->|
         j_begin                       j_end
         *---------------*------------>|
        j_prev=j_mid     j_next
        |                |
         |<-----j₋------|
        077777777777777773

         1                                    n′
         |------*---------------------------->*
                j_begin                       j_end
                |-------------*--------------->*
                              j_prev   j_mid+1=j_next
                              |                |
                               |------j₊----->|
                              266666666666666660
                              366666666666666660

         1                                    n′
         *----------------------------------->*
         j_begin                              j_end
         *------------------------------------>*
        j_prev                                 j_next
        |                                      |
         |-----------------j₊---------------->|
        0888888888888888888888888888888888888880

        =#

        # Precalculate the range of j
        j_begin, j_end = _find_j_begin_end_I(P, P′, i, j_begin, j_end)
        j_range = j_begin:j_end

        # Rule-0: outside of j_range
        j_prev = j_begin-1
        Aᵖᵢⱼ_prev = zero(U)
        for j_next in j_range
            # Rule-1: zero
            isdegenerate_I(P′,j_next) && continue
            # Rule-2: right limit
            if k′[j_next] == k′[j_next+p′]
                I[s], J[s] = i, j_next
                V[s] = Aᵖᵢⱼ_prev = bsplinebasis₊₀(P,i,k′[j_next+1])
                s += 1
                j_prev = j_next
            # Rule-3: left limit (or both limit)
            elseif k′[j_next+1] == k′[j_next+p′]
                j_mid = (j_prev == 0 ? j_prev : (j_prev + j_next) ÷ 2)
                # Rule-6: right recursion
                for j₊ in (j_prev+1):j_mid
                    I[s], J[s] = i, j₊
                    V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * _ΔAᵖ_I(Aᵖ⁻¹,K,K′,i,j₊) / p′
                    s += 1
                end
                I[s], J[s] = i, j_next
                V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_next = bsplinebasis₋₀I(P,i,k′[j_next+1])
                s += 1
                # Rule-7: left recursion
                for j₋ in reverse((j_mid+1):(j_next-1))
                    I[s], J[s] = i, j₋
                    V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * _ΔAᵖ_I(Aᵖ⁻¹,K,K′,i,j₋+1) / p′
                    s += 1
                end
                j_prev = j_next
            end
        end
        # Rule-0: outside of j_range
        j_next = j_end + 1
        if j_next == n′+1
            j_mid = j_next - 1
            if j_prev == 0
                # Rule-8: right recursion with postprocess
                # We can't find Aᵖᵢ₁ or Aᵖᵢₙ′ directly (yet!), so we need Δ-shift.
                I[s], J[s] = i, 1
                V[s] = Aᵖᵢⱼ_prev = zero(U)
                s += 1
                for j₊ in 2:n′
                    I[s], J[s] = i, j₊
                    V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * _ΔAᵖ_I(Aᵖ⁻¹,K,K′,i,j₊) / p′
                    s += 1
                end
                t_mid = (maximum(domain(P))+minimum(domain(P))) / (2*one(U))
                Δ = bsplinebasis(P, i, t_mid) - dot(view(V, s-n′:s-1), bsplinebasis.(P′, 1:n′, t_mid))
                V[s-n′:s-1] .+= Δ
                continue
            end
        elseif j_prev == 0
            j_mid = j_prev
        else
            j_mid = (j_prev + j_next) ÷ 2
        end
        Aᵖᵢⱼ_next = zero(U)
        # Rule-6: right recursion
        for j₊ in (j_prev+1):j_mid
            I[s], J[s] = i, j₊
            V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * _ΔAᵖ_I(Aᵖ⁻¹,K,K′,i,j₊) / p′
            s += 1
        end
        # Rule-7: left recursion
        for j₋ in reverse((j_mid+1):(j_next-1))
            I[s], J[s] = i, j₋
            V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * _ΔAᵖ_I(Aᵖ⁻¹,K,K′,i,j₋+1) / p′
            s += 1
        end
    end

    Aᵖ = sparse(view(I,1:s-1), view(J,1:s-1), view(V,1:s-1), n, n′)
    return Aᵖ
end

## Uniform B-spline space
function _changebasis_R(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    j = findfirst(==(k[1]), _vec(k′))
    A = spzeros(U, n, n′)
    for i in 1:n
        A[i, j+r*(i-1):j+(r-1)*(p+1)+r*(i-1)] = block
    end
    return A
end
function _changebasis_I(P::UniformBSplineSpace{0,T}, P′::UniformBSplineSpace{0,T′}) where {T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    n = dim(P)
    n′ = dim(P′)
    A = spzeros(U, n, n′)
    for i in 1:n
        b = r*i
        a = b-r+1
        rr = a:b
        if rr ⊆ 1:n′
            A[i,rr] .= one(U)
        else
            A[i,max(a,1):min(b,n′)] .= one(U)
        end
    end
    return A
end
function _changebasis_I(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    A = spzeros(U, n, n′)
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
function _changebasis_R(dP::BSplineDerivativeSpace{r,<:BSplineSpace{p, T}}, P′::BSplineSpace{p′,T′}) where {r,p,p′,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(dP)
    n = dim(dP)
    A = SparseMatrixCSC{U,Int32}(I,n,n)
    for _r in 0:r-1
        _p = p - _r
        _A = spzeros(U, n+_r, n+1+_r)
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
    P′ = bsplinespace(dP′)
    return _changebasis_R(dP, P′)
end
function _changebasis_R(dP::BSplineDerivativeSpace{r}, dP′::BSplineDerivativeSpace{r′}) where {r,r′}
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
function _changebasis_R(P::BSplineSpace, dP′::BSplineDerivativeSpace{0})
    P′ = bsplinespace(dP′)
    return _changebasis_R(P, P′)
end

function changebasis(P::AbstractFunctionSpace, P′::AbstractFunctionSpace)
    P ⊆ P′ && return _changebasis_R(P, P′)
    P ⊑ P′ && return _changebasis_I(P, P′)
    throw(DomainError((P, P′),"P ⊆ P′ or P ⊑ P′ must hold."))
end
