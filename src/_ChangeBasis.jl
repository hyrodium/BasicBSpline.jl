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
    n′_exact = exactdim(P′)
    k = knotvector(P)
    k′ = knotvector(P′)
    I = Vector{Int32}(undef, n′_exact)
    J = Vector{Int32}(undef, n′_exact)
    s = 1
    j_begin = 1
    for i in 1:n
        isdegenerate(P, i) && continue
        for j in j_begin:n′
            k′[j] == k[i] && (j_begin = j; break)
        end
        local j_end
        for j in j_begin:n′
            k′[j+p′+1] == k[i+1] && (j_end = j + p′; break)
        end
        for j in j_begin:j_end
            isdegenerate(P′, j) && continue
            I[s] = i
            J[s] = j
            s += 1
        end
        j_begin = j_end + 1
    end
    A⁰ = sparse(view(I,1:s-1), view(J,1:s-1), fill(one(U), s-1), n, n′)
    return A⁰
end

function _changebasis_R(P::BSplineSpace{p,T,KnotVector{T}}, P′::BSplineSpace{p′,T′,KnotVector{T′}}) where {p,p′,T,T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    k = knotvector(P)
    k′ = knotvector(P′)
    n = dim(P)
    n′ = dim(P′)
    K′ = [k′[i+p′] - k′[i] for i in 1:n′+1]
    K = U[ifelse(k[i+p] ≠ k[i], U(1 / (k[i+p] - k[i])), zero(U)) for i in 1:n+1]
    Aᵖ⁻¹ = _changebasis_R(_lower(P), _lower(P′))  # (n+1) × (n′+1) matrix
    n_nonzero = exactdim(P′)*(p+1)  # This is a upper bound of the number of non-zero elements of Aᵖ (rough estimation).
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
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #     |----*----------------*------>|
        #          j_prev  j_mid    j_next
        #          |       |        |
        #           |--j₊->||<-j₋--|
        #          266666666777777773
        #          366666666777777773
        #
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #    *|---------------*------------>|
        #    j_prev  j_mid    j_next
        #    |       |        |
        #     |--j₊->||<-j₋--|
        #    066666666777777773
        #
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #     |-------------*-------------->|*
        #                   j_prev  j_mid    j_next
        #                   |       |        |
        #                    |--j₊->||<-j₋--|
        #                   266666666777777770
        #                   366666666777777770
        #
        #  1                                    n′
        #  *-----------------------------*----->|
        #  j_begin                       j_end
        #  *----------------*----------->|
        #  j_prev  j_mid    j_next
        #  |       |        |
        #   |--j₊->||<-j₋--|
        #  466666666777777773
        #
        #  1                                    n′
        #  |------*---------------------------->*
        #         j_begin                       j_end
        #         |------------*--------------->*
        #                      j_prev  j_mid    j_next
        #                      |       |        |
        #                       |--j₊->||<-j₋--|
        #                      266666666777777775
        #                      366666666777777775

        # Precalculate the range of j
        Pi = BSplineSpace{p}(view(k, i:i+p+1))
        j_end::Int = findnext(j->Pi ⊆ BSplineSpace{p′}(view(k′, j_begin:j+p′+1)), 1:n′, j_end)
        j_begin::Int = findprev(j->Pi ⊆ BSplineSpace{p′}(view(k′, j:j_end+p′+1)), 1:n′, j_end)
        j_range = j_begin:j_end
        j_prev = j_begin-1
        # flag = 0
        Aᵖᵢⱼ_prev = zero(U)

        for j_next in j_range
            # Rule-1: zero
            if k′[j_next] == k′[j_next+p′+1]
                # flag = 1
                continue
            # Rule-2: right limit
            elseif k′[j_next] == k′[j_next+p′]
                I[s] = i
                J[s] = j_next
                V[s] = Aᵖᵢⱼ_prev = bsplinebasis₊₀(P,i,k′[j_next+1])
                s += 1
                j_prev = j_next
                # flag = 2
            # Rule-3: left limit (or both limit)
            elseif k′[j_next+1] == k′[j_next+p′]
                j_mid = (j_prev + j_next) ÷ 2
                # Rule-6: right recursion
                for j₊ in (j_prev+1):j_mid
                    I[s] = i
                    J[s] = j₊
                    V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * K′[j₊] * (K[i] * Aᵖ⁻¹[i, j₊] - K[i+1] * Aᵖ⁻¹[i+1, j₊]) / p′
                    s += 1
                end
                I[s] = i
                J[s] = j_next
                V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_next = bsplinebasis₋₀(P,i,k′[j_next+1])
                s += 1
                # Rule-7: left recursion
                for j₋ in reverse((j_mid+1):(j_next-1))
                    I[s] = i
                    J[s] = j₋
                    V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * K′[j₋+1] * (K[i] * Aᵖ⁻¹[i, j₋+1] - K[i+1] * Aᵖ⁻¹[i+1, j₋+1]) / p′
                    s += 1
                end
                j_prev = j_next
                # flag = 3
            # Rule-4: left boundary
            elseif j_next == 1
                j_prev = j_next
                I[s] = i
                J[s] = j_next
                V[s] = Aᵖᵢⱼ_prev = p * K′[j_next] * (K[i] * Aᵖ⁻¹[i, j_next] - K[i+1] * Aᵖ⁻¹[i+1, j_next]) / p′
                s += 1
                # flag = 4
            # Rule-5: right boundary
            elseif j_next == n′
                j_mid = (j_prev + j_next) ÷ 2
                # Rule-6: right recursion
                for j₊ in (j_prev+1):j_mid
                    I[s] = i
                    J[s] = j₊
                    V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * K′[j₊] * (K[i] * Aᵖ⁻¹[i, j₊] - K[i+1] * Aᵖ⁻¹[i+1, j₊]) / p′
                    s += 1
                end
                I[s] = i
                J[s] = j_next
                V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_next = -p * K′[j_next+1] * (K[i] * Aᵖ⁻¹[i, j_next+1] - K[i+1] * Aᵖ⁻¹[i+1, j_next+1]) / p′
                s += 1
                # Rule-7: left recursion
                for j₋ in reverse((j_mid+1):(j_next-1))
                    I[s] = i
                    J[s] = j₋
                    V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * K′[j₋+1] * (K[i] * Aᵖ⁻¹[i, j₋+1] - K[i+1] * Aᵖ⁻¹[i+1, j₋+1]) / p′
                    s += 1
                end
                j_prev = j_next
                # flag = 5
            end
        end
        j_next = j_end + 1
        j_mid = (j_prev + j_next) ÷ 2
        # Rule-6: right recursion
        for j₊ in (j_prev+1):j_mid
            I[s] = i
            J[s] = j₊
            V[s] = Aᵖᵢⱼ_prev = Aᵖᵢⱼ_prev + p * K′[j₊] * (K[i] * Aᵖ⁻¹[i, j₊] - K[i+1] * Aᵖ⁻¹[i+1, j₊]) / p′
            s += 1
        end
        Aᵖᵢⱼ_next = zero(U)
        # Rule-7: left recursion
        for j₋ in reverse((j_mid+1):(j_next-1))
            I[s] = i
            J[s] = j₋
            V[s] = Aᵖᵢⱼ_next = Aᵖᵢⱼ_next - p * K′[j₋+1] * (K[i] * Aᵖ⁻¹[i, j₋+1] - K[i+1] * Aᵖ⁻¹[i+1, j₋+1]) / p′
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

    A = sparse(U, I, n, n)

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

function _changebasis_I(P::BSplineSpace{p,T}, P′::BSplineSpace{p′,T′}) where {p,p′,T,T′}
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

## Uniform B-spline space
function _changebasis_R(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    j = findfirst(==(k[1]), _vec(k′))
    A = spzeros(StaticArrays.arithmetic_closure(T), n, n′)
    for i in 1:n
        A[i, j+r*(i-1):j+(r-1)*(p+1)+r*(i-1)] = block
    end
    return A
end
function _changebasis_I(P::UniformBSplineSpace{p,T}, P′::UniformBSplineSpace{p,T′}) where {p,T,T′}
    k = knotvector(P)
    k′ = knotvector(P′)
    r = round(Int, step(k)/step(k′))
    block = [r_nomial(p+1,i,r) for i in 0:(r-1)*(p+1)]/r^p
    n = dim(P)
    n′ = dim(P′)
    A = spzeros(StaticArrays.arithmetic_closure(T), n, n′)
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
    A = sparse(U, I, n, n)
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
