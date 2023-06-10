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
    I = fill(0, n′_exact)
    J = fill(0, n′_exact)
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

"""
Generate flags for `changebasis` function.

0 -> default
1 -> zero            k′[i] == k′[i+p′+1]
2 -> left limit      k′[i] < k′[i+1] == k′[i+p′+1]
3 -> right limit     k′[i+1] == k′[i+p′] ( i.e. k′[i] == k′[i+p′] < k′[i+p′+1] or k′[i] < k′[i+1] == k′[i+p′] < k′[i+p′+1] )
4 -> left boundary   i == 1
5 -> right boundary  i == n′
6 -> left recursion
7 -> right recursion
"""
function _generate_flags(P′::BSplineSpace{p′}) where p′
    k′ = knotvector(P′)
    n′ = dim(P′)
    flags = zeros(Int, n′)
    for i in 1:n′
        if k′[i] == k′[i+p′+1]
            flags[i] = 1
        elseif k′[i] < k′[i+1] == k′[i+p′+1]
            flags[i] = 2
        elseif k′[i+1] == k′[i+p′]
            flags[i] = 3
        end
    end
    if flags[1] == 0
        flags[1] = 4
    end
    if flags[end] == 0
        flags[end] = 5
    end
    local j
    for i in 2:n′-1
        if flags[i-1] ≠ 0 && flags[i-1] ≠ 6 && flags[i] == 0
            flags[i] = 6
            j = i
        end
        if flags[i+1] ≠ 0 && flags[i] == 0
            flags[i] = 7
            jj = (i + j) ÷ 2
            ii = jj + 1
            flags[j:jj] .= 6
            flags[ii:i] .= 7
        end
    end
    return flags
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
    Aᵖ = spzeros(U, n, n′)  # n × n′ matrix

    flags = _generate_flags(P′)
    if flags[1] == 4
        Aᵖ[1, 1] = p * K′[1] * (K[1] * Aᵖ⁻¹[1, 1] - K[1+1] * Aᵖ⁻¹[1+1, 1]) / p′
    end
    if flags[n′] == 5
        Aᵖ[n, n′] = -p * K′[n′+1] * (K[n] * Aᵖ⁻¹[n, n′+1] - K[n+1] * Aᵖ⁻¹[n+1, n′+1]) / p′
    end
    for i in 1:n
        # Skip for degenerated basis
        isdegenerate_R(P,i) && continue

        # The following indices `j*` have the following relationships.
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #     |----*----------------*------>|
        #          j_tmp   j_mid    j
        #          |       |        |
        #           |--j₊->||<-j₋--|
        #          366666666777777773
        #          366666666777777772
        #
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #    *|---------------*------------>|
        #    j_tmp   j_mid    j
        #    |       |        |
        #     |--j₊->||<-j₋--|
        #    066666666777777773
        #    066666666777777772
        #
        #  1                                    n′
        #  |--*-----------------------------*-->|
        #     j_begin                       j_end
        #     |-------------*-------------->|*
        #                   j_tmp   j_mid    j
        #                   |       |        |
        #                    |--j₊->||<-j₋--|
        #                   366666666777777770
        #
        #  1                                    n′
        #  *-----------------------------*----->|
        #  j_begin                       j_end
        # *|---------------*------------>|
        # j_tmp   j_mid    j
        # |       |        |
        #  |--j₊->||<-j₋--|
        #  46666666777777773
        #  46666666777777772
        #
        #  1                                    n′
        #  |------*---------------------------->*
        #         j_begin                       j_end
        #         |-------------*-------------->|*
        #                       j_tmp   j_mid    j
        #                       |       |        |
        #                        |--j₊->||<-j₋--|
        #                       36666666677777775

        # Precalculate the range of j
        # TODO: this implementation can be replaced with more effecient way.
        j_begin::Int = findlast(j->BSplineSpace{p}(k[i:i+p+1]) ⊆ BSplineSpace{p′}(k′[j:n′+p′+1]), 1:n′)
        j_end::Int = findnext(j->BSplineSpace{p}(k[i:i+p+1]) ⊆ BSplineSpace{p′}(k′[j_begin:j+p′+1]), 1:n′, j_begin)
        J = j_begin:j_end

        # Elements that can be calculated with left or right limit of B-spline basis functions
        for j in J
            if flags[j] == 2
                Aᵖ[i, j] = bsplinebasis₋₀(P,i,k′[j+1])
            elseif flags[j] == 3
                Aᵖ[i, j] = bsplinebasis₊₀(P,i,k′[j+1])
            end
        end
        for j in J
            if flags[j] == 6
                Aᵖ[i, j] = Aᵖ[i, j-1] + p * K′[j] * (K[i] * Aᵖ⁻¹[i, j] - K[i+1] * Aᵖ⁻¹[i+1, j]) / p′
            end
        end
        for j in reverse(J)
            if flags[j] == 7
                Aᵖ[i, j] = Aᵖ[i, j+1] - p * K′[j+1] * (K[i] * Aᵖ⁻¹[i, j+1] - K[i+1] * Aᵖ⁻¹[i+1, j+1]) / p′
            end
        end
    end

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
