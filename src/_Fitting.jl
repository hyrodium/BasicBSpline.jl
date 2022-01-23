# Fitting

@doc raw"""
Calculate a matrix
```math
A_{ij}=\int_{I} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt
```
"""
function innerproduct_I(P::BSplineSpace{p}) where p
    n = dim(P)
    A = zeros(n,n)
    k = knotvector(P)
    l = length(k)
    @inbounds nodes, weights = SVector{p+1}.(gausslegendre(p+1))
    for m in 1:l-2p-1
        @inbounds t1 = k[m+p]
        @inbounds t2 = k[m+p+1]
        width = t2-t1
        iszero(width) && continue
        dnodes = (width * nodes .+ (t1+t2)) / 2
        bbs = hcat(bsplinebasisall.(P,m,dnodes)...)
        for i in 1:n
            for q in 0:p
                j = i + q
                n < j && continue
                j ≤ m+p || continue
                m+p+1 ≤ i+p+1 || continue
                @inbounds A[i,j] += sum(bbs[i-m+1,:].*bbs[j-m+1,:].*weights)*width/2
            end
        end
    end
    return Symmetric(A)
end

@doc raw"""
Calculate a matrix
```math
A_{ij}=\int_{\mathbb{R}} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt
```
"""
function innerproduct_R(P::BSplineSpace{p}) where p
    n = dim(P)
    A = innerproduct_I(P).data
    k = knotvector(P)
    nodes, weights = SVector{p+1}.(gausslegendre(p+1))
    A1 = zeros(p,p)
    A2 = zeros(p,p)
    for i in 1:p, j in 1:p
        j < i && continue
        for m in 1:p-j+1
            t1 = k[j+m-1]
            t2 = k[j+m]
            width = t2-t1
            iszero(width) && continue
            dnodes = (width * nodes .+ (t1+t2)) / 2
            F = bsplinebasis.(P, i, dnodes) .* bsplinebasis.(P, j, dnodes)
            A1[i,j] += sum(F .* weights)*width/2
        end
        for m in p-j+2:p-j+i+1
            t1 = k[j-p+n+m-1]
            t2 = k[j-p+n+m]
            width = t2-t1
            iszero(width) && continue
            dnodes = (width * nodes .+ (t1+t2)) / 2
            F = bsplinebasis.(P, i-p+n, dnodes) .* bsplinebasis.(P, j-p+n, dnodes)
            A2[i,j] += sum(F .* weights)*width/2
        end
    end
    A[1:p,1:p] += A1
    A[n-p+1:n,n-p+1:n] += A2
    return Symmetric(A)
end

function _f_b_int_R(func, i1, P1::BSplineSpace{p1}, gl1) where {p1}
    k1 = knotvector(P1)
    F(t1) = bsplinebasis(P1, i1, t1) * func(t1)

    f1,l1 = i1, i1+p1

    S1 = integrate(F, k1[f1]..k1[f1+1], gl1)
    for j1 in f1+1:l1
        S1 += integrate(F, k1[j1]..k1[j1+1], gl1)
    end
    return S1
end

function _f_b_int_R(func, i1, i2, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}}, gl1, gl2) where {p1,p2}
    P1, P2 = P
    k1, k2 = knotvector(P1), knotvector(P2)
    F(t1, t2) = bsplinebasis(P1, i1, t1) * bsplinebasis(P2, i2, t2) * func(t1, t2)

    f1, l1 = i1, i1+p1
    f2, l2 = i2, i2+p2

    S2 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], gl1, gl2)
    for j2 in f2+1:l2
        S2 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], gl1, gl2)
    end
    S1 = S2
    for j1 in f1+1:l1
        S2 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], gl1, gl2)
        for j2 in f2+1:l2
            S2 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], gl1, gl2)
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_I(func, i1, i2, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}}, gl1, gl2) where {p1,p2}
    P1, P2 = P
    k1, k2 = knotvector(P1), knotvector(P2)
    m1, m2 = length(k1), length(k2)
    F(t1, t2) = bsplinebasis(P1, i1, t1) * bsplinebasis(P2, i2, t2) * func(t1, t2)

    f1, l1 = max(i1, p1+1), min(i1+p1, m1-p1-1)
    f2, l2 = max(i2, p2+1), min(i2+p2, m2-p2-1)

    S2 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], gl1, gl2)
    for j2 in f2+1:l2
        S2 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], gl1, gl2)
    end
    S1 = S2
    for j1 in f1+1:l1
        S2 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], gl1, gl2)
        for j2 in f2+1:l2
            S2 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], gl1, gl2)
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_R(func, i1, i2, i3, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}, <:AbstractBSplineSpace{p3}}, gl1, gl2, gl3) where {p1,p2,p3}
    P1, P2, P3 = P
    k1, k2, k3 = knotvector(P1), knotvector(P2), knotvector(P3)
    F(t1, t2, t3) = bsplinebasis(P1, i1, t1) * bsplinebasis(P2, i2, t2) * bsplinebasis(P3, i3, t3) * func(t1, t2, t3)

    f1, l1 = i1, i1+p1
    f2, l2 = i2, i2+p2
    f3, l3 = i3, i3+p3

    S3 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
    for j3 in f3+1:l3
        S3 += integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
    end
    S2 = S3
    for j2 in f2+1:l2
        S3 = integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
        end
        S2 += S3
    end
    S1 = S2
    for j1 in f1+1:l1
        S3 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
        end
        S2 = S3
        for j2 in f2+1:l2
            S3 = integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
            for j3 in f3+1:l3
                S3 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
            end
            S2 += S3
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_I(func, i1, i2, i3, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}, <:AbstractBSplineSpace{p3}}, gl1, gl2, gl3) where {p1,p2,p3}
    P1, P2, P3 = P
    k1, k2, k3 = knotvector(P1), knotvector(P2), knotvector(P3)
    m1, m2, m3 = length(k1), length(k2), length(k3)
    F(t1, t2, t3) = bsplinebasis(P1, i1, t1) * bsplinebasis(P2, i2, t2) * bsplinebasis(P3, i3, t3) * func(t1, t2, t3)

    f1, l1 = max(i1, p1+1), min(i1+p1, m1-p1-1)
    f2, l2 = max(i2, p2+1), min(i2+p2, m2-p2-1)
    f3, l3 = max(i3, p3+1), min(i3+p3, m3-p3-1)

    S3 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
    for j3 in f3+1:l3
        S3 += integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
    end
    S2 = S3
    for j2 in f2+1:l2
        S3 = integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
        end
        S2 += S3
    end
    S1 = S2
    for j1 in f1+1:l1
        S3 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
        end
        S2 = S3
        for j2 in f2+1:l2
            S3 = integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], gl1, gl2, gl3)
            for j3 in f3+1:l3
                S3 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], gl1, gl2, gl3)
            end
            S2 += S3
        end
        S1 += S2
    end
    return S1
end

function innerproduct_R(func, Ps::Tuple{<:BSplineSpace{p1}}) where {p1}
    P1, = Ps
    n1 = dim(P1)
    nip1 = p1 + 1
    gl1 = GaussLegendre(nip1)
    b = [_f_b_int_R(func, i1, P1, gl1) for i1 in 1:n1]
    return b
end

function innerproduct_I(func, Ps::Tuple{<:BSplineSpace{p₁}}) where {p₁}
    P₁, = Ps
    n₁ = dim(P₁)
    k₁ = knotvector(P₁)
    l₁ = length(k₁)
    nodes₁, weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))

    sample_point = func(leftendpoint(domain(P₁)))
    b = Array{typeof(sample_point),1}(undef, n₁)
    fill!(b, zero(sample_point))

    for m₁ in 1:l₁-2p₁-1
        ta₁ = k₁[m₁+p₁]
        tb₁ = k₁[m₁+p₁+1]
        w₁ = tb₁-ta₁
        iszero(w₁) && continue
        dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
        bbs₁ = hcat(bsplinebasisall.(P₁,m₁,dnodes₁)...)
        F = func.(dnodes₁)
        for q₁ in 0:p₁
            i₁ = m₁ + q₁
            b[i₁] += sum(bbs₁[i₁-m₁+1,:].*F.*weights₁)*w₁/2
        end
    end
    return b
end

function innerproduct_R(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}}) where {p1,p2}
    n1, n2 = dim.(P)
    nip1 = p1 + 1
    nip2 = p2 + 1
    gl1 = GaussLegendre(nip1)
    gl2 = GaussLegendre(nip2)
    b = [_f_b_int_R(func, i1, i2, P, gl1, gl2) for i1 in 1:n1, i2 in 1:n2]
    return b
end

function innerproduct_I(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}}) where {p1,p2}
    n1, n2 = dim.(P)
    nip1 = p1 + 1
    nip2 = p2 + 1
    gl1 = GaussLegendre(nip1)
    gl2 = GaussLegendre(nip2)
    b = [_f_b_int_I(func, i1, i2, P, gl1, gl2) for i1 in 1:n1, i2 in 1:n2]
    return b
end

function innerproduct_R(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}, <:AbstractBSplineSpace{p3}}) where {p1,p2,p3}
    n1, n2, n3 = dim.(P)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nip3 = p3 + 1
    gl1 = GaussLegendre(nip1)
    gl2 = GaussLegendre(nip2)
    gl3 = GaussLegendre(nip3)
    b = [_f_b_int_R(func, i1, i2, i3, P, gl1, gl2, gl3) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
    return b
end

function innerproduct_I(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}, <:AbstractBSplineSpace{p3}}) where {p1,p2,p3}
    P1, P2, P3 = P
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nip3 = p3 + 1
    gl1 = GaussLegendre(nip1)
    gl2 = GaussLegendre(nip2)
    gl3 = GaussLegendre(nip3)
    b = [_f_b_int_I(func, i1, i2, i3, P, gl1, gl2, gl3) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
    return b
end

"""
* func: Real -> ℝ-vector space
"""
function fittingcontrolpoints(func, P::Tuple{<:AbstractBSplineSpace{p1}}; domain=:I) where {p1}
    P1, = P
    if domain == :I
        b = innerproduct_I(func, P)
        A = innerproduct_I(P1)
    elseif domain == :R
        b = innerproduct_R(func, P)
        A = innerproduct_R(P1)
    end
    return inv(A) * b
end

"""
* func: (Real,Real) -> ℝ-vector space
"""
function fittingcontrolpoints(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}}; domain=:I) where {p1,p2}
    P1, P2 = P
    n1, n2 = dim.(P)
    if domain == :I
        A1, A2 = innerproduct_I(P1), innerproduct_I(P2)
        b = innerproduct_I(func, P)
    elseif domain == :R
        A1, A2 = innerproduct_R(P1), innerproduct_R(P2)
        b = innerproduct_R(func, P)
    end
    A = [A1[i1, j1] * A2[i2, j2] for i1 in 1:n1, i2 in 1:n2, j1 in 1:n1, j2 in 1:n2]
    _A = reshape(A, n1 * n2, n1 * n2)
    _b = reshape(b, n1 * n2)
    return reshape(inv(_A) * _b, n1, n2)
end

"""
* func: (Real,Real,Real) -> ℝ-vector space
"""
function fittingcontrolpoints(func, P::Tuple{<:AbstractBSplineSpace{p1}, <:AbstractBSplineSpace{p2}, <:AbstractBSplineSpace{p3}}; domain=:I) where {p1,p2,p3}
    P1, P2, P3 = P
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    if domain == :I
        A1, A2, A3 = innerproduct_I(P1), innerproduct_I(P2), innerproduct_I(P3)
        b = innerproduct_I(func, P)
    elseif domain == :R
        A1, A2, A3 = innerproduct_R(P1), innerproduct_R(P2), innerproduct_R(P3)
        b = innerproduct_R(func, P)
    end
    A = [A1[i1, j1] * A2[i2, j2] * A3[i3, j3] for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3, j1 in 1:n1, j2 in 1:n2, j3 in 1:n3]
    _A = reshape(A, n1 * n2 * n3, n1 * n2 * n3)
    _b = reshape(b, n1 * n2 * n3)
    return reshape(inv(_A) * _b, n1, n2, n3)
end
