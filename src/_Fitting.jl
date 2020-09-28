# Fitting

@doc raw"""
```math
\int_{[k_{j+n-1}, k_{j+n}]} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt
```
Assumption:
* i ≤ j
* 1 ≤ n ≤ p-j+i+1
"""
function _b_b_int(P::AbstractBSplineSpace, i, j, n, nip, nodes, weights)
    p = degree(P)
    I = knots(P)[j+n-1]..knots(P)[j+n]

    f(t) = bsplinebasis(i, P, t) * bsplinebasis(j, P, t)
    return integrate(f, I, nip, nodes, weights)
end

@doc raw"""
```math
\begin{align}
&\int_{\mathbb{R}} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt \\
={}&
\begin{cases}
\displaystyle \int_{[k_{j}, k_{i+p+1}]} B_{(i,p,k)}(t)  B_{(j,p,k)}(t) dt & (i \le j) \\
\displaystyle \int_{[k_{i}, k_{j+p+1}]} B_{(i,p,k)}(t)  B_{(j,p,k)}(t) dt & (j \le i)
\end{cases}
\end{align}
```
"""
function _b_b_int_R(P::AbstractBSplineSpace, i, j, nip, nodes, weights)
    p = degree(P)
    k = knots(P)
    Δ = j - i
    if Δ < -p
        return 0.0
    elseif Δ > p
        return 0.0
    elseif Δ ≥ 0  # i ≤ j
        s = 0.0
        for n in 1:p-j+i+1
            s += _b_b_int(P, i, j, n, nip, nodes, weights)
        end
        return s
    else  # j < i
        s = 0.0
        for n in 1:p-i+j+1
            s += _b_b_int(P, j, i, n, nip, nodes, weights)
        end
        return s
    end
end

@doc raw"""
```math
\begin{align}
&\int_{I} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt \\
\end{align}
```
"""
function _b_b_int_I(P::AbstractBSplineSpace, i, j, nip, nodes, weights)
    p = degree(P)
    k = knots(P)
    n = dim(P)
    Δ = j - i
    if Δ < -p
        return 0.0
    elseif Δ > p
        return 0.0
    elseif Δ ≥ 0  # i ≤ j
        s = 0.0
        for m in max(p-j+2,1):min(n-j+1,p-j+i+1)
            s += _b_b_int(P, i, j, m, nip, nodes, weights)
        end
        return s
    else  # j < i
        s = 0.0
        for m in max(p-i+2,1):min(n-i+1,p-i+j+1)
            s += _b_b_int(P, j, i, m, nip, nodes, weights)
        end
        return s
    end
end

function innerproduct_R(P::AbstractBSplineSpace)
    p = degree(P)
    n = dim(P)
    nip = p + 1
    nodes, weights = gausslegendre(nip)
    return [_b_b_int_R(P, i, j, nip, nodes, weights) for i in 1:n, j in 1:n]
end

function innerproduct_I(P::AbstractBSplineSpace)
    p = degree(P)
    n = dim(P)
    nip = p + 1
    nodes, weights = gausslegendre(nip)
    return [_b_b_int_I(P, i, j, nip, nodes, weights) for i in 1:n, j in 1:n]
end

function _f_b_int_R(func, i1, P1::AbstractBSplineSpace, nip1, nodes1, weights1)
    k1 = knots(P1)
    p1 = degree(P1)
    F(t1) = bsplinebasis(i1, P1, t1) * func([t1])

    f1,l1 = i1, i1+p1

    S1 = integrate(F, k1[f1]..k1[f1+1], nip1, nodes1, weights1)
    for j1 in f1+1:l1
        S1 += integrate(F, k1[j1]..k1[j1+1], nip1, nodes1, weights1)
    end
    return S1
end

function _f_b_int_I(func, i1, P1::AbstractBSplineSpace, nip1, nodes1, weights1)
    k1 = knots(P1)
    m1 = length(k1)
    p1 = degree(P1)
    F(t1) = bsplinebasis(i1, P1, t1) * func([t1])

    f1,l1 = max(i1, p1+1), min(i1+p1, m1-p1-1)

    S1 = integrate(F, k1[f1]..k1[f1+1], nip1, nodes1, weights1)
    for j1 in f1+1:l1
        S1 += integrate(F, k1[j1]..k1[j1+1], nip1, nodes1, weights1)
    end
    return S1
end

function _f_b_int_R(func, i1, i2, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, nip1, nip2, nodes1, nodes2, weights1, weights2)
    k1, k2 = knots(P1), knots(P2)
    p1, p2 = degree(P1), degree(P2)
    F(t1, t2) = bsplinebasis(i1, P1, t1) * bsplinebasis(i2, P2, t2) * func([t1, t2])

    f1, l1 = i1, i1+p1
    f2, l2 = i2, i2+p2

    S2 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
    for j2 in f2+1:l2
        S2 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
    end
    S1 = S2
    for j1 in f1+1:l1
        S2 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
        for j2 in f2+1:l2
            S2 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_I(func, i1, i2, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, nip1, nip2, nodes1, nodes2, weights1, weights2)
    k1, k2 = knots(P1), knots(P2)
    m1, m2 = length(k1), length(k2)
    p1, p2 = degree(P1), degree(P2)
    F(t1, t2) = bsplinebasis(i1, P1, t1) * bsplinebasis(i2, P2, t2) * func([t1, t2])

    f1, l1 = max(i1, p1+1), min(i1+p1, m1-p1-1)
    f2, l2 = max(i2, p2+1), min(i2+p2, m2-p2-1)

    S2 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
    for j2 in f2+1:l2
        S2 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
    end
    S1 = S2
    for j1 in f1+1:l1
        S2 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
        for j2 in f2+1:l2
            S2 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], nip1, nip2, nodes1, nodes2, weights1, weights2)
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_R(func, i1, i2, i3, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace, nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    k1, k2, k3 = knots(P1), knots(P2), knots(P3)
    p1, p2, p3 = degree(P1), degree(P2), degree(P3)
    F(t1, t2, t3) = bsplinebasis(i1, P1, t1) * bsplinebasis(i2, P2, t2) * bsplinebasis(i3, P3, t3) * func([t1, t2, t3])

    f1, l1 = i1, i1+p1
    f2, l2 = i2, i2+p2
    f3, l3 = i3, i3+p3

    S3 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    for j3 in f3+1:l3
        S3 += integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    end
    S2 = S3
    for j2 in f2+1:l2
        S3 = integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        end
        S2 += S3
    end
    S1 = S2
    for j1 in f1+1:l1
        S3 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        end
        S2 = S3
        for j2 in f2+1:l2
            S3 = integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
            for j3 in f3+1:l3
                S3 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
            end
            S2 += S3
        end
        S1 += S2
    end
    return S1
end

function _f_b_int_I(func, i1, i2, i3, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace, nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    k1, k2, k3 = knots(P1), knots(P2), knots(P3)
    m1, m2, m3 = length(k1), length(k2), length(k3)
    p1, p2, p3 = degree(P1), degree(P2), degree(P3)
    F(t1, t2, t3) = bsplinebasis(i1, P1, t1) * bsplinebasis(i2, P2, t2) * bsplinebasis(i3, P3, t3) * func([t1, t2, t3])

    f1, l1 = max(i1, p1+1), min(i1+p1, m1-p1-1)
    f2, l2 = max(i2, p2+1), min(i2+p2, m2-p2-1)
    f3, l3 = max(i3, p3+1), min(i3+p3, m3-p3-1)

    S3 = integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    for j3 in f3+1:l3
        S3 += integrate(F, k1[f1]..k1[f1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    end
    S2 = S3
    for j2 in f2+1:l2
        S3 = integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[f1]..k1[f1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        end
        S2 += S3
    end
    S1 = S2
    for j1 in f1+1:l1
        S3 = integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        for j3 in f3+1:l3
            S3 += integrate(F, k1[j1]..k1[j1+1], k2[f2]..k2[f2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
        end
        S2 = S3
        for j2 in f2+1:l2
            S3 = integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[f3]..k3[f3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
            for j3 in f3+1:l3
                S3 += integrate(F, k1[j1]..k1[j1+1], k2[j2]..k2[j2+1], k3[j3]..k3[j3+1], nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
            end
            S2 += S3
        end
        S1 += S2
    end
    return S1
end

function innerproduct_R(func, P1::AbstractBSplineSpace)
    p1 = degree(P1)
    n1 = dim(P1)
    nip1 = p1 + 1
    nodes1, weights1 = gausslegendre(nip1)
    b = [_f_b_int_R(func, i1, P1, nip1, nodes1, weights1) for i1 in 1:n1]
end

function innerproduct_I(func, P1::AbstractBSplineSpace)
    p1 = degree(P1)
    n1 = dim(P1)
    nip1 = p1 + 1
    nodes1, weights1 = gausslegendre(nip1)
    b = [_f_b_int_I(func, i1, P1, nip1, nodes1, weights1) for i1 in 1:n1]
end

function innerproduct_R(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace)
    p1, p2 = degree(P1), degree(P2)
    n1, n2 = dim(P1), dim(P2)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nodes1, weights1 = gausslegendre(nip1)
    nodes2, weights2 = gausslegendre(nip2)
    b = [_f_b_int_R(func, i1, i2, P1, P2, nip1, nip2, nodes1, nodes2, weights1, weights2) for i1 in 1:n1, i2 in 1:n2]
end

function innerproduct_I(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace)
    p1, p2 = degree(P1), degree(P2)
    n1, n2 = dim(P1), dim(P2)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nodes1, weights1 = gausslegendre(nip1)
    nodes2, weights2 = gausslegendre(nip2)
    b = [_f_b_int_I(func, i1, i2, P1, P2, nip1, nip2, nodes1, nodes2, weights1, weights2) for i1 in 1:n1, i2 in 1:n2]
end

function innerproduct_R(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace)
    p1, p2, p3 = degree(P1), degree(P2), degree(P3)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nip3 = p3 + 1
    nodes1, weights1 = gausslegendre(nip1)
    nodes2, weights2 = gausslegendre(nip2)
    nodes3, weights3 = gausslegendre(nip3)
    b = [_f_b_int_R(func, i1, i2, i3, P1, P2, P3, nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
end

function innerproduct_I(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace)
    p1, p2, p3 = degree(P1), degree(P2), degree(P3)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    nip1 = p1 + 1
    nip2 = p2 + 1
    nip3 = p3 + 1
    nodes1, weights1 = gausslegendre(nip1)
    nodes2, weights2 = gausslegendre(nip2)
    nodes3, weights3 = gausslegendre(nip3)
    b = [_f_b_int_I(func, i1, i2, i3, P1, P2, P3, nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3) for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3]
end

"""
* func: Array{Real,1} -> ℝ-linear space
"""
function fittingcontrolpoints_1dim(func, P1::AbstractBSplineSpace)
    b = innerproduct_I(func, P1)
    A = innerproduct_I(P1)
    return inv(A) * b
end

function fittingcontrolpoints_2dim(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace)
    n1, n2 = dim(P1), dim(P2)
    A1, A2 = innerproduct_I(P1), innerproduct_I(P2)
    A = [A1[i1, j1] * A2[i2, j2] for i1 in 1:n1, i2 in 1:n2, j1 in 1:n1, j2 in 1:n2]
    b = innerproduct_I(func, P1, P2)
    _A = reshape(A, n1 * n2, n1 * n2)
    _b = reshape(b, n1 * n2)
    return reshape(inv(_A) * _b, n1, n2)
end

function fittingcontrolpoints_3dim(func, P1::AbstractBSplineSpace, P2::AbstractBSplineSpace, P3::AbstractBSplineSpace)
    n1, n2, n3 = dim(P1), dim(P2), dim(P3)
    A1, A2, A3 = innerproduct_I(P1), innerproduct_I(P2), innerproduct_I(P3)
    A = [A1[i1, j1] * A2[i2, j2] * A3[i3, j3] for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3, j1 in 1:n1, j2 in 1:n2, j3 in 1:n3]
    b = innerproduct_I(func, P1, P2, P3)
    _A = reshape(A, n1 * n2 * n3, n1 * n2 * n3)
    _b = reshape(b, n1 * n2 * n3)
    return reshape(inv(_A) * _b, n1, n2, n3)
end


"""
Approximate given function by linear combination of B-spline functions.
This function returns its control points.
"""
function fittingcontrolpoints(func, Ps::Array{<:AbstractBSplineSpace,1})
    # TODO: currently, this function only supports for 1-dim and 2-dim B-spline manifold.
    d = length(Ps)
    if d == 1
        return fittingcontrolpoints_1dim(func, Ps[1])
    elseif d == 2
        return fittingcontrolpoints_2dim(func, Ps[1], Ps[2])
    elseif d == 3
        return fittingcontrolpoints_3dim(func, Ps[1], Ps[2], Ps[3])
    else
        error("fittingcontrolpoints supports only 1 to 3 dimensions.")
    end
end
