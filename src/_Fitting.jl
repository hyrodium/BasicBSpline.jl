# Fitting

@doc raw"""
Calculate a matrix
```math
A_{ij}=\int_{I} B_{(i,p,k)}(t) B_{(j,p,k)}(t) dt
```
"""
function innerproduct_I(P::AbstractBSplineSpace{p}) where p
    Base.depwarn("`innterproduct_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_I)
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
function innerproduct_R(P::AbstractBSplineSpace{p}) where p
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
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

function innerproduct_I(func, Ps::Tuple{AbstractBSplineSpace{p₁}}) where {p₁}
    Base.depwarn("`innterproduct_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_I)
    P₁, = Ps
    n₁ = dim(P₁)
    k₁ = knotvector(P₁)
    l₁ = length(k₁)
    nodes₁, weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))

    sample_point = func(float(leftendpoint(domain(P₁))))
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

function innerproduct_I(func, Ps::Tuple{AbstractBSplineSpace{p₁},AbstractBSplineSpace{p₂}}) where {p₁,p₂}
    Base.depwarn("`innterproduct_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_I)
    P₁,P₂ = Ps
    n₁,n₂ = dim.(Ps)
    k₁,k₂ = knotvector.(Ps)
    l₁,l₂ = length(k₁), length(k₂)
    nodes₁, weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))
    nodes₂, weights₂ = SVector{p₂+1}.(gausslegendre(p₂+1))

    sample_point = func(float(leftendpoint(domain(P₁))),float(leftendpoint(domain(P₂))))
    b = Array{typeof(sample_point),2}(undef, n₁, n₂)
    fill!(b, zero(sample_point))

    for m₁ in 1:l₁-2p₁-1
        ta₁ = k₁[m₁+p₁]
        tb₁ = k₁[m₁+p₁+1]
        w₁ = tb₁-ta₁
        iszero(w₁) && continue
        dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
        bbs₁ = hcat(bsplinebasisall.(P₁,m₁,dnodes₁)...)
        for m₂ in 1:l₂-2p₂-1
            ta₂ = k₂[m₂+p₂]
            tb₂ = k₂[m₂+p₂+1]
            w₂ = tb₂-ta₂
            iszero(w₂) && continue
            dnodes₂ = (w₂ * nodes₂ .+ (ta₂+tb₂)) / 2
            bbs₂ = hcat(bsplinebasisall.(P₂,m₂,dnodes₂)...)
            F = func.(dnodes₁,dnodes₂')
            for q₁ in 0:p₁
                i₁ = m₁ + q₁
                for q₂ in 0:p₂
                    i₂ = m₂ + q₂
                    b[i₁,i₂] += sum(F.*bbs₁[i₁-m₁+1,:].*weights₁.*(bbs₂[i₂-m₂+1,:].*weights₂)')*w₁*w₂/4
                end
            end
        end
    end
    return b
end

function innerproduct_I(func, Ps::Tuple{AbstractBSplineSpace{p₁},AbstractBSplineSpace{p₂},AbstractBSplineSpace{p₃}}) where {p₁,p₂,p₃}
    Base.depwarn("`innterproduct_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_I)
    P₁,P₂,P₃ = Ps
    n₁,n₂,n₃ = dim.(Ps)
    k₁,k₂,k₃ = knotvector.(Ps)
    l₁,l₂,l₃ = length(k₁), length(k₂), length(k₃)
    nodes₁, weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))
    nodes₂, weights₂ = SVector{p₂+1}.(gausslegendre(p₂+1))
    nodes₃, weights₃ = SVector{p₃+1}.(gausslegendre(p₃+1))

    sample_point = func(float(leftendpoint(domain(P₁))),float(leftendpoint(domain(P₂))),float(leftendpoint(domain(P₃))))
    b = Array{typeof(sample_point),3}(undef, n₁, n₂, n₃)
    fill!(b, zero(sample_point))

    for m₁ in 1:l₁-2p₁-1
        ta₁ = k₁[m₁+p₁]
        tb₁ = k₁[m₁+p₁+1]
        w₁ = tb₁-ta₁
        iszero(w₁) && continue
        dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
        bbs₁ = hcat(bsplinebasisall.(P₁,m₁,dnodes₁)...)
        for m₂ in 1:l₂-2p₂-1
            ta₂ = k₂[m₂+p₂]
            tb₂ = k₂[m₂+p₂+1]
            w₂ = tb₂-ta₂
            iszero(w₂) && continue
            dnodes₂ = (w₂ * nodes₂ .+ (ta₂+tb₂)) / 2
            bbs₂ = hcat(bsplinebasisall.(P₂,m₂,dnodes₂)...)
            for m₃ in 1:l₃-2p₃-1
                ta₃ = k₃[m₃+p₃]
                tb₃ = k₃[m₃+p₃+1]
                w₃ = tb₃-ta₃
                iszero(w₃) && continue
                dnodes₃ = (w₃ * nodes₃ .+ (ta₃+tb₃)) / 2
                bbs₃ = hcat(bsplinebasisall.(P₃,m₃,dnodes₃)...)
                # TODO: This can be potentially faster
                for j in 1:p₃+1
                    F = func.(dnodes₁,dnodes₂',dnodes₃[j])
                    for q₁ in 0:p₁
                        i₁ = m₁ + q₁
                        for q₂ in 0:p₂
                            i₂ = m₂ + q₂
                            for q₃ in 0:p₃
                                i₃ = m₃ + q₃
                                b[i₁,i₂,i₃] += sum(F.*bbs₁[i₁-m₁+1,:].*weights₁.*(bbs₂[i₂-m₂+1,:].*weights₂)')*bbs₃[i₃-m₃+1,j]*weights₃[j]*w₁*w₂*w₃/8
                            end
                        end
                    end
                end
            end
        end
    end
    return b
end

function innerproduct_R(func, Ps::Tuple{AbstractBSplineSpace{p₁}}) where {p₁}
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    P₁, = Ps
    n₁ = dim(P₁)
    k₁ = knotvector(P₁)
    l₁ = length(k₁)
    nodes₁,weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))
    b = innerproduct_I(func,Ps)
    for i₁ in 1:n₁
        F(t₁) = bsplinebasis(P₁, i₁, t₁) * func(t₁)
        for j₁ in 1:l₁-1
            i₁ ≤ j₁ ≤ i₁+p₁ || continue
            1+p₁ ≤ j₁ ≤ l₁-p₁-1 && continue
            ta₁ = k₁[j₁]
            tb₁ = k₁[j₁+1]
            w₁ = tb₁-ta₁
            iszero(w₁) && continue
            dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
            b[i₁] += sum(F.(dnodes₁).*weights₁)*w₁/2
        end
    end
    return b
end

function innerproduct_R(func, Ps::Tuple{AbstractBSplineSpace{p₁},AbstractBSplineSpace{p₂}}) where {p₁,p₂}
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    P₁,P₂ = Ps
    n₁,n₂ = dim.(Ps)
    k₁,k₂ = knotvector.(Ps)
    l₁,l₂ = length(k₁), length(k₂)
    nodes₁,weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))
    nodes₂,weights₂ = SVector{p₂+1}.(gausslegendre(p₂+1))
    b = innerproduct_I(func,Ps)
    for i₁ in 1:n₁, i₂ in 1:n₂
        F(t₁,t₂) = bsplinebasis(P₁,i₁,t₁) * bsplinebasis(P₂,i₂,t₂) * func(t₁,t₂)
        for j₁ in 1:l₁-1
            i₁ ≤ j₁ ≤ i₁+p₁ || continue
            ta₁ = k₁[j₁]
            tb₁ = k₁[j₁+1]
            w₁ = tb₁-ta₁
            iszero(w₁) && continue
            dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
            for j₂ in 1:l₂-1
                i₂ ≤ j₂ ≤ i₂+p₂ || continue
                (1+p₁ ≤ j₁ ≤ l₁-p₁-1 && 1+p₂ ≤ j₂ ≤ l₂-p₂-1) && continue
                ta₂ = k₂[j₂]
                tb₂ = k₂[j₂+1]
                w₂ = tb₂-ta₂
                iszero(w₂) && continue
                dnodes₂ = (w₂ * nodes₂ .+ (ta₂+tb₂)) / 2
                b[i₁,i₂] += sum(F.(dnodes₁,dnodes₂').*weights₁.*weights₂')*w₁*w₂/4
            end
        end
    end
    return b
end

function innerproduct_R(func, Ps::Tuple{AbstractBSplineSpace{p₁},AbstractBSplineSpace{p₂},AbstractBSplineSpace{p₃}}) where {p₁,p₂,p₃}
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    P₁,P₂,P₃ = Ps
    n₁,n₂,n₃ = dim.(Ps)
    k₁,k₂,k₃ = knotvector.(Ps)
    l₁,l₂,l₃ = length(k₁), length(k₂), length(k₃)
    nodes₁,weights₁ = SVector{p₁+1}.(gausslegendre(p₁+1))
    nodes₂,weights₂ = SVector{p₂+1}.(gausslegendre(p₂+1))
    nodes₃,weights₃ = SVector{p₃+1}.(gausslegendre(p₃+1))
    b = innerproduct_I(func,Ps)
    for i₁ in 1:n₁, i₂ in 1:n₂, i₃ in 1:n₃
        F(t₁,t₂,t₃) = bsplinebasis(P₁,i₁,t₁) * bsplinebasis(P₂,i₂,t₂)  * bsplinebasis(P₃,i₃,t₃) * func(t₁,t₂,t₃)
        for j₁ in 1:l₁-1
            i₁ ≤ j₁ ≤ i₁+p₁ || continue
            ta₁ = k₁[j₁]
            tb₁ = k₁[j₁+1]
            w₁ = tb₁-ta₁
            iszero(w₁) && continue
            dnodes₁ = (w₁ * nodes₁ .+ (ta₁+tb₁)) / 2
            for j₂ in 1:l₂-1
                i₂ ≤ j₂ ≤ i₂+p₂ || continue
                ta₂ = k₂[j₂]
                tb₂ = k₂[j₂+1]
                w₂ = tb₂-ta₂
                iszero(w₂) && continue
                dnodes₂ = (w₂ * nodes₂ .+ (ta₂+tb₂)) / 2
                for j₃ in 1:l₃-1
                    i₃ ≤ j₃ ≤ i₃+p₃ || continue
                    # TODO: This can be potentially faster
                    (1+p₁ ≤ j₁ ≤ l₁-p₁-1 && 1+p₂ ≤ j₂ ≤ l₂-p₂-1 && 1+p₃ ≤ j₃ ≤ l₃-p₃-1) && continue
                    ta₃ = k₃[j₃]
                    tb₃ = k₃[j₃+1]
                    w₃ = tb₃-ta₃
                    iszero(w₃) && continue
                    dnodes₃ = (w₃ * nodes₃ .+ (ta₃+tb₃)) / 2
                    for j in 1:p₃+1
                        b[i₁,i₂,i₃] += weights₃[j]*sum(F.(dnodes₁,dnodes₂',dnodes₃[j]).*weights₁.*weights₂')*w₁*w₂*w₃/8
                    end
                end
            end
        end
    end
    return b
end

function innerproduct_R(P::UniformBSplineSpace{p,T,<:AbstractUnitRange}) where {p,T}
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    U = StaticArrays.arithmetic_closure(T)
    n = dim(P)
    A = zeros(U,n,n)
    m = factorial(2p+1)

    for q in 0:p
        a = U(eulertriangle(2p+1,q))/m
        for i in 1:n-(p-q)
            A[i,i+(p-q)] = a
        end
    end
    return Symmetric(A)
end

function innerproduct_R(P::UniformBSplineSpace{p,T}) where {p,T}
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    U = StaticArrays.arithmetic_closure(T)
    d = step(_vec(knotvector(P)))
    n = dim(P)
    A = zeros(U,n,n)
    m = factorial(2p+1)

    for q in 0:p
        a = U(eulertriangle(2p+1,q)*d)/m
        for i in 1:n-(p-q)
            A[i,i+(p-q)] = a
        end
    end
    return Symmetric(A)
end

for (fname_fit, fname_inner) in ((:fittingcontrolpoints_I, :innerproduct_I), (:fittingcontrolpoints_R, :innerproduct_R))
    # 1-dim
    @eval function $fname_fit(func, P::NTuple{1, AbstractBSplineSpace})
        Base.depwarn("`$(fname_fit)` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", fname_fit)
        P1, = P
        b = $fname_inner(func, P)
        A = $fname_inner(P1)
        return _leftdivision(A, b)
    end

    # 2-dim
    @eval function $fname_fit(func, P::NTuple{2, AbstractBSplineSpace})
        Base.depwarn("`$(fname_fit)` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", fname_fit)
        P1, P2 = P
        n1, n2 = dim.(P)
        A1, A2 = $fname_inner(P1), $fname_inner(P2)
        b = $fname_inner(func, P)
        A = [A1[i1, j1] * A2[i2, j2] for i1 in 1:n1, i2 in 1:n2, j1 in 1:n1, j2 in 1:n2]
        _A = reshape(A, n1 * n2, n1 * n2)
        _b = reshape(b, n1 * n2)
        return reshape(_leftdivision(_A, _b), n1, n2)
    end

    # 3-dim
    @eval function $fname_fit(func, P::NTuple{3, AbstractBSplineSpace})
        Base.depwarn("`$(fname_fit)` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", fname_fit)
        P1, P2, P3 = P
        n1, n2, n3 = dim(P1), dim(P2), dim(P3)
        A1, A2, A3 = $fname_inner(P1), $fname_inner(P2), $fname_inner(P3)
        b = $fname_inner(func, P)
        A = [A1[i1, j1] * A2[i2, j2] * A3[i3, j3] for i1 in 1:n1, i2 in 1:n2, i3 in 1:n3, j1 in 1:n1, j2 in 1:n2, j3 in 1:n3]
        _A = reshape(A, n1 * n2 * n3, n1 * n2 * n3)
        _b = reshape(b, n1 * n2 * n3)
        return reshape(_leftdivision(_A, _b), n1, n2, n3)
    end
end

@doc raw"""
Fitting controlpoints with least squares method.

    fittingcontrolpoints(func, Ps::Tuple)

This function will calculate ``\bm{a}_i`` to minimize the following integral.
```math
\int_I \left\|f(t)-\sum_i B_{(i,p,k)}(t) \bm{a}_i\right\|^2 dt
```
Similarly, for the two-dimensional case, minimize the following integral.
```math
\int_{I^1 \times I^2} \left\|f(t^1, t^2)-\sum_{i,j} B_{(i,p^1,k^1)}(t^1)B_{(j,p^2,k^2)}(t^2) \bm{a}_{ij}\right\|^2 dt^1dt^2
```
Currently, this function supports up to three dimensions.

# Examples
```jldoctest
julia> f(t) = SVector(cos(t),sin(t),t);

julia> P = BSplineSpace{3}(KnotVector(range(0,2π,30)) + 3*KnotVector([0,2π]));

julia> a = fittingcontrolpoints(f, P);

julia> M = BSplineManifold(a, P);

julia> norm(M(1) - f(1)) < 1e-5
true
```
"""
function fittingcontrolpoints(func, P::NTuple{Dim, AbstractBSplineSpace}) where Dim
    Base.depwarn("`fittingcontrolpoints` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :fittingcontrolpoints)
    fittingcontrolpoints_I(func,P)
end

function innerproduct_R(func, P::Vararg{AbstractBSplineSpace})
    Base.depwarn("`innterproduct_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_R)
    innerproduct_R(func,P)
end
function innerproduct_I(func, P::Vararg{AbstractBSplineSpace})
    Base.depwarn("`innterproduct_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :innterproduct_I)
    innerproduct_I(func,P)
end
function fittingcontrolpoints(func, P::Vararg{AbstractBSplineSpace})
    Base.depwarn("`fittingcontrolpoints` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :fittingcontrolpoints)
    fittingcontrolpoints_I(func,P)
end
function fittingcontrolpoints_I(func, P::Vararg{AbstractBSplineSpace})
    Base.depwarn("`fittingcontrolpoints_I` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :fittingcontrolpoints_I)
    fittingcontrolpoints_I(func,P)
end
function fittingcontrolpoints_R(func, P::Vararg{AbstractBSplineSpace})
    Base.depwarn("`fittingcontrolpoints_R` is deprecated. Please use BasicBSplineFitting.jl from the next breaking release v0.9.0", :fittingcontrolpoints_R)
    fittingcontrolpoints_R(func,P)
end
