# Fitting

"""
Assumption:
* i ≤ j
* 1 ≤ n ≤ p-j+i+1
"""
function _bsplineintegrate(P::AbstractBSplineSpace, i, j, n, nip, nodes, weights)
    p = degree(P)
    I = knots(P)[j+n-1]..knots(P)[j+n]

    f(t) = bsplinebasis(i, P, t) * bsplinebasis(j, P, t)
    return integrate(f, I, nip, nodes, weights)
end

function _bsplineintegrate(P::AbstractBSplineSpace, i, j, nip, nodes, weights)
    p = degree(P)
    k = knots(P)
    Δ = j - i
    if Δ < -p
        return 0.0
    elseif Δ > p
        return 0.0
    elseif Δ ≥ 0
        s = 0.0
        for n in 1:p-j+i+1
            s += _bsplineintegrate(P, i, j, n, nip, nodes, weights)
        end
        return s
    else
        s = 0.0
        for n in 1:p-i+j+1
            s += _bsplineintegrate(P, j, i, n, nip, nodes, weights)
        end
        return s
    end
end

function innerproduct(P::AbstractBSplineSpace)
    p = degree(P)
    n = dim(P)
    nip = p + 1
    nodes, weights = gausslegendre(nip)
    return [_bsplineintegrate(P, i, j, nip, nodes, weights) for i in 1:n, j in 1:n]
end

function fittingcontrolpoints_1dim(func::Function, Ps::Array{<:AbstractBSplineSpace,1})
    d = 1 # length(Ps)
    k = knots.(Ps)
    p = degree.(Ps)
    nip = maximum(p) + 1
    nodes, weights = gausslegendre(nip)
    function β(I)
        min = [I[i] for i in 1:d]
        max = [I[i] + degree(Ps[i]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = zero(func([0.0]))
            for i1 in rng[1]
                S += GaussianQuadrature_1dim(t -> (bsplinebasis(I[1], Ps[1], t[1]) * func(t)), [k[1][i1]..k[1][i1+1]], nodes, weights)
            end
            return S
        end
    end
    n = dim.(Ps)
    X = Array{Float64}(undef, n...)
    b = [β(I) for I in CartesianIndices(X)]
    A = innerproduct(Ps[1])
    b = reshape(b, prod(n))
    return reshape(inv(A) * b, n...)
end


function fittingcontrolpoints_2dim(func::Function, Ps::Array{<:AbstractBSplineSpace,1})
    d = 2 # length(Ps)
    k = knots.(Ps)
    p = degree.(Ps)
    nip = maximum(p) + 1
    nodes, weights = gausslegendre(nip)
    function β(I)
        min = [I[i] for i in 1:d]
        max = [I[i] + degree(Ps[i]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = zero(func([0, 0.0]))
            for i1 in rng[1], i2 in rng[2]
                S += GaussianQuadrature_2dim(
                    t -> (bsplinebasis(I[1], Ps[1], t[1]) * bsplinebasis(I[2], Ps[2], t[2]) * func(t)),
                    [k[1][i1]..k[1][i1+1], k[2][i2]..k[2][i2+1]],
                    nodes,
                    weights,
                )
            end
            return S
        end
    end
    n = dim.(Ps)
    X = Array{Float64}(undef, n...)
    b = [β(I) for I in CartesianIndices(X)]
    n1, n2 = dim.(Ps)
    A1, A2 = innerproduct.(Ps)
    A = [A1[i1, j1] * A2[i2, j2] for i1 in 1:n1, i2 in 1:n2, j1 in 1:n1, j2 in 1:n2]
    A = reshape(A, prod(n), prod(n))
    b = reshape(b, prod(n))
    return reshape(inv(A) * b, n...)
end


"""
Approximate given function by linear combination of B-spline functions.
This function returns its control points.
"""
function fittingcontrolpoints(func::Function, Ps::Array{<:AbstractBSplineSpace,1})
    # TODO: currently, this function only supports for 1-dim and 2-dim B-spline manifold.
    d = length(Ps)
    if d == 1
        return fittingcontrolpoints_1dim(func, Ps)
    elseif d == 2
        return fittingcontrolpoints_2dim(func, Ps)
    end
end
