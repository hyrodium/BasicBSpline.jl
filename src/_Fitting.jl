# Numerical Integration

"""
fast, but only for 1-dim
"""
function GaussianQuadrature_1dim(func::Function, D::Array{ClosedInterval{T},1} where T<:Real , nodes, weights)
    d = length(D) # d must be equal to 2
    nip = length(nodes)
    dnodes = [(width(I1)*nodes.+sum(extrema(I1)))/2 for I1 in D]
    widths = width.(D)
    S = zero(func(mean.(D)))
    for i1 ∈ 1:nip
        S += weights[i1]*func([dnodes[1][i1]])
    end
    return S*prod(widths)/2^d
end

"""
fast, but only for 2-dim
"""
function GaussianQuadrature_2dim(func::Function, D::Array{ClosedInterval{T},1} where T<:Real , nodes, weights)
    d = length(D) # d must be equal to 2
    nip = length(nodes)
    dnodes = [(width(I1)*nodes.+sum(extrema(I1)))/2 for I1 in D]
    widths = width.(D)
    S = zero(func(mean.(D)))
    for i1 ∈ 1:nip, i2 ∈ 1:nip
        S += weights[i1]*weights[i2]*func([dnodes[1][i1],dnodes[2][i2]])
    end
    return S*prod(widths)/2^d
end

function fittingcontrolpoints_1dim(func::Function, Ps::Array{T,1} where T<:AbstractBSplineSpace)
    d = 1 # length(Ps)
    k = knots.(Ps)
    p = degree.(Ps)
    nip = maximum(p) + 1
    nodes, weights = gausslegendre(nip)
    function α(I,J)
        min = [minimum([I[i], J[i]]) for i in 1:d]
        max = [maximum([I[i]+degree(Ps[i]), J[i]+degree(Ps[i])]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = 0.0
            for i1 in rng[1]
                S += GaussianQuadrature_1dim(t->(bsplinebasis(I[1],Ps[1],t[1]) * bsplinebasis(J[1],Ps[1],t[1])), [k[1][i1]..k[1][i1+1]], nodes, weights)
            end
            return S
        end
    end
    function β(I)
        min = [I[i] for i in 1:d]
        max = [I[i]+degree(Ps[i]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = zero(func([0.0]))
            for i1 in rng[1]
                S += GaussianQuadrature_1dim(t->(bsplinebasis(I[1],Ps[1],t[1]) * func(t)), [k[1][i1]..k[1][i1+1]], nodes, weights)
            end
            return S
        end
    end
    n = dim.(Ps)
    X = Array{Float64}(undef, n...)
    A = [α(I,J) for I in CartesianIndices(X), J in CartesianIndices(X)]
    b = [β(I) for I in CartesianIndices(X)]
    A = reshape(A,prod(n),prod(n))
    b = reshape(b,prod(n))
    return reshape(inv(A)*b, n...)
end


function fittingcontrolpoints_2dim(func::Function, Ps::Array{T,1} where T<:AbstractBSplineSpace)
    d = 2 # length(Ps)
    k = knots.(Ps)
    p = degree.(Ps)
    nip = maximum(p) + 1
    nodes, weights = gausslegendre(nip)
    function α(I,J)
        min = [minimum([I[i], J[i]]) for i in 1:d]
        max = [maximum([I[i]+degree(Ps[i]), J[i]+degree(Ps[i])]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = 0.0
            for i1 in rng[1], i2 in rng[2]
                S += GaussianQuadrature_2dim(t->(bsplinebasis(I[1],Ps[1],t[1]) * bsplinebasis(I[2],Ps[2],t[2]) * bsplinebasis(J[1],Ps[1],t[1]) * bsplinebasis(J[2],Ps[2],t[2])), [k[1][i1]..k[1][i1+1], k[2][i2]..k[2][i2+1]], nodes, weights)
            end
            return S
        end
    end
    function β(I)
        min = [I[i] for i in 1:d]
        max = [I[i]+degree(Ps[i]) for i in 1:d]
        rng = [min[i]:max[i] for i in 1:d]
        if prod(length.(rng)) == 0
            return 0.0
        else
            S = zero(func([0,0.0]))
            for i1 in rng[1], i2 in rng[2]
                S += GaussianQuadrature_2dim(t->(bsplinebasis(I[1],Ps[1],t[1]) * bsplinebasis(I[2],Ps[2],t[2]) * func(t)), [k[1][i1]..k[1][i1+1], k[2][i2]..k[2][i2+1]], nodes, weights)
            end
            return S
        end
    end
    n = dim.(Ps)
    X = Array{Float64}(undef, n...)
    A = [α(I,J) for I in CartesianIndices(X), J in CartesianIndices(X)]
    b = [β(I) for I in CartesianIndices(X)]
    A = reshape(A,prod(n),prod(n))
    b = reshape(b,prod(n))
    return reshape(inv(A)*b, n...)
end


"""
Approximate given function by linear combination of B-spline functions.
This function returns its control points.
TODO: currently, this function only supports for 1-dim and 2-dim B-spline manifold.
"""
function fittingcontrolpoints(func::Function, Ps::Array{T,1} where T<:AbstractBSplineSpace)
    d = length(Ps)

    if d == 1
        a = fittingcontrolpoints_1dim(func, Ps)
    elseif d == 2
        a = fittingcontrolpoints_2dim(func, Ps)
    end

    return a
end
