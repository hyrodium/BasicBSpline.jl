# Numerical Integration

"""
Integrate on interval
"""
function integrate(func::Function, I::ClosedInterval{<:Real}, nip, nodes, weights)
    dnodes = (width(I) * nodes .+ sum(extrema(I))) / 2
    S = weights[1] * func(dnodes[1])
    for i in 2:nip
        S += weights[i] * func(dnodes[i])
    end
    return S * width(I) / 2
end

"""
fast, but only for 1-dim
"""
function GaussianQuadrature_1dim(func::Function, D, nodes, weights)
    d = length(D) # d must be equal to 2
    nip = length(nodes)
    dnodes = [(width(I1) * nodes .+ sum(extrema(I1))) / 2 for I1 in D]
    widths = width.(D)
    S = zero(func(mean.(D)))
    for i1 in 1:nip
        S += weights[i1] * func([dnodes[1][i1]])
    end
    return S * prod(widths) / 2^d
end

"""
fast, but only for 2-dim
"""
function GaussianQuadrature_2dim(func::Function, D, nodes, weights)
    d = length(D) # d must be equal to 2
    nip = length(nodes)
    dnodes = [(width(I1) * nodes .+ sum(extrema(I1))) / 2 for I1 in D]
    widths = width.(D)
    S = zero(func(mean.(D)))
    for i1 in 1:nip, i2 in 1:nip
        S += weights[i1] * weights[i2] * func([dnodes[1][i1], dnodes[2][i2]])
    end
    return S * prod(widths) / 2^d
end
