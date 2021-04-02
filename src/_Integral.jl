# Numerical Integration

struct GaussLegendre
    nip::Int  # number of integral points
    nodes::Vector{Float64}
    weights::Vector{Float64}
    function GaussLegendre(nip::Int)
        nodes, weights = gausslegendre(nip)
        return new(nip, nodes, weights)
    end
end

"""
Integrate on interval
"""
function integrate(func, I1::ClosedInterval{<:Real}, gl1)
    nip1, nodes1, weights1 = gl1.nip, gl1.nodes, gl1.weights
    dnodes1 = (IntervalSets.width(I1) * nodes1 .+ sum(extrema(I1))) / 2

    S1 = weights1[1] * func(dnodes1[1])
    for i1 in 2:nip1
        S1 += weights1[i1] * func(dnodes1[i1])
    end

    return S1 * IntervalSets.width(I1) / 2
end

"""
Integrate on rectangular region.
"""
function integrate(func, I1::ClosedInterval{<:Real}, I2::ClosedInterval{<:Real}, gl1, gl2)
    nip1, nodes1, weights1 = gl1.nip, gl1.nodes, gl1.weights
    nip2, nodes2, weights2 = gl2.nip, gl2.nodes, gl2.weights
    dnodes1 = (IntervalSets.width(I1) * nodes1 .+ sum(extrema(I1))) / 2
    dnodes2 = (IntervalSets.width(I2) * nodes2 .+ sum(extrema(I2))) / 2

    S2 = weights2[1] * func(dnodes1[1], dnodes2[1])
    for i2 in 2:nip2
        S2 += weights2[i2] * func(dnodes1[1], dnodes2[i2])
    end
    S1 = weights1[1] * S2
    for i1 in 2:nip1
        S2 = weights2[1] * func(dnodes1[i1], dnodes2[1])
        for i2 in 2:nip2
            S2 += weights2[i2] * func(dnodes1[i1], dnodes2[i2])
        end
        S1 += weights1[i1] * S2
    end

    return S1 * IntervalSets.width(I1) * IntervalSets.width(I2) / 4
end

"""
Integrate on rectangular solid
"""
function integrate(func, I1::ClosedInterval{<:Real}, I2::ClosedInterval{<:Real}, I3::ClosedInterval{<:Real}, gl1, gl2, gl3)
    nip1, nodes1, weights1 = gl1.nip, gl1.nodes, gl1.weights
    nip2, nodes2, weights2 = gl2.nip, gl2.nodes, gl2.weights
    nip3, nodes3, weights3 = gl3.nip, gl3.nodes, gl3.weights
    dnodes1 = (IntervalSets.width(I1) * nodes1 .+ sum(extrema(I1))) / 2
    dnodes2 = (IntervalSets.width(I2) * nodes2 .+ sum(extrema(I2))) / 2
    dnodes3 = (IntervalSets.width(I3) * nodes3 .+ sum(extrema(I3))) / 2

    S3 = weights3[1] * func(dnodes1[1],dnodes2[1],dnodes3[1])
    for i3 in 2:nip3
        S3 += weights3[i3] * func(dnodes1[1],dnodes2[1],dnodes3[i3])
    end
    S2 = weights2[1] * S3
    for i2 in 2:nip2
        S3 = weights3[1] * func(dnodes1[1],dnodes2[i2],dnodes3[1])
        for i3 in 2:nip3
            S3 += weights3[i3] * func(dnodes1[1],dnodes2[i2],dnodes3[i3])
        end
        S2 += weights2[i2] * S3
    end
    S1 = weights1[1] * S2
    for i1 in 2:nip1
        S3 = weights3[1] * func(dnodes1[i1],dnodes2[1],dnodes3[1])
        for i3 in 2:nip3
            S3 += weights3[i3] * func(dnodes1[i1],dnodes2[1],dnodes3[i3])
        end
        S2 = weights2[1] * S3
        for i2 in 2:nip2
            S3 = weights3[1] * func(dnodes1[i1],dnodes2[i2],dnodes3[1])
            for i3 in 2:nip3
                S3 += weights3[i3] * func(dnodes1[i1],dnodes2[i2],dnodes3[i3])
            end
            S2 += weights2[i2] * S3
        end
        S1 += weights1[i1] * S2
    end

    return S1 * IntervalSets.width(I1) * IntervalSets.width(I2)  * IntervalSets.width(I3) / 8
end
