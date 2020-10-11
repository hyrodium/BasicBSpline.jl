# Numerical Integration

"""
Integrate on interval
"""
function integrate(func, I1::ClosedInterval{<:Real}, nip1, nodes1, weights1)
    dnodes1 = (width(I1) * nodes1 .+ sum(extrema(I1))) / 2

    S1 = weights1[1] * func(dnodes1[1])
    for i1 in 2:nip1
        S1 += weights1[i1] * func(dnodes1[i1])
    end

    return S1 * width(I1) / 2
end

function integrate(func, I1::ClosedInterval{<:Real}, I2::ClosedInterval{<:Real}, nip1, nip2, nodes1, nodes2, weights1, weights2)
    dnodes1 = (width(I1) * nodes1 .+ sum(extrema(I1))) / 2
    dnodes2 = (width(I2) * nodes2 .+ sum(extrema(I2))) / 2

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

    return S1 * width(I1) * width(I2) / 4
end

function integrate(func, I1::ClosedInterval{<:Real}, I2::ClosedInterval{<:Real}, I3::ClosedInterval{<:Real}, nip1, nip2, nip3, nodes1, nodes2, nodes3, weights1, weights2, weights3)
    dnodes1 = (width(I1) * nodes1 .+ sum(extrema(I1))) / 2
    dnodes2 = (width(I2) * nodes2 .+ sum(extrema(I2))) / 2
    dnodes3 = (width(I3) * nodes3 .+ sum(extrema(I3))) / 2

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

    return S1 * width(I1) * width(I2)  * width(I3) / 8
end
