function bsplinebasis(i::Int,P::FastBSplineSpace{0},t::Real)
    k_1, k_2, k_end = P.vector[i], P.vector[i+1], P.vector[end]
    B_1 = Float64((k_1 ≤ t < k_2) | (k_1 ≤ t ≤ k_2 == k_end))

    return B_1
end

function bsplinebasis(i::Int,P::FastBSplineSpace{1},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[end]
    B_1, B_2 = Float64(k_1 ≤ t < k_2), Float64((k_2 ≤ t < k_3) | (k_2 ≤ t ≤ k_3 == k_end))

    K_1, K_2 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2)
    B_1 = K_1*B_1+(1-K_2)*B_2

    return B_1
end

function bsplinebasis(i::Int,P::FastBSplineSpace{2},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_4, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[i+3], P.vector[end]
    B_1, B_2, B_3 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64((k_3 ≤ t < k_4) | (k_3 ≤ t ≤ k_4 == k_end))

    K_1, K_2, K_3 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3)
    B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

    K_1, K_2 = (t-k_1)÷(k_3-k_1), (t-k_2)÷(k_4-k_2)
    B_1 = K_1*B_1+(1-K_2)*B_2

    return B_1
end

function bsplinebasis(i::Int,P::FastBSplineSpace{3},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_4, k_5, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[i+3], P.vector[i+4], P.vector[end]
    B_1, B_2, B_3, B_4 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64(k_3 ≤ t < k_4), Float64((k_4 ≤ t < k_5) | (k_4 ≤ t ≤ k_5 == k_end))

    K_1, K_2, K_3, K_4 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3), (t-k_4)÷(k_5-k_4)
    B_1, B_2, B_3 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3, K_3*B_3+(1-K_4)*B_4

    K_1, K_2, K_3 = (t-k_1)÷(k_3-k_1), (t-k_2)÷(k_4-k_2), (t-k_3)÷(k_5-k_3)
    B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

    K_1, K_2 = (t-k_1)÷(k_4-k_1), (t-k_2)÷(k_5-k_2)
    B_1 = K_1*B_1+(1-K_2)*B_2

    return B_1
end

function bsplinebasis′(i::Int,P::FastBSplineSpace{0},t::Real)
    return 0.0
end

function bsplinebasis′(i::Int,P::FastBSplineSpace{1},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[end]
    B_1, B_2 = Float64(k_1 ≤ t < k_2), Float64((k_2 ≤ t < k_3) | (k_2 ≤ t ≤ k_3 == k_end))

    K_1, K_2 = 1÷(k_2-k_1), 1÷(k_3-k_2)
    B_1 = K_1*B_1-K_2*B_2

    return B_1
end

function bsplinebasis′(i::Int,P::FastBSplineSpace{2},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_4, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[i+3], P.vector[end]
    B_1, B_2, B_3 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64((k_3 ≤ t < k_4) | (k_3 ≤ t ≤ k_4 == k_end))

    K_1, K_2, K_3 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3)
    B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

    K_1, K_2 = 1÷(k_3-k_1), 1÷(k_4-k_2)
    B_1 = 2*(K_1*B_1-K_2*B_2)

    return B_1
end

function bsplinebasis′(i::Int,P::FastBSplineSpace{3},t::Real)
    ÷(a,b) = ifelse(b==0, 0.0, a/b)
    k_1, k_2, k_3, k_4, k_5, k_end = P.vector[i], P.vector[i+1], P.vector[i+2], P.vector[i+3], P.vector[i+4], P.vector[end]
    B_1, B_2, B_3, B_4 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64(k_3 ≤ t < k_4), Float64((k_4 ≤ t < k_5) | (k_4 ≤ t ≤ k_5 == k_end))

    K_1, K_2, K_3, K_4 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3), (t-k_4)÷(k_5-k_4)
    B_1, B_2, B_3 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3, K_3*B_3+(1-K_4)*B_4

    K_1, K_2, K_3 = (t-k_1)÷(k_3-k_1), (t-k_2)÷(k_4-k_2), (t-k_3)÷(k_5-k_3)
    B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

    K_1, K_2 = 1÷(k_4-k_1), 1÷(k_5-k_2)
    B_1 = 3*(K_1*B_1-K_2*B_2)

    return B_1
end


for p in 0:MAX_DEGREE
    @eval begin

        function bsplinesupport(i::Int64, P::FastBSplineSpace{$p})
            k = P.vector
            return k[i]..k[i+$p+1]
        end

        function bsplinesupport(P::FastBSplineSpace{$p})
            k = P.vector
            return [k[i]..k[i+$p+1] for i ∈ 1:dim(P)]
        end

    end
end
