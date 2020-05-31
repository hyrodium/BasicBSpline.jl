using BasicBSpline
using BenchmarkTools

k=Knots(rand(8))

P=BSplineSpace(3,k)
fP=FastBSplineSpace(3,k)

dim(P)
dim(fP)

bsplinebasis(1,P,0.2)
bsplinebasis(1,fP,0.2)
BasicBSpline.bsplinebasis3(1,fP,0.2)
splinebasis(1,P,0.2)


@benchmark bsplinebasis(1,P,0.4)
@benchmark bsplinebasis(1,fP,0.4)
@benchmark BasicBSpline.bsplinebasis3(1,fP,0.4)
@benchmark splinebasis(1,P,0.4)









FastBSplineSpace{2}([1,2,3])


function tess(fP::FastBSplineSpace{3})
    bsplinebasis(1,fP,0.4)
end

@benchmark tess(fP)




@benchmark sum(collect(1:5))


@benchmark sum(1:5)


@benchmark collect(1:5)

@benchmark sum([1,2,3,4,5])


@benchmark bsplinebasis1(1,fP,0.4)+bsplinebasis1(1,fP,0.3)



function splinebasis(i::Int,P::BSplineSpace,t::Real)
    p = P.degree
    k = P.knots.vector
    ÷(a,b) = ifelse(b==0, 0.0, a/b)

    if p == 0
        k_1, k_2, k_end = k[i], k[i+1], k[end]
        B_1 = Float64((k_1 ≤ t < k_2) | (k_1 ≤ t ≤ k_2 == k_end))

        return B_1

    elseif p == 1
        k_1, k_2, k_3, k_end = k[i], k[i+1], k[i+2], k[end]
        B_1, B_2 = Float64(k_1 ≤ t < k_2), Float64((k_2 ≤ t < k_3) | (k_2 ≤ t ≤ k_3 == k_end))

        K_1, K_2 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2)
        B_1 = K_1*B_1+(1-K_2)*B_2

        return B_1

    elseif p == 2
        k_1, k_2, k_3, k_4, k_end = k[i], k[i+1], k[i+2], k[i+3], k[end]
        B_1, B_2, B_3 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64((k_3 ≤ t < k_4) | (k_3 ≤ t ≤ k_4 == k_end))

        K_1, K_2, K_3 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3)
        B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

        K_1, K_2 = (t-k_1)÷(k_3-k_1), (t-k_2)÷(k_4-k_2)
        B_1 = K_1*B_1+(1-K_2)*B_2

        return B_1

    elseif p == 3
        k_1, k_2, k_3, k_4, k_5, k_end = k[i], k[i+1], k[i+2], k[i+3], k[i+4], k[end]
        B_1, B_2, B_3, B_4 = Float64(k_1 ≤ t < k_2), Float64(k_2 ≤ t < k_3), Float64(k_3 ≤ t < k_4), Float64((k_4 ≤ t < k_5) | (k_4 ≤ t ≤ k_5 == k_end))

        K_1, K_2, K_3, K_4 = (t-k_1)÷(k_2-k_1), (t-k_2)÷(k_3-k_2), (t-k_3)÷(k_4-k_3), (t-k_4)÷(k_5-k_4)
        B_1, B_2, B_3 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3, K_3*B_3+(1-K_4)*B_4

        K_1, K_2, K_3 = (t-k_1)÷(k_3-k_1), (t-k_2)÷(k_4-k_2), (t-k_3)÷(k_5-k_3)
        B_1, B_2 = K_1*B_1+(1-K_2)*B_2, K_2*B_2+(1-K_3)*B_3

        K_1, K_2 = (t-k_1)÷(k_4-k_1), (t-k_2)÷(k_5-k_2)
        B_1 = K_1*B_1+(1-K_2)*B_2

        return B_1
    end
end
