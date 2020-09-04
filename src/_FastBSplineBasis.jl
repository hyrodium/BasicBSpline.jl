function swi(c::Char, i::Int)
    string(c)*"_"*string(i)
end

function _code_K(p,q)
    join([swi('K',j) for j in 1:p-q+2],',')*" = "*join(["(t-"*swi('k',j)*")÷("*swi('k',j+q)*"-"*swi('k',j)*")" for j in 1:p-q+2],',')
end

function _code_K′(p)
    join([swi('K',j) for j in 1:2],',')*" = "*join(["1÷("*swi('k',j+p)*"-"*swi('k',j)*")" for j in 1:2],',')
end

function _code_B(p,q)
    join([swi('B',j) for j in 1:p-q+1],',')*" = "*join([swi('K',j)*"*"*swi('B',j)*"+(1-"*swi('K',j+1)*")*"*swi('B',j+1) for j in 1:p-q+1],',')
end

for p in 0:MAX_DEGREE
    @eval function bsplinebasis₊₀(i::Integer, P::FastBSplineSpace{$p}, t::Real)
        ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
        $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
        $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" ≤ t < "*swi('k',j+1)*")" for j in 1:p+1],',')))
        $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p],'\n')*" end"))
        return B_1
    end
    @eval function bsplinebasis₋₀(i::Integer, P::FastBSplineSpace{$p}, t::Real)
        ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
        $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
        $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" < t ≤ "*swi('k',j+1)*")" for j in 1:p+1],',')))
        $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p],'\n')*" end"))
        return B_1
    end
    @eval function bsplinebasis(i::Integer, P::FastBSplineSpace{$p}, t::Real)
        ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
        k_end = P.knotvector[end]
        $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
        $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" ≤ t < "*swi('k',j+1)*")" for j in 1:p+1],',')))
        $(Meta.parse(swi('B',p+1)*" += Float64(t == "*swi('k',p+2)*" == k_end)"))
        $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p],'\n')*" end"))
        return B_1
    end
end

function bsplinebasis′₊₀(i::Integer,P::FastBSplineSpace{0},t::Real)
    return 0.0
end
function bsplinebasis′₋₀(i::Integer,P::FastBSplineSpace{0},t::Real)
    return 0.0
end
function bsplinebasis′(i::Integer,P::FastBSplineSpace{0},t::Real)
    return 0.0
end

for p in 1:MAX_DEGREE
    @eval function bsplinebasis′₊₀(i::Integer, P::FastBSplineSpace{$p}, t::Real)
        ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
        $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
        $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" ≤ t < "*swi('k',j+1)*")" for j in 1:p+1],',')))
        $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p-1],'\n')*" end"))
        $(Meta.parse(_code_K′(p)))
        B_1 = $(p)*(K_1*B_1-K_2*B_2)
        return B_1
    end
    @eval function bsplinebasis′₋₀(i::Integer, P::FastBSplineSpace{$p}, t::Real)
        ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
        $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
        $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" < t ≤ "*swi('k',j+1)*")" for j in 1:p+1],',')))
        $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p-1],'\n')*" end"))
        $(Meta.parse(_code_K′(p)))
        B_1 = $(p)*(K_1*B_1-K_2*B_2)
        return B_1
    end
    # @eval function bsplinebasis′(i::Integer, P::FastBSplineSpace{$p}, t::Real)
    #     ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
    #     k_end = P.knotvector[end]
    #     $(Meta.parse(join([swi('k',j) for j in 1:p+2],',')*" = "*join(["P.knotvector[i+"*string(j)*"]" for j in 0:p+1],',')))
    #     $(Meta.parse(join([swi('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*swi('k',j)*" ≤ t < "*swi('k',j+1)*")" for j in 1:p+1],',')))
    #     $(Meta.parse(swi('B',p+1)*" += Float64(t == "*swi('k',p+2)*" == k_end)"))
    #     $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p-1],'\n')*" end"))
    #     B_1 = $(p)*(K_1*B_1-K_2*B_2)
    #     return B_1
    # end
end

"""
TODO: faster.....
"""
function bsplinebasis′(i::Integer, P::FastBSplineSpace, t::Real)::Float64
    ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)

    p = degree(P)
    k = knots(P)

    B_1, B_2 = bsplinebasis(i, FastBSplineSpace(p-1, k), t), bsplinebasis(i+1, FastBSplineSpace(p-1, k), t)

    K_1, K_2 = 1÷(k[i+p]-k[i]), 1÷(k[i+p+1]-k[i+1])
    B_1 = p*(K_1*B_1-K_2*B_2)

    return B_1
end
