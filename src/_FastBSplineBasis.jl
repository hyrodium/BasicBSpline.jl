# Faster B-spline basis function

function _s(c::Char, i::Int)
    string(c) * "_" * string(i)
end

function _code_K(p, q)
    join([_s('K', j) for j in 1:p-q+2], ',') *
    " = " *
    join(["(t-" * _s('k', j) * ")÷(" * _s('k', j + q) * "-" * _s('k', j) * ")" for j in 1:p-q+2], ',')
end

function _code_K′(p)
    join([_s('K', j) for j in 1:2], ',') *
    " = " *
    join(["1÷(" * _s('k', j + p) * "-" * _s('k', j) * ")" for j in 1:2], ',')
end

function _code_B(p, q)
    join([_s('B', j) for j in 1:p-q+1], ',') *
    " = " *
    join([_s('K', j) * "*" * _s('B', j) * "+(1-" * _s('K', j + 1) * ")*" * _s('B', j + 1) for j in 1:p-q+1], ',')
end

for p in 0:MAX_DEGREE
    @eval function bsplinebasis₊₀(P::FastBSplineSpace{$p}, i::Integer, t::Real)
        ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)
        $(Meta.parse(join([_s('k', j) for j in 1:p+2], ',') * " = " * join(["P.knots[i+" * string(j) * "]" for j in 0:p+1], ',')))
        $(Meta.parse(
            join([_s('B', j) for j in 1:p+1], ',') * " = " * join(["Float64(" * _s('k', j) * " ≤ t < " * _s('k', j + 1) * ")" for j in 1:p+1], ','),
        ))
        $(Meta.parse("begin " * join([_code_K(p, q) * "\n" * _code_B(p, q) for q in 1:p], '\n') * " end"))
        return B_1
    end
    @eval function bsplinebasis₋₀(P::FastBSplineSpace{$p}, i::Integer, t::Real)
        ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)
        $(Meta.parse(join([_s('k', j) for j in 1:p+2], ',') * " = " * join(["P.knots[i+" * string(j) * "]" for j in 0:p+1], ',')))
        $(Meta.parse(
            join([_s('B', j) for j in 1:p+1], ',') * " = " * join(["Float64(" * _s('k', j) * " < t ≤ " * _s('k', j + 1) * ")" for j in 1:p+1], ','),
        ))
        $(Meta.parse("begin " * join([_code_K(p, q) * "\n" * _code_B(p, q) for q in 1:p], '\n') * " end"))
        return B_1
    end
    @eval function bsplinebasis(P::FastBSplineSpace{$p}, i::Integer, t::Real)
        ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)
        k_end = P.knots[end]
        $(Meta.parse(join([_s('k', j) for j in 1:p+2], ',') * " = " * join(["P.knots[i+" * string(j) * "]" for j in 0:p+1], ',')))
        $(Meta.parse(
            join([_s('B', j) for j in 1:p+1], ',') * " = " * join(["Float64(" * _s('k', j) * " ≤ t < " * _s('k', j + 1) * ")" for j in 1:p+1], ','),
        ))
        $(Meta.parse(_s('B', p + 1) * " += Float64(t == " * _s('k', p + 2) * " == k_end)"))
        $(Meta.parse("begin " * join([_code_K(p, q) * "\n" * _code_B(p, q) for q in 1:p], '\n') * " end"))
        return B_1
    end
end

bsplinebasis′₊₀(P::FastBSplineSpace{0}, i::Integer, t::Real) = 0.0
bsplinebasis′₋₀(P::FastBSplineSpace{0}, i::Integer, t::Real) = 0.0
bsplinebasis′(P::FastBSplineSpace{0}, i::Integer, t::Real) = 0.0

for p in 1:MAX_DEGREE
    @eval function bsplinebasis′₊₀(P::FastBSplineSpace{$p}, i::Integer, t::Real)
        ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)
        $(Meta.parse(join([_s('k', j) for j in 1:p+2], ',') * " = " * join(["P.knots[i+" * string(j) * "]" for j in 0:p+1], ',')))
        $(Meta.parse(
            join([_s('B', j) for j in 1:p+1], ',') * " = " * join(["Float64(" * _s('k', j) * " ≤ t < " * _s('k', j + 1) * ")" for j in 1:p+1], ','),
        ))
        $(Meta.parse("begin " * join([_code_K(p, q) * "\n" * _code_B(p, q) for q in 1:p-1], '\n') * " end"))
        $(Meta.parse(_code_K′(p)))
        B_1 = $p * (K_1 * B_1 - K_2 * B_2)
        return B_1
    end
    @eval function bsplinebasis′₋₀(P::FastBSplineSpace{$p}, i::Integer, t::Real)
        ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)
        $(Meta.parse(join([_s('k', j) for j in 1:p+2], ',') * " = " * join(["P.knots[i+" * string(j) * "]" for j in 0:p+1], ',')))
        $(Meta.parse(
            join([_s('B', j) for j in 1:p+1], ',') * " = " * join(["Float64(" * _s('k', j) * " < t ≤ " * _s('k', j + 1) * ")" for j in 1:p+1], ','),
        ))
        $(Meta.parse("begin " * join([_code_K(p, q) * "\n" * _code_B(p, q) for q in 1:p-1], '\n') * " end"))
        $(Meta.parse(_code_K′(p)))
        B_1 = $p * (K_1 * B_1 - K_2 * B_2)
        return B_1
    end
    # @eval function bsplinebasis′(P::FastBSplineSpace{$p}, i::Integer, t::Real)
    #     ÷(a,b) = ifelse(b == 0.0, 0.0, a/b)
    #     k_end = P.knots[end]
    #     $(Meta.parse(join([_s('k',j) for j in 1:p+2],',')*" = "*join(["P.knots[i+"*string(j)*"]" for j in 0:p+1],',')))
    #     $(Meta.parse(join([_s('B',j) for j in 1:p+1],',')*" = "*join(["Float64("*_s('k',j)*" ≤ t < "*_s('k',j+1)*")" for j in 1:p+1],',')))
    #     $(Meta.parse(_s('B',p+1)*" += Float64(t == "*_s('k',p+2)*" == k_end)"))
    #     $(Meta.parse("begin "*join([_code_K(p,q)*"\n"*_code_B(p,q) for q in 1:p-1],'\n')*" end"))
    #     $(Meta.parse(_code_K′(p)))
    #     B_1 = $p*(K_1*B_1-K_2*B_2)
    #     return B_1
    # end
end

function bsplinebasis′(P::FastBSplineSpace, i::Integer, t::Real)::Float64
    # TODO: faster
    ÷(a, b) = ifelse(b == 0.0, 0.0, a / b)

    p = degree(P)
    k = knots(P)

    B_1, B_2 = bsplinebasis(FastBSplineSpace(p-1, k), i, t), bsplinebasis(FastBSplineSpace(p-1, k), i+1, t)

    K_1, K_2 = 1 ÷ (k[i+p] - k[i]), 1 ÷ (k[i+p+1] - k[i+1])
    B_1 = p * (K_1 * B_1 - K_2 * B_2)

    return B_1
end

# TODO: use @eval macro to produce code
"""
Returns the value of ``B_{(i,0,k)}(t)``.
Assumption:
* ```k_{i} ≤ t < k_{i+1}``
"""
function _bsb0(::Union{Vector{<:Real},Knots}, ::Real, ::Integer)
    return 1.0
end

"""
Returns the values of ``B_{(i,1,k)}(t), B_{(i+1,1,k)}(t)``.
Assumption:
* ``k_{i} ≤ t < k_{i+1}``
"""
function _bsb1(k::Union{Vector{<:Real},Knots}, t::Real, i::Integer)
    B1 = (k[i+1]-t)/(k[i+1]-k[i])
    B2 = (t-k[i])/(k[i+1]-k[i])
    return B1, B2
end

"""
Returns the values of ``B_{(i,2,k)}(t), B_{(i+1,2,k)}(t), B_{(i+2,2,k)}(t)``.
Assumption:
* ``k_{i} ≤ t < k_{i+1}``
"""
@inline function _bsb2(k::Union{Vector{<:Real},Knots}, t::Real, i::Integer)
    B = _bsb1(k, t, i)

    B1 = (k[i+1]-t)/(k[i+1]-k[i-1]) * B[1]
    B2 = (t-k[i-1])/(k[i+1]-k[i-1]) * B[1] +
         (k[i+2]-t)/(k[i+2]-k[i]) * B[2]
    B3 = (t-k[i])/(k[i+2]-k[i]) * B[2]
    return B1, B2, B3
end

"""
Returns the values of ``B_{(i,3,k)}(t), B_{(i+1,3,k)}(t), B_{(i+2,3,k)}(t), B_{(i+3,3,k)}(t)``.
Assumption:
* ``k_{i} ≤ t < k_{i+1}``
"""
@inline function _bsb3(k::Union{Vector{<:Real},Knots}, t::Real, i::Integer)
    B = _bsb2(k, t, i)

    B1 = (k[i+1]-t)/(k[i+1]-k[i-2]) * B[1]
    B2 = (t-k[i-2])/(k[i+1]-k[i-2]) * B[1] +
         (k[i+2]-t)/(k[i+2]-k[i-1]) * B[2]
    B3 = (t-k[i-1])/(k[i+2]-k[i-1]) * B[2] +
         (k[i+3]-t)/(k[i+3]-k[i]) * B[3]
    B4 = (t-k[i])/(k[i+3]-k[i]) * B[3]
    return B1, B2, B3, B4
end

"""
Returns the values of ``B_{(i,4,k)}(t), B_{(i+1,4,k)}(t), B_{(i+2,4,k)}(t), B_{(i+3,4,k)}(t), B_{(i+4,4,k)}(t)``.
Assumption:
* ``k_{i} ≤ t < k_{i+1}``
"""
@inline function _bsb4(k::Union{Vector{<:Real},Knots}, t::Real, i::Integer)
    B = _bsb3(k, t, i)

    B1 = (k[i+1]-t)/(k[i+1]-k[i-3]) * B[1]
    B2 = (t-k[i-3])/(k[i+1]-k[i-3]) * B[1] +
         (k[i+2]-t)/(k[i+2]-k[i-2]) * B[2]
    B3 = (t-k[i-2])/(k[i+2]-k[i-2]) * B[2] +
         (k[i+3]-t)/(k[i+3]-k[i-1]) * B[3]
    B4 = (t-k[i-1])/(k[i+3]-k[i-1]) * B[3] +
         (k[i+4]-t)/(k[i+4]-k[i]) * B[4]
    B5 = (t-k[i])/(k[i+4]-k[i]) * B[4]
    return B1, B2, B3, B4, B5
end

"""
Returns the values of ``B_{(i,5,k)}(t), B_{(i+1,5,k)}(t), B_{(i+2,5,k)}(t), B_{(i+3,5,k)}(t), B_{(i+4,5,k)}(t), B_{(i+5,5,k)}(t)``.
Assumption:
* ``k_{i} ≤ t < k_{i+1}``
"""
@inline function _bsb5(k::Union{Vector{<:Real},Knots}, t::Real, i::Integer)
    B = _bsb4(k, t, i)

    B1 = (k[i+1]-t)/(k[i+1]-k[i-4]) * B[1]
    B2 = (t-k[i-4])/(k[i+1]-k[i-4]) * B[1] +
         (k[i+2]-t)/(k[i+2]-k[i-3]) * B[2]
    B3 = (t-k[i-3])/(k[i+2]-k[i-3]) * B[2] +
         (k[i+3]-t)/(k[i+3]-k[i-2]) * B[3]
    B4 = (t-k[i-2])/(k[i+3]-k[i-2]) * B[3] +
         (k[i+4]-t)/(k[i+4]-k[i-1]) * B[4]
    B5 = (t-k[i-1])/(k[i+4]-k[i-1]) * B[4] +
         (k[i+5]-t)/(k[i+5]-k[i]) * B[5]
    B6 = (t-k[i])/(k[i+5]-k[i]) * B[5]
    return B1, B2, B3, B4, B5, B6
end

@inline function _bsplinebasis(P::FastBSplineSpace{0}, t::Real, i::Integer)
    return _bsb0(P.knots,t,i)
end

@inline function _bsplinebasis(P::FastBSplineSpace{1}, t::Real, i::Integer)
    return _bsb1(P.knots,t,i)
end

@inline function _bsplinebasis(P::FastBSplineSpace{2}, t::Real, i::Integer)
    return _bsb2(P.knots,t,i)
end

@inline function _bsplinebasis(P::FastBSplineSpace{3}, t::Real, i::Integer)
    return _bsb3(P.knots,t,i)
end

@inline function _bsplinebasis(P::FastBSplineSpace{4}, t::Real, i::Integer)
    return _bsb4(P.knots,t,i)
end

@inline function _bsplinebasis(P::FastBSplineSpace{5}, t::Real, i::Integer)
    return _bsb5(P.knots,t,i)
end
