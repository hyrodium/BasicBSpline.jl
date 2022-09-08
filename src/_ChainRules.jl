const BSPLINESPACE_INFO = """
derivatives of B-spline basis functions with respect to BSplineSpace not implemented currently.
"""

# bsplinebasis
function ChainRulesCore.frule((_, ΔP, Δi, Δt), ::typeof(bsplinebasis), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis(P,i,t)
    ∂B_∂P = @not_implemented BSPLINESPACE_INFO
    # ∂B_∂i = NoTangent()
    ∂B_∂t = bsplinebasis′(P,i,t)
    return (B, ∂B_∂P*ΔP + ∂B_∂t*Δt)
end
function ChainRulesCore.rrule(::typeof(bsplinebasis), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis(P,i,t)
    # project_t = ProjectTo(t)  # Not sure we need this ProjectTo.
    function bsplinebasis_pullback(ΔB)
        P̄ = @not_implemented BSPLINESPACE_INFO
        ī = NoTangent()
        t̄ = bsplinebasis′(P,i,t) * ΔB
        return (NoTangent(), P̄, ī, t̄)
    end
    return (B, bsplinebasis_pullback)
end

# bsplinebasis₊₀
function ChainRulesCore.frule((_, ΔP, Δi, Δt), ::typeof(bsplinebasis₊₀), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis₊₀(P,i,t)
    ∂B_∂P = @not_implemented BSPLINESPACE_INFO
    # ∂B_∂i = NoTangent()
    ∂B_∂t = bsplinebasis′₊₀(P,i,t)
    return (B, ∂B_∂P*ΔP + ∂B_∂t*Δt)
end
function ChainRulesCore.rrule(::typeof(bsplinebasis₊₀), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis₊₀(P,i,t)
    # project_t = ProjectTo(t)  # Not sure we need this ProjectTo.
    function bsplinebasis_pullback(ΔB)
        P̄ = @not_implemented BSPLINESPACE_INFO
        ī = NoTangent()
        t̄ = bsplinebasis′₊₀(P,i,t) * ΔB
        return (NoTangent(), P̄, ī, t̄)
    end
    return (B, bsplinebasis_pullback)
end

# bsplinebasis₋₀
function ChainRulesCore.frule((_, ΔP, Δi, Δt), ::typeof(bsplinebasis₋₀), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis₋₀(P,i,t)
    ∂B_∂P = @not_implemented BSPLINESPACE_INFO
    # ∂B_∂i = NoTangent()
    ∂B_∂t = bsplinebasis′₋₀(P,i,t)
    return (B, ∂B_∂P*ΔP + ∂B_∂t*Δt)
end
function ChainRulesCore.rrule(::typeof(bsplinebasis₋₀), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasis₋₀(P,i,t)
    # project_t = ProjectTo(t)  # Not sure we need this ProjectTo.
    function bsplinebasis_pullback(ΔB)
        P̄ = @not_implemented BSPLINESPACE_INFO
        ī = NoTangent()
        t̄ = bsplinebasis′₋₀(P,i,t) * ΔB
        return (NoTangent(), P̄, ī, t̄)
    end
    return (B, bsplinebasis_pullback)
end

# bsplinebasisall
function ChainRulesCore.frule((_, ΔP, Δi, Δt), ::typeof(bsplinebasisall), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasisall(P,i,t)
    ∂B_∂P = @not_implemented BSPLINESPACE_INFO
    # ∂B_∂i = NoTangent()
    ∂B_∂t = bsplinebasisall(derivative(P),i,t)
    return (B, ∂B_∂P*ΔP + ∂B_∂t*Δt)
end
function ChainRulesCore.rrule(::typeof(bsplinebasisall), P::AbstractFunctionSpace, i::Integer, t::Real)
    B = bsplinebasisall(P,i,t)
    # project_t = ProjectTo(t)  # Not sure we need this ProjectTo.
    function bsplinebasis_pullback(ΔB)
        P̄ = @not_implemented BSPLINESPACE_INFO
        ī = NoTangent()
        t̄ = bsplinebasisall(derivative(P),i,t)' * ΔB
        return (NoTangent(), P̄, ī, t̄)
    end
    return (B, bsplinebasis_pullback)
end
