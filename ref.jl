using Revise
using BasicBSpline
using Plots
gr()

## define
p = 1
p′ = 2

k = Knots([1,2,3.4,4,5])
k′ = Knots([-1,0.3,2,2.5,3.4,3.4,4,5.2,6])
P = BSplineSpace(p, k)
P′ = BSplineSpace(p′, k′)

k̄ = k[2:end-1]
k̄′ = k′[2:end-1]
P̄ = BSplineSpace(p-1, k̄)
P̄′ = BSplineSpace(p′-1, k̄′)

P ⊑ P′
P̄ ⊑ P̄′

## plotting

plot(
    plot([t->bsplinebasis₊₀(i,P,t) for i in 1:dim(P)], -1, 6, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(i,P′,t) for i in 1:dim(P′)], -1, 6, ylims=(0,1), legend=false),
    layout=(2,1),
    link=:x
)

plot(
    plot([t->bsplinebasis₊₀(i,P̄,t) for i in 1:dim(P̄)], -1, 6, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(i,P̄′,t) for i in 1:dim(P̄′)], -1, 6, ylims=(0,1), legend=false),
    layout=(2,1),
    link=:x
)


## refinement (p-1)
Aᵖ⁻¹ = BasicBSpline.changebasis_I(P̄, P̄′)
plot(
    plot([t->bsplinebasis₊₀(i,P̄,t) for i in 1:dim(P̄)], -1, 6, ylims=(0,1), legend=false),
    plot([t->sum(Aᵖ⁻¹[i,j]*bsplinebasis₊₀(j,P̄′,t) for j in 1:dim(P̄′)) for i in 1:dim(P̄)], -1, 6, ylims=(0,1), legend=false),
    layout=(2,1),
    link=:x
)

## refinement (p)
Aᵖ_true = [1.607 0.8215 0.3214 0 0 0;-0.607 0.1785 0.6786 1 0.5 -1;0 0 0 0 0.5 2]
plot(
    plot([t->bsplinebasis₊₀(i,P,t) for i in 1:dim(P)], -1, 6, ylims=(0,1), legend=false),
    plot([t->sum(Aᵖ_true[i,j]*bsplinebasis₊₀(j,P′,t) for j in 1:dim(P′)) for i in 1:dim(P)], -1, 6, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(j,P′,t) for j in 1:dim(P′)], -1, 6, ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)

## refinement (p)
n = dim(P)
n′ = dim(P′)

K′ = [k′[i+p′]-k′[i] for i ∈ 1:n′+1]
K = [ifelse(k[i+p]≠k[i], 1/(k[i+p]-k[i]), 0.0) for i ∈ 1:n+1]

Δ = vcat(
    (p/p′)*[K′[j+1]*(-K[i+1]*Aᵖ⁻¹[i,j]) for i ∈ 1:1, j ∈ 1:n′-1],
    (p/p′)*[K′[j+1]*(K[i]*Aᵖ⁻¹[i-1,j]-K[i+1]*Aᵖ⁻¹[i,j]) for i ∈ 2:n-1, j ∈ 1:n′-1],
    (p/p′)*[K′[j+1]*(K[i]*Aᵖ⁻¹[i-1,j]) for i ∈ n:n, j ∈ 1:n′-1]
)

Ap = zeros(3,6)
Ap[:,end] = [0,-1,2]
for j in 1:5
    Ap[:,end-j]=Ap[:,end-j+1]-Δ[:,end-j+1]
end
Ap

plot(
    plot([t->bsplinebasis₊₀(i,P,t) for i in 1:dim(P)], -1, 6, ylims=(0,1), legend=false),
    plot([t->sum(Ap[i,j]*bsplinebasis₊₀(j,P′,t) for j in 1:dim(P′)) for i in 1:dim(P)], -1, 6, ylims=(0,1), legend=false),
    plot([t->bsplinebasis₊₀(j,P′,t) for j in 1:dim(P′)], -1, 6, ylims=(0,1), legend=false),
    layout=(3,1),
    link=:x
)

## initial value

たぶん,足して1が使える!
