# Mathematical properties of B-spline

## Introduction
[B-spline](https://en.wikipedia.org/wiki/B-spline) is a mathematical object, and it has a lot of application(e.g. [NURBS](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline), [IGA](https://en.wikipedia.org/wiki/Isogeometric_analysis)).

In this page, we'll explain the mathematical definition and property of B-spline with Julia code.

Before running the following code, do not forget importing the package:
```julia
using BasicBSpline
```


## Knot vector
A finite sequence
```math
k=(k_1, \dots, k_l)
```
is called **knot vector** if the sequence is
broad monotonic increase, i.e. (``k_{i} \le k_{i+1}``)

```julia
k=Knots([1,2,3])
```

We denotes a number of knots by *sharp* symbol like this:
```math
\sharp k=\sharp(k_1, \dots, k_l) =l
```

```julia
♯(k) # 3
length(k) # 3
```

We introduce additional operator ``+`` and product operator ``\cdot``
```math
\begin{aligned}
k^{(1)}+k^{(2)}
&=(k^{(1)}_1, \dots, k^{(1)}_l) + (k^{(2)}_1, \dots, k^{(2)}_l) \\
&=(\text{sort of union of} \  k^{(1)} \ \text{and} \  k^{(2)} \text{)} \\
nk&=\underbrace{k+\cdots+k}_{n}
\end{aligned}
```
For example, ``(1,2,3)+(2,4,5)=(1,2,2,3,5)``, ``2\cdot (2,3)=(2,2,3,3)``

```julia
k = Knots([1,2,3]) + Knots([2,4,5]) # Knots([1,2,2,3,5])
2Knots([2,3]) # Knots([2,2,3,3])
```

Unique operator
```math
\begin{aligned}
\widehat{k}
&=(\text{remove duplicates of} \  k) \\
\end{aligned}
```
For example, ``\widehat{(1,2,2,3)}=(1,2,3)``.

```julia
unique(Knots([1,2,2,3])) # Knots([1,2,3])
```

## B-spline space
Space of Piecewise polynominal.
```math
\mathcal{P}[p,k]
=..
```
where ``p\ge 0`` is called polynominal degree of space, and ``k`` is a knot vector.

```julia
p = 2
k = [1,3,5,6,8,9]
BSplineSpace(p,k)
𝒫(p,k) # same as above, for legibility
```

A B-spline space is said to be **proper** if its degree and knots satisfies following property:
```math
\begin{aligned}
k_{i}&<k_{i+p+1} & (1 \le i \le l-p-1)
\end{aligned}
```

```julia
isproper(𝒫(2,Knots([1,3,5,6,8,9]))) # true
isproper(𝒫(1,Knots([1,3,3,3,8,9]))) # false
```

The B-spline space is linear space, and if a B-spline space is proper, its dimension is calculated by:
```math
\dim(\mathcal{P}[p,k])=\sharp k -p -1
```

```julia
dim(𝒫(2,Knots([1,3,5,6,8,9]))) # 3
```

## Inclusive relationship between B-spline spaces
For proper B-spline spaces, the following relationship holds.
```math
\mathcal{P}[p,k]
\subseteq \mathcal{P}[p',k']
\Leftrightarrow (m=p'-p \ge 0 \ \text{and} \ k+m\widehat{k}\subseteq k')
```

```julia
P1 = 𝒫(2,Knots([1,3,5])
P2 = 𝒫(2,Knots([1,3,5,6,8,9])
P3 = 𝒫(3,Knots([1,1,3,3,5,5])
P1 ⊆ P2 # true
P1 ⊆ P3 # true
P2 ⊆ P3 # false
```


## B-spline basis function
B-spline basis function is defined by Cox–de Boor recursion formula.
```math
\begin{aligned}
{B}_{(i,p,k)}(t)
&=
\frac{t-k_{i}}{k_{i+p}-k_{i}}{B}_{(i,p-1,k)}(t)
+\frac{k_{i+p+1}-t}{k_{i+p+1}-k_{i+1}}{B}_{(i+1,p-1,k)}(t) \\
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}\le t< k_{i+1})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

You can choose the first terms in different ways.

```math
\begin{aligned}
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i} < t \le k_{i+1}\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

```math
\begin{aligned}
{B}_{(i,0,k)}(t)
&=
\begin{cases}
    &1\quad (k_{i}\le t<k_{i+1}<k_{l})\\
    &1\quad (k_{i}\le t\le k_{i+1}=k_{l})\\
    &0\quad (\text{otherwise})
\end{cases}
\end{aligned}
```

```julia
using Plots
gr()
...
```


## Support of B-spline basis function
If a B-spline space``\mathcal{P}[p,k]`` is proper, the support of its basis function is calculated as follows:
```math
\operatorname{supp}(B_{(i,p,k)})=[k_{i},k_{i+p+1}]
```


```julia
i=2
k=Knots([5,12,13,13,14])
p=2
P=𝒫(p,k)
BSplineSupport(P) # [5..13, 12..14]
BSplineSupport(i,P) # 12..14
```

## Derivative of B-spline basis function
```math
\begin{aligned}
\dot{B}_{(i,p,k)}(t)
&=\frac{d}{dt}B_{(i,p,k)}(t) \\
&=p\left(\frac{1}{k_{i+p}-k_{i}}B_{(i,p-1,k)}(t)-\frac{1}{k_{i+p+1}-k_{i+1}}B_{(i+1,p-1,k)}(t)\right)
\end{aligned}
```

Note that``\dot{B}_{(i,p,k)}\in\mathcal{P}[p-1,k]``.

## Partition of unity


```math
\begin{aligned}
\sum_{i}B_{(i,p,k)}(t) &= 1 & (k_{p+1}< t <k_{l-p}) \\
0 \le B_{(i,p,k)}(t) &\le 1
\end{aligned}
```




## Multi-dimentional B-spline
```math
B_{i^1,\dots,i^d}(t^1,\dots,t^d)
=B_{(i^1,p^1,k^1)}(t^1)\cdots B_{(i^d,p^d,k^d)}(t^d)
```


## B-spline manifold
```math
\bm{p}(t^1,\dots,t^d;\bm{a}_{i^1,\dots,i^d})
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
where ``\bm{a}_{i^1,\dots,i^d}`` is called **control points**

We will also write ``\bm{p}(t^1,\dots,t^d; \bm{a})``, ``\bm{p}(t^1,\dots,t^d)``, ``\bm{p}(t; \bm{a})`` or ``\bm{p}(t)`` for simplicity.

### B-spline curve

### B-spline surface

## Affine commutativity
If ``T`` is a affine transform ``T:\mathbb{R}^d\to\mathbb{R}^d``, then ..
```math
T(\bm{p}(t; \bm{a}))
=\bm{p}(t; T(a))
```

## Refinement
