var documenterSearchIndex = {"docs":
[{"location":"detail/#Docstrings-1","page":"Docstrings","title":"Docstrings","text":"","category":"section"},{"location":"detail/#All-docstrings-1","page":"Docstrings","title":"All docstrings","text":"","category":"section"},{"location":"detail/#","page":"Docstrings","title":"Docstrings","text":"Modules = [BasicBSpline]","category":"page"},{"location":"detail/#BasicBSpline.BSplineManifold","page":"Docstrings","title":"BasicBSpline.BSplineManifold","text":"B-spline manifold for general polynomial degree\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.BSplineManifold-Tuple{AbstractBSplineManifold}","page":"Docstrings","title":"BasicBSpline.BSplineManifold","text":"convert AbstractBSplineManifold to BSplineManifold\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineSpace","page":"Docstrings","title":"BasicBSpline.BSplineSpace","text":"Construct B-spline space from given polynominal degree and knot vector.\n\nmathcalPpk\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.BSplineSpace-Tuple{AbstractBSplineSpace}","page":"Docstrings","title":"BasicBSpline.BSplineSpace","text":"convert AbstractBSplineSpace to BSplineSpace\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.FastBSplineManifold","page":"Docstrings","title":"BasicBSpline.FastBSplineManifold","text":"B-spline manifold for lower polynomial degree TODO: make the field bsplinespaces to be conposite type, not abstract type, for performance\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.FastBSplineManifold-Tuple{AbstractBSplineManifold}","page":"Docstrings","title":"BasicBSpline.FastBSplineManifold","text":"convert AbstractBSplineManifold to FastBSplineManifold\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.FastBSplineSpace","page":"Docstrings","title":"BasicBSpline.FastBSplineSpace","text":"B-spline space for lower polynomial degree\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.FastBSplineSpace-Tuple{AbstractBSplineSpace}","page":"Docstrings","title":"BasicBSpline.FastBSplineSpace","text":"convert AbstractBSplineSpace to FastBSplineSpace\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.Knots","page":"Docstrings","title":"BasicBSpline.Knots","text":"Construct knot vector from given array.\n\nk=(k_1dotsk_l)\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.FittingControlPoints-Tuple{Function,Array{T,1} where T<:AbstractBSplineSpace}","page":"Docstrings","title":"BasicBSpline.FittingControlPoints","text":"Approximate given function by linear combination of B-spline functions. This function returns its control points. TODO: currently, this function only supports for 1-dim and 2-dim B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis-Tuple{Array{BSplineSpace,1},Array{#s16,1} where #s16<:Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis","text":"Multi-dimensional B-spline basis function.\n\nB_i^1dotsi^d(t^1dotst^d)\n=B_(i^1p^1k^1)(t^1)cdots B_(i^dp^dk^d)(t^d)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis-Tuple{BSplineSpace,Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis","text":"B-spline basis function. Modified version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile tk_i+1k_l)\n    1quad (k_ile tle k_i+1=k_l)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis-Tuple{CartesianIndex,Array{BSplineSpace,1},Array{#s16,1} where #s16<:Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis","text":"Multi-dimensional B-spline basis function.\n\nB_i^1dotsi^d(t^1dotst^d)\n=B_(i^1p^1k^1)(t^1)cdots B_(i^dp^dk^d)(t^d)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis-Tuple{Integer,BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.bsplinebasis","text":"i-th B-spline basis function. Modified version.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis′-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.bsplinebasis′","text":"1st derivative of B-spline basis function. Modified version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis′-Tuple{Integer,FastBSplineSpace,Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis′","text":"TODO: faster.....\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis′₊₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.bsplinebasis′₊₀","text":"1st derivative of B-spline basis function. Right-sided limit version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis′₋₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.bsplinebasis′₋₀","text":"1st derivative of B-spline basis function. Left-sided limit version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis₊₀-Tuple{BSplineSpace,Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis₊₀","text":"B-spline basis function. Right-sided limit version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile t k_i+1)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis₊₀-Tuple{Integer,BSplineSpace,Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis₊₀","text":"i-th B-spline basis function. Right-sided limit version.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis₋₀-Tuple{BSplineSpace,Real}","page":"Docstrings","title":"BasicBSpline.bsplinebasis₋₀","text":"B-spline basis function. Left-sided limit version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_i tle k_i+1)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinebasis₋₀-Tuple{Integer,BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.bsplinebasis₋₀","text":"i-th B-spline basis function. Left-sided limit version.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.bsplinesupport-Tuple{Integer,AbstractBSplineSpace}","page":"Docstrings","title":"BasicBSpline.bsplinesupport","text":"Return support of i-th B-spline basis function.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.dim-Tuple{AbstractBSplineManifold}","page":"Docstrings","title":"BasicBSpline.dim","text":"Calculate the dimension of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.dim-Tuple{AbstractBSplineSpace}","page":"Docstrings","title":"BasicBSpline.dim","text":"Return dimention of a B-spline space.\n\ndim(mathcalPpk)\n=sharp k - p -1\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.mapping-Tuple{BSplineManifold,Array{#s16,1} where #s16<:Real}","page":"Docstrings","title":"BasicBSpline.mapping","text":"Calculate the mapping of B-spline manifold for given parameter.\n\nbmp(t^1dotst^d)\n=sum_i^1dotsi^dB_i^1dotsi^d(t^1dotst^d) bma_i^1dotsi^d\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.refinement-Tuple{BSplineManifold,Array{BSplineSpace,1}}","page":"Docstrings","title":"BasicBSpline.refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.refinement-Tuple{BSplineManifold}","page":"Docstrings","title":"BasicBSpline.refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.refinement-Tuple{FastBSplineManifold,Array{T,1} where T<:FastBSplineSpace}","page":"Docstrings","title":"BasicBSpline.refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.refinement-Tuple{FastBSplineManifold}","page":"Docstrings","title":"BasicBSpline.refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.𝒫-Tuple{Int64,Knots}","page":"Docstrings","title":"BasicBSpline.𝒫","text":"Retrun FastBSplineSpace if ≤ MAX_DEGREE, or BSplineSpace if not.\n\nmathcalPpk\n\n\n\n\n\n","category":"method"},{"location":"detail/#Base.:⊆-Tuple{AbstractBSplineSpace,AbstractBSplineSpace}","page":"Docstrings","title":"Base.:⊆","text":"Check inclusive relationship between B-spline spaces.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.GaussianQuadrature_1dim-Tuple{Function,Array{IntervalSets.Interval{:closed,:closed,T},1} where T<:Real,Any,Any}","page":"Docstrings","title":"BasicBSpline.GaussianQuadrature_1dim","text":"fast, but only for 1-dim\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.GaussianQuadrature_2dim-Tuple{Function,Array{IntervalSets.Interval{:closed,:closed,T},1} where T<:Real,Any,Any}","page":"Docstrings","title":"BasicBSpline.GaussianQuadrature_2dim","text":"fast, but only for 2-dim\n\n\n\n\n\n","category":"method"},{"location":"future/#Future-work-1","page":"Future work","title":"Future work","text":"","category":"section"},{"location":"future/#","page":"Future work","title":"Future work","text":"Add support for Drawing function\nThis may be in other repository.\nAdd support for NURBS\nThis may be in other repository.\nAdd example for IGA\n...","category":"page"},{"location":"contributing/#Contributing-1","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"The main contributer Hyrodium is not native English speaker. So, English corrections would be really helpful. Of course, other code improvement are welcomed!","category":"page"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"Feel free to open issue, and pull request!","category":"page"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"Issues\nPull requests","category":"page"},{"location":"math/#Mathematical-properties-of-B-spline-1","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"","category":"section"},{"location":"math/#Introduction-1","page":"Mathematical properties of B-spline","title":"Introduction","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B-spline is a mathematical object, and it has a lot of application(e.g. NURBS, IGA).","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"In this page, we'll explain the mathematical definition and property of B-spline with Julia code.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Before running the following code, do not forget importing the package:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using BasicBSpline","category":"page"},{"location":"math/#Notice-1","page":"Mathematical properties of B-spline","title":"Notice","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"A book \"Geometric Modeling with Splines\" by Elaine Cohen, Richard F. Riesenfeld, Gershon Elber is really recommended.\nBut, some of notations in this page are my original, but these are well-considered results.","category":"page"},{"location":"math/#Knot-vector-1","page":"Mathematical properties of B-spline","title":"Knot vector","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tip: Def.  Knot vector\nA finite sequencek = (k_1 dots k_l)is called knot vector if the sequence is broad monotonic increase, i.e. k_i le k_i+1.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"[fig]","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k = Knots([1,2,3])\nk = Knots(1:3)\nk = Knots(1,2,3)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We denotes a number of knots by sharp symbol like this:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"sharp k = sharp(k_1 dots k_l) =l","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k = Knots([4,5,6])\n♯(k) # 3\nlength(k) # 3","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We introduce additional operator + and product operator cdot","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nk^(1)+k^(2)\n=(k^(1)_1 dots k^(1)_l) + (k^(2)_1 dots k^(2)_l) \n=(textsort of union of   k^(1)  textand   k^(2) text) \nmcdot k=underbracek+cdots+k_m\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For example, (123)+(245)=(12235), 2cdot (23)=(2233).","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Knots([1,2,3]) + Knots([2,4,5]) # Knots([1,2,2,3,5])\n2 * Knots([2,3]) # Knots([2,2,3,3])","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Deleting duplicates operator","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nwidehatk\n=(textremove duplicates of   k) \nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For example, widehat(1223)=(123).","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"unique(Knots([1,2,2,3])) # Knots([1,2,3])","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For Given knot vector k, the following function mathfrakn_kmathbbRtomathbbZ represents the counts of ..","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"mathfrakn_k(t) = sharpi mid k_i=t ","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For example, if k=(1223), then mathfrakn_k(03)=0, mathfrakn_k(1)=1, mathfrakn_k(2)=2.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k = Knots([1,2,2,3])\n𝔫(k,0.3) # 0\n𝔫(k,1.0) # 1\n𝔫(k,2.0) # 2","category":"page"},{"location":"math/#B-spline-space-1","page":"Mathematical properties of B-spline","title":"B-spline space","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Before defining B-spline space, we'll define polynomial space with degree p.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tip: Def.  Polynomial space\nPolynomial space with degree p.mathcalPp\n=leftfmathbbRtomathbbR  tmapsto a_0+a_1t^1+cdots+a_pt^p   left  \n    a_iin mathbbR\n    right\nrightThis space mathcalPp is a p+1-dimensional linear space.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Note that tmapsto t^i_0 le i le p is a basis of mathcalPp, and also the set of Bernstein polynomial B_(ip)_i is a basis of mathcalPp.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(ip)(t)\n=binompi-1t^i-1(1-t)^p-i+1\n(i=1 dots p+1)\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Where binompi-1 is a binomial coefficient.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tip: Def.  B-spline space\nFor given polynomial degree ple 0 and knot vector k=(k_1dotsk_l), B-spline space mathcalPpk is defined as follows:mathcalPpk\n=leftfmathbbRtomathbbR   left  \n    begingathered\n        f((-infty k_1))=f(k_l+infty))=left0right \n        exists tildefinmathcalPp f_k_i k_i+1) = tildef_k_i k_i+1)  \n        forall t in mathbbR exists delta  0 f_(t-deltat+delta)in C^p-mathfrakn_k(t)\n    endgathered right\nright","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Note that a element of the space mathcalPpk is piecewise polynomial.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"[fig]","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"p = 2\nk = Knots([1,3,5,6,8,9])\nBSplineSpace(p,k)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"A B-spline space is said to be proper if its degree and knots satisfies following property:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nk_ik_i+p+1  (1 le i le l-p-1)\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"isproper(𝒫(2,Knots([1,3,5,6,8,9]))) # true\nisproper(𝒫(1,Knots([1,3,3,3,8,9]))) # false","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"The B-spline space is linear space, and if a B-spline space is proper, its dimension is calculated by:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"dim(mathcalPpk)=sharp k -p -1","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"dim(𝒫(2,Knots([1,3,5,6,8,9]))) # 3","category":"page"},{"location":"math/#B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tip: Def.  B-spline space\nB-spline basis function is defined by Cox–de Boor recursion formula.beginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile t k_i+1)\n    0quad (textotherwise)\nendcases\nendalignedIf the denominator is ..","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Basis of B-spline space\nThe set of functions B_(ipk)_i is a basis of B-spline space mathcalPpk.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np = 2\nk = Knots(1:8)\nP = BSplineSpace(p,k)\nplot([t->bsplinebasis₊₀(i,P,t) for i in 1:dim(P)], 1, 8, ylims=(0,1.05))","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"You can choose the first terms in different ways.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(i0k)(t)\n=\nbegincases\n    1quad (k_i  t le k_i+1) \n    0quad (textotherwise)\nendcases\nendaligned","category":"page"},{"location":"math/#Support-of-B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"Support of B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Support of B-spline basis function\nIf a B-spline spacemathcalPpk is proper, the support of its basis function is calculated as follows:operatornamesupp(B_(ipk))=k_ik_i+p+1","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"[fig]","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"i = 2\nk = Knots([5,12,13,13,14])\np = 2\nP = 𝒫(p,k)\nbsplinesupport(P) # [5..13, 12..14]\nbsplinesupport(i,P) # 12..14","category":"page"},{"location":"math/#Derivative-of-B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"Derivative of B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Derivative of B-spline basis function\nThe derivative of B-spline basis function can be expressed as follows:beginaligned\ndotB_(ipk)(t)\n=fracddtB_(ipk)(t) \n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\nendalignedNote that dotB_(ipk)inmathcalPp-1k.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np = 2\nk = Knots(1:8)\nP = BSplineSpace(p,k)\nplot([t->bsplinebasis′₊₀(i,P,t) for i in 1:dim(P)], 1, 8, ylims=(0,1.05))","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#Partition-of-unity-1","page":"Mathematical properties of B-spline","title":"Partition of unity","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Partition of unity\nbeginaligned\nsum_iB_(ipk)(t) = 1  (k_p+1 le t  k_l-p) \n0 le B_(ipk)(t) le 1\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np = 2\nk = Knots(1:8)\nP = BSplineSpace(p,k)\nplot(t->sum(bsplinebasis₊₀(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"To satisfy the partition of unity on whole interval 18, sometimes more knots will be inserted to the endpoints of the interval.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np = 2\nk = Knots(1:8) + p * Knots([1,8])\nP = BSplineSpace(p,k)\nplot(t->sum(bsplinebasis₊₀(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"But, the sum sum_i B_(ipk)(t) is not equal to 1 if t=8. Therefore, to satisfy partition of unity on closed interval k_p+1 k_l-p, the definition of first terms of B-spline basis functions are sometimes replaced:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(i0k)(t)\n=\nbegincases\n    1quad (k_i le tk_i+1k_l)\n    1quad (k_i le t le k_i+1=k_l)\n    0quad (textotherwise)\nendcases\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np = 2\nk = Knots(1:8) + p * Knots([1,8])\nP = BSplineSpace(p,k)\nplot(t->sum(bsplinebasis(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#Inclusion-relation-between-B-spline-spaces-1","page":"Mathematical properties of B-spline","title":"Inclusion relation between B-spline spaces","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Support of B-spline basis function\nFor proper B-spline spaces, the following relationship holds.mathcalPpk\nsubseteq mathcalPpk\nLeftrightarrow (m=p-p ge 0  textand  k+mwidehatksubseteq k)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(as linear subspace..)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"P1 = 𝒫(1,Knots([1,3,5,8]))\nP2 = 𝒫(1,Knots([1,3,5,6,8,9]))\nP3 = 𝒫(2,Knots([1,1,3,3,5,5,8,8]))\nP1 ⊆ P2 # true\nP1 ⊆ P3 # true\nP2 ⊆ P3 # false\nP2 ⊈ P3 # true","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Here are plots of the B-spline basis functions of the spaces P1, P2, P3.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\nP1 = 𝒫(1,Knots([1,3,5,8]))\nP2 = 𝒫(1,Knots([1,3,5,6,8,9]))\nP3 = 𝒫(2,Knots([1,1,3,3,5,5,8,8]))\nplot(\n    plot([t->bsplinebasis₊₀(i,P1,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),\n    plot([t->bsplinebasis₊₀(i,P2,t) for i in 1:dim(P2)], 1, 9, ylims=(0,1.05), legend=false),\n    plot([t->bsplinebasis₊₀(i,P3,t) for i in 1:dim(P3)], 1, 9, ylims=(0,1.05), legend=false),\n    layout=(3,1),\n    link=:x\n)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"This means, there exists a n times n matrix A which holds:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(ipk)\n=sum_jA_ij B_(jpk) \nn=dim(mathcalPpk) \nn=dim(mathcalPpk)\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"You can calculate the change of basis matrix A with changebasis.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"A12 = changebasis(P1,P2)\nA13 = changebasis(P1,P3)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\nplot(\n    plot([t->bsplinebasis₊₀(i,P1,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),\n    plot([t->sum(A12[i,j]*bsplinebasis₊₀(j,P2,t) for j in 1:dim(P2)) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),\n    plot([t->sum(A13[i,j]*bsplinebasis₊₀(j,P3,t) for j in 1:dim(P3)) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),\n    layout=(3,1),\n    link=:x\n)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"(Image: )","category":"page"},{"location":"math/#Multi-dimensional-B-spline-1","page":"Mathematical properties of B-spline","title":"Multi-dimensional B-spline","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tensor product","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B_i^1dotsi^d(t^1dotst^d)\n=B_(i^1p^1k^1)(t^1)cdots B_(i^dp^dk^d)(t^d)","category":"page"},{"location":"math/#B-spline-manifold-1","page":"Mathematical properties of B-spline","title":"B-spline manifold","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B-spline manifold is a parametric representation of a shape.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"tip: Def.  B-spline manifold\nFor given d-dimensional B-spline basis functions B_i^1dotsi^d and given points bma_i^1dotsi^d in mathbbR^hatd, B-spline manifold is defined by following equality:bmp(t^1dotst^dbma_i^1dotsi^d)\n=sum_i^1dotsi^dB_i^1dotsi^d(t^1dotst^d) bma_i^1dotsi^dWhere bma_i^1dotsi^d are called control points.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We will also write bmp(t^1dotst^d bma), bmp(t^1dotst^d), bmp(t bma) or bmp(t) for simplicity.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"P1 = 𝒫(1,Knots([0,0,1,1]))\nP2 = 𝒫(1,Knots([1,1,2,3,3]))\nn1 = dim(P1) # 2\nn2 = dim(P2) # 3\n𝒂 = [[i, j] for i in 1:n1, j in 1:n2]  # n1 × n2 array of d̂ array.\nM = BSplineManifold([P1, P2], 𝒂)","category":"page"},{"location":"math/#B-spline-curve-1","page":"Mathematical properties of B-spline","title":"B-spline curve","text":"","category":"section"},{"location":"math/#B-spline-surface-1","page":"Mathematical properties of B-spline","title":"B-spline surface","text":"","category":"section"},{"location":"math/#Affine-commutativity-1","page":"Mathematical properties of B-spline","title":"Affine commutativity","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"info: Thm.  Affine commutativity\nIf T is a affine transform mathbbR^dtomathbbR^d, then the following equality holds.T(bmp(t bma))\n=bmp(t T(bma))","category":"page"},{"location":"math/#Refinement-1","page":"Mathematical properties of B-spline","title":"Refinement","text":"","category":"section"},{"location":"math/#h-refinemnet-1","page":"Mathematical properties of B-spline","title":"h-refinemnet","text":"","category":"section"},{"location":"math/#p-refinemnet-1","page":"Mathematical properties of B-spline","title":"p-refinemnet","text":"","category":"section"},{"location":"#BasicBSpline.jl-1","page":"Home","title":"BasicBSpline.jl","text":"","category":"section"},{"location":"#Summary-1","page":"Home","title":"Summary","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package provides basic (mathematical) operations for B-spline.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The package Interpolations.jl says:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Currently this package's support is best for B-splines and also supports irregular grids. However, the API has been designed with intent to support more options. Pull-requests are more than welcome! It should be noted that the API may continue to evolve over time.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"As mentioned before, this package treats mathematical aspect of B-spline, so the difference between these packages is a main purpose.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"If you are interested in Interpolations, Interpolations.jl would be helpful.\nIf you would like to deal with raw B-spline functions, this package would be the best for you.  For example:\nB-spline curve\nB-spline surface\nNURBS (Non-Uniform Rational B-Spline)\nIGA (Isogeometric Analysis)","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"add https://github.com/hyrodium/BasicBSpline.jl","category":"page"},{"location":"#Example-Images-1","page":"Home","title":"Example Images","text":"","category":"section"},{"location":"#Example-of-B-spline-function-1","page":"Home","title":"Example of B-spline function","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Note that this package do not support image export. (Image: )","category":"page"}]
}
