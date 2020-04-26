var documenterSearchIndex = {"docs":
[{"location":"detail/#Docstrings-1","page":"Docstrings","title":"Docstrings","text":"","category":"section"},{"location":"detail/#All-docstrings-1","page":"Docstrings","title":"All docstrings","text":"","category":"section"},{"location":"detail/#","page":"Docstrings","title":"Docstrings","text":"Modules = [BasicBSpline]","category":"page"},{"location":"detail/#BasicBSpline.BSplineSpace","page":"Docstrings","title":"BasicBSpline.BSplineSpace","text":"Construct B-spline space from given polynominal degree and knot vector.\n\nmathcalPpk\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.Knots","page":"Docstrings","title":"BasicBSpline.Knots","text":"Construct knot vector from given array.\n\nk=(k_1dotsk_l)\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.𝒫","page":"Docstrings","title":"BasicBSpline.𝒫","text":"Same as BSplineSpace.\n\nmathcalPpk\n\n\n\n\n\n","category":"type"},{"location":"detail/#BasicBSpline.BSplineBasis-Tuple{Array{BSplineSpace,1},Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis","text":"Multi-dimentional B-spline basis function.\n\nB_i^1dotsi^d(t^1dotst^d)\n=B_(i^1p^1k^1)(t^1)cdots B_(i^dp^dk^d)(t^d)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis","text":"B-spline basis function. Modified version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile tk_i+1k_l)\n    1quad (k_ile tle k_i+1=k_l)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis-Tuple{Int64,BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis","text":"i-th B-spline basis function. Modified version.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis′-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis′","text":"1st derivative of B-spline basis function. Modified version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis′₊₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis′₊₀","text":"1st derivative of B-spline basis function. Right-sided limit version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis′₋₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis′₋₀","text":"1st derivative of B-spline basis function. Left-sided limit version.\n\ndotB_(ipk)(t)\n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis₊₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis₊₀","text":"B-spline basis function. Right-sided limit version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile t k_i+1)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineBasis₋₀-Tuple{BSplineSpace,Any}","page":"Docstrings","title":"BasicBSpline.BSplineBasis₋₀","text":"B-spline basis function. Left-sided limit version.\n\nbeginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_i tle k_i+1)\n    0quad (textotherwise)\nendcases\nendaligned\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.BSplineSupport-Tuple{Int64,BSplineSpace}","page":"Docstrings","title":"BasicBSpline.BSplineSupport","text":"Return support of i-th B-spline basis function.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.Mapping-Tuple{BSplineManifold,Array{Float64,1}}","page":"Docstrings","title":"BasicBSpline.Mapping","text":"Calculate the mapping of B-spline manifold for given parameter.\n\nbmp(t^1dotst^d)\n=sum_i^1dotsi^dB_i^1dotsi^d(t^1dotst^d) bma_i^1dotsi^d\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.Refinement-Tuple{BSplineManifold,Array{BSplineSpace,1}}","page":"Docstrings","title":"BasicBSpline.Refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.Refinement-Tuple{BSplineManifold}","page":"Docstrings","title":"BasicBSpline.Refinement","text":"Refinement of B-spline manifold.\n\n\n\n\n\n","category":"method"},{"location":"detail/#BasicBSpline.dim-Tuple{BSplineSpace}","page":"Docstrings","title":"BasicBSpline.dim","text":"Return dimention of a B-spline space.\n\ndim(mathcalPpk)\n=sharp k - p -1\n\n\n\n\n\n","category":"method"},{"location":"detail/#Base.:⊆-Tuple{BSplineSpace,BSplineSpace}","page":"Docstrings","title":"Base.:⊆","text":"Check inclusive relationship between B-spline spaces.\n\n\n\n\n\n","category":"method"},{"location":"future/#Future-work-1","page":"Future work","title":"Future work","text":"","category":"section"},{"location":"future/#","page":"Future work","title":"Future work","text":"Add support for Drawing function\nThis may be in other repository.\nAdd support for NURBS\nThis may be in other repository.\nAdd example for IGA\n...","category":"page"},{"location":"contributing/#Contributing-1","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"The main contributer Hyrodium is not native English speaker. So, English corrections would be really helpful. Of course, other code improvement are welcomed!","category":"page"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"Feel free to open issue, and pull request!","category":"page"},{"location":"contributing/#","page":"Contributing","title":"Contributing","text":"Issues\nPull requests","category":"page"},{"location":"math/#Mathematical-properties-of-B-spline-1","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"","category":"section"},{"location":"math/#Introduction-1","page":"Mathematical properties of B-spline","title":"Introduction","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B-spline is a mathematical object, and it has a lot of application(e.g. NURBS, IGA).","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"In this page, we'll explain the mathematical definition and property of B-spline with Julia code.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Before running the following code, do not forget importing the package:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using BasicBSpline","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Notice: Some of notations of B-spline are my original, but these are well-considered results.","category":"page"},{"location":"math/#Knot-vector-1","page":"Mathematical properties of B-spline","title":"Knot vector","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"A finite sequence","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k=(k_1 dots k_l)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"is called knot vector if the sequence is broad monotonic increase, i.e. (k_i le k_i+1)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k=Knots([1,2,3])","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We denotes a number of knots by sharp symbol like this:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"sharp k=sharp(k_1 dots k_l) =l","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"♯(k) # 3\nlength(k) # 3","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We introduce additional operator + and product operator cdot","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nk^(1)+k^(2)\n=(k^(1)_1 dots k^(1)_l) + (k^(2)_1 dots k^(2)_l) \n=(textsort of union of   k^(1)  textand   k^(2) text) \nnk=underbracek+cdots+k_n\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For example, (123)+(245)=(12235), 2cdot (23)=(2233)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"k = Knots([1,2,3]) + Knots([2,4,5]) # Knots([1,2,2,3,5])\n2Knots([2,3]) # Knots([2,2,3,3])","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Unique operator","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nwidehatk\n=(textremove duplicates of   k) \nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For example, widehat(1223)=(123).","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"unique(Knots([1,2,2,3])) # Knots([1,2,3])","category":"page"},{"location":"math/#B-spline-space-1","page":"Mathematical properties of B-spline","title":"B-spline space","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Space of Piecewise polynominal.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"mathcalPpk\n=","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"where pge 0 is called polynominal degree of space, and k is a knot vector.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"p = 2\nk = [1,3,5,6,8,9]\nBSplineSpace(p,k)\n𝒫(p,k) # same as above, for legibility","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"A B-spline space is said to be proper if its degree and knots satisfies following property:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nk_ik_i+p+1  (1 le i le l-p-1)\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"isproper(𝒫(2,Knots([1,3,5,6,8,9]))) # true\nisproper(𝒫(1,Knots([1,3,3,3,8,9]))) # false","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"The B-spline space is linear space, and if a B-spline space is proper, its dimension is calculated by:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"dim(mathcalPpk)=sharp k -p -1","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"dim(𝒫(2,Knots([1,3,5,6,8,9]))) # 3","category":"page"},{"location":"math/#Inclusive-relationship-between-B-spline-spaces-1","page":"Mathematical properties of B-spline","title":"Inclusive relationship between B-spline spaces","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"For proper B-spline spaces, the following relationship holds.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"mathcalPpk\nsubseteq mathcalPpk\nLeftrightarrow (m=p-p ge 0  textand  k+mwidehatksubseteq k)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"P1 = 𝒫(2,Knots([1,3,5])\nP2 = 𝒫(2,Knots([1,3,5,6,8,9])\nP3 = 𝒫(3,Knots([1,1,3,3,5,5])\nP1 ⊆ P2 # true\nP1 ⊆ P3 # true\nP2 ⊆ P3 # false","category":"page"},{"location":"math/#B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B-spline basis function is defined by Cox–de Boor recursion formula.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(ipk)(t)\n=\nfract-k_ik_i+p-k_iB_(ip-1k)(t)\n+frack_i+p+1-tk_i+p+1-k_i+1B_(i+1p-1k)(t) \nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile t k_i+1)\n    0quad (textotherwise)\nendcases\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"The set of functions B_(ipk)_i is a basis of B-spline space mathcalPpk.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np=2\nk=Knots(1:8)\nP=BSplineSpace(p,k)\nplot([t->BSplineBasis(i,P,t) for i in 1:dim(P)], 1,8)","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"You can choose the first terms in different ways.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(i0k)(t)\n=\nbegincases\n    1quad (k_i  t le k_i+1\n    0quad (textotherwise)\nendcases\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nB_(i0k)(t)\n=\nbegincases\n    1quad (k_ile tk_i+1k_l)\n    1quad (k_ile tle k_i+1=k_l)\n    0quad (textotherwise)\nendcases\nendaligned","category":"page"},{"location":"math/#Support-of-B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"Support of B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"If a B-spline spacemathcalPpk is proper, the support of its basis function is calculated as follows:","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"operatornamesupp(B_(ipk))=k_ik_i+p+1","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"i=2\nk=Knots([5,12,13,13,14])\np=2\nP=𝒫(p,k)\nBSplineSupport(P) # [5..13, 12..14]\nBSplineSupport(i,P) # 12..14","category":"page"},{"location":"math/#Derivative-of-B-spline-basis-function-1","page":"Mathematical properties of B-spline","title":"Derivative of B-spline basis function","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\ndotB_(ipk)(t)\n=fracddtB_(ipk)(t) \n=pleft(frac1k_i+p-k_iB_(ip-1k)(t)-frac1k_i+p+1-k_i+1B_(i+1p-1k)(t)right)\nendaligned","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"Note thatdotB_(ipk)inmathcalPp-1k.","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"using Plots\ngr()\np=2\nk=Knots(1:8)\nP=BSplineSpace(p,k)\nplot([t->BSplineBasis(i,P,t) for i in 1:dim(P)], 1,8)","category":"page"},{"location":"math/#Partition-of-unity-1","page":"Mathematical properties of B-spline","title":"Partition of unity","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"beginaligned\nsum_iB_(ipk)(t) = 1  (k_p+1 t k_l-p) \n0 le B_(ipk)(t) le 1\nendaligned","category":"page"},{"location":"math/#Multi-dimentional-B-spline-1","page":"Mathematical properties of B-spline","title":"Multi-dimentional B-spline","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"B_i^1dotsi^d(t^1dotst^d)\n=B_(i^1p^1k^1)(t^1)cdots B_(i^dp^dk^d)(t^d)","category":"page"},{"location":"math/#B-spline-manifold-1","page":"Mathematical properties of B-spline","title":"B-spline manifold","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"bmp(t^1dotst^dbma_i^1dotsi^d)\n=sum_i^1dotsi^dB_i^1dotsi^d(t^1dotst^d) bma_i^1dotsi^d","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"where bma_i^1dotsi^d is called control points","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"We will also write bmp(t^1dotst^d bma), bmp(t^1dotst^d), bmp(t bma) or bmp(t) for simplicity.","category":"page"},{"location":"math/#B-spline-curve-1","page":"Mathematical properties of B-spline","title":"B-spline curve","text":"","category":"section"},{"location":"math/#B-spline-surface-1","page":"Mathematical properties of B-spline","title":"B-spline surface","text":"","category":"section"},{"location":"math/#Affine-commutativity-1","page":"Mathematical properties of B-spline","title":"Affine commutativity","text":"","category":"section"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"If T is a affine transform TmathbbR^dtomathbbR^d, then ..","category":"page"},{"location":"math/#","page":"Mathematical properties of B-spline","title":"Mathematical properties of B-spline","text":"T(bmp(t bma))\n=bmp(t T(a))","category":"page"},{"location":"math/#Refinement-1","page":"Mathematical properties of B-spline","title":"Refinement","text":"","category":"section"},{"location":"#BasicBSpline.jl-1","page":"Home","title":"BasicBSpline.jl","text":"","category":"section"},{"location":"#Summary-1","page":"Home","title":"Summary","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package provides basic (mathematical) operations for B-spline.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"The package Interpolations.jl says:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Currently this package's support is best for B-splines and also supports irregular grids. However, the API has been designed with intent to support more options. Pull-requests are more than welcome! It should be noted that the API may continue to evolve over time.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"As mentioned before, this package treats mathematical aspect of B-spline, so the difference between these packages is a main purpose.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"If you are interested in Interpolations, Interpolations.jl would be helpful.\nIf you would like to deal with raw B-spline functions, this package would be the best for you.  For example:\nB-spline curve\nB-spline surface\nNURBS (Non-Uniform Rational B-Spline)\nIGA (Isogeometric Analysis)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"add https://github.com/hyrodium/BasicBSpline.jl","category":"page"}]
}
