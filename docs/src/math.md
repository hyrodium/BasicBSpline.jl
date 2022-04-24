# Mathematical properties of B-spline

## Introduction
[B-spline](https://en.wikipedia.org/wiki/B-spline) is a mathematical object, and it has a lot of application. (e.g. [NURBS](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline), [IGA](https://en.wikipedia.org/wiki/Isogeometric_analysis))

In this page, we'll explain the mathematical definition and property of B-spline with Julia code.
Before running the code in the following section, you need to import packages:
```@example
using BasicBSpline
using Plots; plotly()
```

## Notice
* A book ["Geometric Modeling with Splines"](https://www.routledge.com/p/book/9780367447243) by Elaine Cohen, Richard F. Riesenfeld, Gershon Elber is really recommended.
* **Some of notations in this page are our original**, but these are well-considered results.
