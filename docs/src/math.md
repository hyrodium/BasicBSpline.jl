# Mathematical properties of B-spline

## Introduction
[B-spline](https://en.wikipedia.org/wiki/B-spline) is a mathematical object, and it has a lot of application. (e.g. [NURBS](https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline), [IGA](https://en.wikipedia.org/wiki/Isogeometric_analysis))

In this page, we'll explain the mathematical definitions and properties of B-spline with Julia code.
Before running the code in the following section, you need to import packages:
```@example
using BasicBSpline
using Plots; plotly()
```

## Notice
Some of notations in this page are our original, but these are well-considered results.

## Learning B-spline
Most of this documentation around B-spline is self-contained.
If you want to learn more, the following resources are recommended.

* ["Geometric Modeling with Splines"](https://www.routledge.com/p/book/9780367447243)
* [Spline Functions: Basic Theory](https://www.cambridge.org/core/books/spline-functions-basic-theory/843475201223F90091FFBDDCBF210BFB)

日本語の文献では以下がおすすめです。
* [NURBS多様体による形状表現](https://hyrodium.github.io/ja/pdf/#NURBS%E5%A4%9A%E6%A7%98%E4%BD%93%E3%81%AB%E3%82%88%E3%82%8B%E5%BD%A2%E7%8A%B6%E8%A1%A8%E7%8F%BE)
* [BasicBSpline.jlを作ったので宣伝です！](https://zenn.dev/hyrodium/articles/5fb08f98d4a918)
