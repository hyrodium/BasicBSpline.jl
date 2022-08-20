using Documenter
using BasicBSpline
using BasicBSplineExporter
using Plots
using Random

gr()
plotly()
Random.seed!(42)

# Setup for doctests in docstrings
DocMeta.setdocmeta!(BasicBSpline, :DocTestSetup, :(using LinearAlgebra, BasicBSpline, StaticArrays))

function generate_indexmd_from_readmemd()
    path_readme = "README.md"
    path_index = "docs/src/index.md"

    # open README.md
    f = open(path_readme)
    text_readme = read(f,String)
    close(f)

    # generate text for index.md
    text_index = replace(text_readme,"![](docs/src/img" => "![](img")

    # save index.md
    open(path_index, "w") do f
        write(f,text_index)
    end
end

generate_indexmd_from_readmemd()

makedocs(;
    modules = [BasicBSpline],
    format = Documenter.HTML(
        ansicolor=true,
        canonical = "https://hyrodium.github.io/BasicBSpline.jl/stable/",
        assets = ["assets/favicon.ico", "assets/custom.css"],
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical properties of B-spline" => [
            "Introduction" => "math.md",
            "Knot vector" => "math-knotvector.md",
            "B-spline space" => "math-bsplinespace.md",
            "B-spline basis function" => "math-bsplinebasis.md",
            "B-spline manifold" => "math-bsplinemanifold.md",
            "Derivative" => "math-derivative.md",
            "Inclusive relationship" => "math-inclusive.md",
            "Refinement" => "math-refinement.md",
            "Fitting" => "math-fitting.md",
            "Rational B-spline manifold" => "math-rationalbsplinemanifold.md",
        ],
        # "Differentiation" => [
        #     "BSplineDerivativeSpace" => "bsplinederivativespace.md",
        #     "ForwardDiff" => "forwarddiff.md",
        #     "ChainRules" => "chainrules.md"
        # ],
        "Visualization" => [
            "Plots.jl" => "plots.md",
            "PlotlyJS.jl" => "plotlyjs.md",
            "BasicBSplineExporter.jl" => "basicbsplineexporter.md",
        ],
        "Interpolations" => "interpolations.md",
        "Private API" => "internal.md",
        "Contributing" => "contributing.md",
    ],
    repo = "https://github.com/hyrodium/BasicBSpline.jl/blob/{commit}{path}#L{line}",
    sitename = "BasicBSpline.jl",
    authors = "hyrodium <hyrodium@gmail.com>",
)

deploydocs(
    repo = "github.com/hyrodium/BasicBSpline.jl",
    push_preview = true
)
