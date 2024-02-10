using Documenter
using BasicBSpline
using BasicBSplineExporter
using BasicBSplineFitting
using InteractiveUtils
using Plots
using Random

gr()
plotly()
Random.seed!(42)

# Setup for doctests in docstrings
DocMeta.setdocmeta!(BasicBSpline, :DocTestSetup, :(using LinearAlgebra, BasicBSpline, StaticArrays, BasicBSplineFitting))
DocMeta.setdocmeta!(BasicBSplineFitting, :DocTestSetup, :(using LinearAlgebra, BasicBSpline, StaticArrays, BasicBSplineFitting))

function generate_indexmd_from_readmemd()
    path_readme = "README.md"
    path_index = "docs/src/index.md"

    # open README.md
    text_readme = read(path_readme, String)

    # generate text for index.md
    text_index = text_readme
    text_index = replace(text_index, "![](docs/src/img" => "![](img")
    text_index = replace(text_index, r"\$\$((.|\n)*?)\$\$" => s"```math\1```")
    text_index = replace(text_index, "(https://hyrodium.github.io/BasicBSpline.jl/dev/math/bsplinebasis/#Differentiability-and-knot-duplications)" => "(@ref differentiability-and-knot-duplications)")
    text_index = replace(text_index, r"https://github.com/hyrodium/BasicBSpline\.jl/assets/.*" => "![](math/differentiability.mp4)")
    text_index = """
    ```@meta
    EditURL = "https://github.com/hyrodium/BasicBSpline.jl/blob/main/README.md"
    ```

    """ * text_index

    # save index.md
    open(path_index, "w") do f
        write(f,text_index)
    end
end

generate_indexmd_from_readmemd()

makedocs(;
    modules = [BasicBSpline, BasicBSplineFitting],
    format = Documenter.HTML(
        ansicolor=true,
        canonical = "https://hyrodium.github.io/BasicBSpline.jl",
        assets = ["assets/favicon.ico", "assets/custom.css"],
        edit_link="main",
        repolink="https://github.com/hyrodium/BasicBSpline.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical properties of B-spline" => [
            "Introduction" => "math/introduction.md",
            "Knot vector" => "math/knotvector.md",
            "B-spline space" => "math/bsplinespace.md",
            "B-spline basis function" => "math/bsplinebasis.md",
            "B-spline manifold" => "math/bsplinemanifold.md",
            "Rational B-spline manifold" => "math/rationalbsplinemanifold.md",
            "Inclusive relationship" => "math/inclusive.md",
            "Refinement" => "math/refinement.md",
            "Derivative" => "math/derivative.md",
            "Fitting" => "math/fitting.md",
        ],
        # "Differentiation" => [
        #     "BSplineDerivativeSpace" => "bsplinederivativespace.md",
        #     "ForwardDiff" => "forwarddiff.md",
        #     "ChainRules" => "chainrules.md"
        # ],
        "Visualization" => [
            "PlotlyJS.jl" => "visualization/plotlyjs.md",
            "BasicBSplineExporter.jl" => "visualization/basicbsplineexporter.md",
        ],
        "Interpolations" => "interpolations.md",
        "API" => "api.md",
    ],
    sitename = "BasicBSpline.jl",
    authors = "hyrodium <hyrodium@gmail.com>",
    warnonly = true,
)

deploydocs(
    repo = "github.com/hyrodium/BasicBSpline.jl",
    push_preview = true,
    devbranch="main",
)
