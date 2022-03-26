using Documenter, BasicBSpline, BasicBSplineExporter, Plots

gr()
plotlyjs()

# Setup for doctests in docstrings
DocMeta.setdocmeta!(BasicBSpline, :DocTestSetup, :(using LinearAlgebra, BasicBSpline))

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
        canonical = "https://hyrodium.github.io/BasicBSpline.jl/stable/",
        assets = ["assets/favicon.ico", "assets/custom.css"],
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical properties of B-spline" => "math.md",
        # "Differentiation" => [
        #     "BSplineDerivativeSpace" => "bsplinederivativespace.md",
        #     "ForwardDiff" => "forwarddiff.md",
        #     "ChainRules" => "chainrules.md"
        # ],
        "Visualization" => [
            "BasicBSplineExporter.jl" => "basicbsplineexporter.md",
            "PlotlyJS.jl" => "plotlyjs.md"
        ],
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
