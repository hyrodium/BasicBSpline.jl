using Documenter, BasicBSpline

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
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Mathematical properties of B-spline" => "math.md",
        "Benchmark" => "benchmark.md",
        "Example 1 - Poisson's equation" => "poisson.md",
        "Example 2 - Linear elasticity" => "elasticity.md",
        "Docstrings" => "detail.md",
        "Future work" => "future.md",
        "Contributing" => "contributing.md",
    ],
    repo = "https://github.com/hyrodium/BasicBSpline.jl/blob/{commit}{path}#L{line}",
    sitename = "BasicBSpline.jl",
    authors = "hyrodium <hyrodium@gmail.com>",
    assets = String[],
)

deploydocs(; repo = "github.com/hyrodium/BasicBSpline.jl")
