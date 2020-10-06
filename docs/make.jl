using Documenter, BasicBSpline

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
