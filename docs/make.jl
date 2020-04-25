using Documenter, BasicBSpline

makedocs(;
    modules=[BasicBSpline],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/hyrodium/BasicBSpline.jl/blob/{commit}{path}#L{line}",
    sitename="BasicBSpline.jl",
    authors="hyrodium <hyrodium@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/hyrodium/BasicBSpline.jl",
)
