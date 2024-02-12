# Run this script after running make.jl

dir_build = joinpath(@__DIR__, "build")
dir_img = joinpath(@__DIR__, "src", "img")

for filename in readdir(dir_build)
    path_src = joinpath(dir_build, filename)
    if startswith(filename, "readme-") && endswith(filename, ".html")
        println(path_src)
        path_png = joinpath(dir_img, chopsuffix(chopprefix(filename, "readme-"), ".html") * ".png")
        cmd = `google-chrome-stable --headless --disable-gpu --screenshot=$path_png $path_src`
        run(cmd)
    end

    if startswith(filename, "readme-") && endswith(filename, ".png")
        println(path_src)
        path_png = joinpath(dir_img, chopprefix(filename, "readme-"))
        cp(path_src, path_png; force=true)
    end
end
