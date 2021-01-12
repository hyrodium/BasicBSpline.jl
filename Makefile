preview:
	julia --project=docs -e 'using Pkg;Pkg.develop(PackageSpec(path=pwd()));Pkg.instantiate();include("docs/make.jl");'
	julia -e 'using LiveServer; serve(dir="docs/build")'
