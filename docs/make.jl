using Documenter, Liga

makedocs(
	assets = ["assets/favicon.ico"],
	sitename = "Liga- Library for Geometric Algebra",
	pages =["Overview" => "index.md",
			  "Getting Started" => "gettingstarted.md",
			  "Guide" => "functions.md"],
#	format = Documenter.HTML(prettyurls = false),
	modules=[Liga]
	)
deploydocs(
	repo = "github.com/evcastelani/Liga.jl.git"
)
