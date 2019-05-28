using Documenter

makedocs(format = :html,
	assets = ["assets/favicon.ico"],
	sitename = "Liga- Library for Geometric Algebra",
	pages =["Overview" => "index.md",
			  "Getting Started" => "gettingstarted.md",
			  "Guide" => "functions.md"],
	#html_prettyurls = false
	)
deploydocs(
	repo = "github.com/evcastelani/Liga.jl.git"
)
