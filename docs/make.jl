#push!(LOAD_PATH,"/home/emerson/evcastelani@uem.br/LiGA/src/")

#push!(LOAD_PATH,"/Users/emersonvitor/Google Drive/LiGA/src/")


using Documenter, Liga

makedocs(
    #modules = [LigaTypes],
    format = :html,
    assets = ["assets/favicon.ico"],
    sitename = "LiGA- Library for Geometric Algebra",
    pages = Any[
        "Overview" => "index.md",
        "Getting Started" =>"gettingstarted.md",
        "Tutorials"=>["types.md","Gkspace.md", "projspace.md","confspace.md","advancedex.md"],
        "Summaries"=> "summary.md"
    ]
)

deploydocs(
    repo = "github.com/evcastelani/Liga.jl.git",
    target = "build",
    julia = "release",
    deps = nothing,
    make = nothing
)