using YBusRevisited
using Documenter

DocMeta.setdocmeta!(YBusRevisited, :DocTestSetup, :(using YBusRevisited); recursive=true)

makedocs(;
    modules=[YBusRevisited],
    authors="Hantao Cui",
    repo="https://github.com/YBusRevisited/YBusRevisited.jl/blob/{commit}{path}#{line}",
    sitename="YBusRevisited.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://YBusRevisited.github.io/YBusRevisited.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/YBusRevisited/YBusRevisited.jl",
    devbranch="main",
)
