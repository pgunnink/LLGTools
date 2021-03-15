using LLG
using Documenter

DocMeta.setdocmeta!(LLG, :DocTestSetup, :(using LLG); recursive=true)

makedocs(;
    modules=[LLG],
    authors="PGunnink <pietergunnink@gmail.com> and contributors",
    repo="https://github.com/pgunnink/LLG.jl/blob/{commit}{path}#{line}",
    sitename="LLG.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
