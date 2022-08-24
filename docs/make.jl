using WienerHopfScalar
using Documenter

DocMeta.setdocmeta!(WienerHopfScalar, :DocTestSetup, :(using WienerHopfScalar); recursive=true)

makedocs(;
    modules=[WienerHopfScalar],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/WienerHopfScalar.jl/blob/{commit}{path}#{line}",
    sitename="WienerHopfScalar.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/WienerHopfScalar.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/WienerHopfScalar.jl",
    devbranch="main",
)
