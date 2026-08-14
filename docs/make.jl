using Jessamine
using Documenter

DocMeta.setdocmeta!(Jessamine, :DocTestSetup, :(using Jessamine); recursive=true)

makedocs(;
    modules=[Jessamine],
    authors="W. G. Mitchener <mitchenerg@charleston.edu> and others",
    sitename="Jessamine.jl",
    format=Documenter.HTML(;
        canonical="https://wgmitchener.github.io/Jessamine.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/wgm-applied-math/Jessamine.jl.git",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#", "dev" =>  "dev"] # Explicitly forces version tracking
)
