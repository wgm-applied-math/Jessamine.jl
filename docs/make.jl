using Documenter
# using DocumenterMarkdown
using Jessamine

makedocs(
    modules = [Jessamine],
    sitename = "Documentation for Jessamine"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/wgm-applied-math/Jessamine.jl.git",
    devbranch = "main",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#", "dev" =>  "dev"] # Explicitly forces version tracking
)
