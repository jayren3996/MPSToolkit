using Documenter
using MPSToolkit

# Make the package available inside `jldoctest` blocks so docstring examples are executed
# and verified during the docs build (`makedocs` runs doctests by default).
DocMeta.setdocmeta!(MPSToolkit, :DocTestSetup, :(using MPSToolkit); recursive=true)

const PRETTY_URLS = get(ENV, "CI", "false") == "true"

makedocs(
  sitename="MPSToolkit.jl",
  authors="Jie Ren and contributors",
  modules=[MPSToolkit],
  clean=true,
  checkdocs=:none,
  format=Documenter.HTML(
    prettyurls=PRETTY_URLS,
    canonical="https://jayren3996.github.io/MPSToolkit/stable",
    edit_link="main",
    assets=["assets/favicon.ico"],
  ),
  pages=[
    "Home" => "index.md",
    "Getting Started" => "getting-started.md",
    "Manual" => [
      "Architecture" => "manual/architecture.md",
      "TEBD" => "manual/tebd.md",
      "TDVP" => "manual/tdvp.md",
      "ScarFinder" => "manual/scarfinder.md",
      "Operator Space" => "manual/operator-space.md",
      "DMT" => "manual/dmt.md",
      "DAOE" => "manual/daoe.md",
      "Chebyshev" => "manual/chebyshev.md",
    ],
    "Examples" => "examples.md",
    "API Reference" => "api.md",
  ],
)

deploydocs(
  repo="github.com/jayren3996/MPSToolkit.git",
  devbranch="main",
)
