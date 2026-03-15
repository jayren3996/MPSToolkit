using Documenter
using MPSToolkit

const PRETTY_URLS = get(ENV, "CI", "false") == "true"

makedocs(
  sitename="MPSToolkit.jl",
  modules=[MPSToolkit],
  clean=true,
  checkdocs=:none,
  format=Documenter.HTML(
    prettyurls=PRETTY_URLS,
    canonical="https://jayren3996.github.io/MPSToolkit.jl/stable",
    edit_link="main",
  ),
  pages=[
    "Home" => "index.md",
    "Getting Started" => "getting-started.md",
    "Manual" => [
      "Architecture" => "manual/architecture.md",
      "Evolution" => "manual/evolution.md",
      "Operator Space" => "manual/operator-space.md",
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
