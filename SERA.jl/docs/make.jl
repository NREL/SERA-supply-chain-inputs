# Add the SERA package to the path.
push!(LOAD_PATH, "../src/")

# Import modules.
using Documenter
using SERA,
      SERA.Supply,
      SERA.Supply.IO,
      SERA.Supply.Types,
      SERA.Util

# Build the documentation.
makedocs(
    sitename="SERA.jl Documentation",
    modules=[
        SERA
        SERA.Supply
        SERA.Supply.IO
        SERA.Supply.Types
        SERA.Util
    ]
)
