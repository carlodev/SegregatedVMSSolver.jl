push!(LOAD_PATH,joinpath("..","src"))

using Documenter, DocumenterCitations, SegregatedVMSSolver

bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)


makedocs(bib;
    sitename = "SegregatedVMSSolver.jl",
    modules = [SegregatedVMSSolver],
    pages = [
        "Introduction" => "index.md",
        "Running Simulation" => "run.md",
        "Post Processing" => "post_proc.md",
        "Boundary Layer Initialization" => "blinit.md",
        "API information" => "api_info.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/carlodev/SegregatedVMSSolver.jl",
    push_preview = true,
)
