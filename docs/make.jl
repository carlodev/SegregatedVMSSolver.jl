# push!(LOAD_PATH,joinpath("..","src"))

using Documenter, DocumenterCitations, SegregatedVMSSolver


bib = CitationBibliography(joinpath(@__DIR__, "src", "docs.bib"); style=:numeric)


makedocs(;plugins=[bib],
    sitename = "SegregatedVMSSolver.jl",
    modules = [SegregatedVMSSolver],
    pages = [
        "Introduction" => "index.md",
        "Parameters Guide" => [
        "User Parameters"=>"user_params.md",
        "Simulation Parameters" => "sim_params.md"],

        "MPI run" => "mpi.md",
        "Examples" =>[
        "Taylor Green" => "Cases/taylorgreen.md",
        "Lid Driven Cavity Flow" => "Cases/liddriven.md",        
        "Cylinder" => "Cases/cylinder.md",
        "Airfoil" => "Cases/airfoil.md",
        "Create your Own Case" => "create_case.md",
        ],
        "Theory"=>[
            "Incompressible Navier Stokes"=>"Theory/NS_inc.md",
            "FEM"=>"Theory/FEMTheory.md",
            "SUPG & VMS Stabilization"=>"Theory/SUPG_VMS_stab.md"],
        "Tools"=>[       
        "Post Processing" => "post_proc.md",
        "Boundary Layer Initialization" => "blinit.md"],
        "API information" => "api_info.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/carlodev/SegregatedVMSSolver.jl",
    push_preview = true,
)


#julia --project make.jl