abstract type UserParameters end
abstract type MeshInfo <: UserParameters end

struct TimeParameters <: UserParameters
    t0::Float64
    dt::Float64
    tF::Float64
    t_endramp::Float64
    time_average::Bool
    time_window::Tuple{Float64,Float64}
end


function TimeParameters(t0::Float64,dt::Float64,tF::Float64,time_window::Tuple{Float64,Float64})
    @assert tF>t0
    TimeParameters(t0,dt,tF,t0,true,time_window)
end

function TimeParameters(t0::Float64,dt::Float64,tF::Float64,t_endramp::Float64)
    @assert tF>t0
    time_window = (t0,t0)
    TimeParameters(t0,dt,tF,t_endramp,false,time_window)
end

function TimeParameters(t0::Float64,dt::Float64,tF::Float64)
    TimeParameters(t0,dt,tF,t0)
end


struct PhysicalParameters <: UserParameters
    Re::Int64
    ν::Float64
    u_in::Float64
    c::Float64
end

function PhysicalParameters(;Re::Int64,c=1.0, u_in=1.0)
    ν = (c*u_in)/Re
    PhysicalParameters(Re,ν,u_in,c)
end

struct SolverParameters <: UserParameters
    θ::Float64
    petsc_options::String
    matrix_freq_update::Int64
    a_err_threshold::Int64
    M::Int64
end

function SolverParameters(;θ=0.5, petsc_options=petsc_options_default(),matrix_freq_update=20,a_err_threshold=200, M=20 )
    SolverParameters(θ, petsc_options,matrix_freq_update,a_err_threshold, M )
end

struct MeshParameters{T<:MeshInfo} <: UserParameters
    rank_partition::Tuple
    D::Int64
    meshinfo::T
end


struct CartesianMeshParams <:MeshInfo
    N::Vector{Int64}
    L::Vector{Float64}
end

struct GmshMeshParams <:MeshInfo
   filename::String
end

function MeshParameters(rank_partition::Tuple, D::Int64; N::Int64,L::Float64)
    @assert length(rank_partition) == D

    Nt = N .* ones(Int64, D)
    Lt = L .* ones(Float64,D)

    meshinfo = CartesianMeshParams(Nt,Lt)
    MeshParameters(rank_partition, D, meshinfo)
end

function MeshParameters(rank_partition::Int64, D::Int64, meshinfo::GmshMeshParams)
    rank_partition_t = (rank_partition,ones(D-1))
    MeshParameters(rank_partition_t, D, meshinfo)
end


function MeshParameters(rank_partition::Union{Int64,Tuple}, D::Int64, filename::String)
    meshinfo = GmshMeshParams(filename)
    MeshParameters(rank_partition, D, meshinfo)
end


struct ExportParameters <: UserParameters
    printmodel::Bool
    printinitial::Bool
    benchmark::Bool
    log_dir::String
    save_sim_dir::String
    export_field::Bool
    name_tags::Vector{String}
    fieldexport::Vector{Vector{String}}
end

function ExportParameters(;printmodel=true,printinitial=true, benchmark=false,
                            log_dir="Log",  save_sim_dir="Results_vtu",name_tags=[""],fieldexport=[[""]])
    export_field = true
    if isempty(name_tags[1])
        export_field = false

    end

    return ExportParameters(printmodel,printinitial,benchmark,log_dir,save_sim_dir,export_field,name_tags,fieldexport)
end


struct RestartParameters <: UserParameters
    restart::Bool
    restartfile::String
end

function RestartParameters()
    RestartParameters(false, " ")
end

function RestartParameters(rfile::String)
    RestartParameters(true, rfile)
end


