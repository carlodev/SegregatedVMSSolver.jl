abstract type UserParameters end
abstract type MeshInfo <: UserParameters end

@with_kw struct TimeParameters <: UserParameters
    t0::Float64 = 0.0
    dt::Float64
    tF::Float64
    t_endramp::Float64 = t0
    time_window::Tuple{Float64,Float64} = (t0,t0)
    time_average::Bool = (time_window[1]==time_window[2]) ? false : true
    @assert t0<tF "t0 = $t0, tF = $tF, condition t0<tF not satisfied"
    @assert time_window[1] >= t0 && time_window[2] <= tF "time_window exceeds time limits: t0 = $t0, tF = $tF, time_window=$(time_window)"
    @assert t0<=t_endramp<tF "$(t_endramp) not valid"
end

function TimeParameters(t0::Float64,dt::Float64,tF::Float64)
    TimeParameters(t0=t0,dt=dt,tF=tF)
end


@with_kw struct PhysicalParameters <: UserParameters
    Re::Int64
    u_in::Float64=1.0
    c::Float64=1.0
    ν::Float64=(c*u_in)/Re
    # @info "Viscosity set ν = $ν"
end


@with_kw struct SolverParameters <: UserParameters
    θ::Float64=0.5
    petsc_options::String=petsc_options_default()
    matrix_freq_update::Int64=20
    a_err_threshold::Int64=200
    M::Int64=20
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


@with_kw struct ExportParameters <: UserParameters
    printmodel::Bool=true
    printinitial::Bool=true
    benchmark::Bool=false
    log_dir::String="Log"
    save_sim_dir::String="Results_vtu"
    name_tags::Vector{String}=[""]
    fieldexport::Vector{Vector{String}}=[[""]]
    export_field::Bool = (isempty(name_tags[1])) ? false : true
    @assert length(name_tags) == length(fieldexport)
end


@with_kw struct RestartParameters <: UserParameters
    restartfile::String = ""
    restart::Bool= (isempty(restartfile)) ? false : true
end





