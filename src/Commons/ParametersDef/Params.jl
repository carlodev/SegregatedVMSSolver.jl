abstract type UserParameters end
abstract type MeshInfo <: UserParameters end
abstract type TurbulenceDomain end

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
    u_in::Vector=[1.0,0.0,0.0]
    u_in_mag::Float64 = LinearAlgebra.norm(u_in)
    c::Float64=1.0
    ν::Float64=(c*u_in_mag)/Re
    # @info "Viscosity set ν = $ν"
end


function PhysicalParameters(Re::Int64)
    PhysicalParameters(Re=Re)
end


@with_kw mutable struct TurbulenceParameters{T<:TurbulenceDomain} <: UserParameters
    TI::Float64=0.0
    Re_stress::Matrix=zeros(3,3)
    Eddies::Vector=zeros(3) #Vector{SemEddy} when turbulence at the inlet is active, this is just a placeholder
    Vboxinfo=nothing
    TurbulenceInlet::Bool= (TI==0.0) ? false : true
    Domain::T=Inlet(false)
    R::Float64 = (typeof(Domain)==Inlet) ? 6.0 : Domain.SEM_Boundary_Distance ## radius of the curvature of the SEM domain
end

@with_kw mutable struct Internal <: TurbulenceDomain
    SEM_Boundary_Distance::Real = 1.0
    Create_SEM_Boundary::Bool = true
end

@with_kw mutable struct Inlet <: TurbulenceDomain
    Create_SEM_Boundary::Bool = false
end

function TurbulenceParameters(TI::Float64, Vboxinfo::VirtualBox, physicalp::PhysicalParameters, Create_SEM_Boundary::Bool=false)
    Domain = Inlet()

    if Create_SEM_Boundary
        Domain = Internal()
    end

    TurbulenceParameters(TI, Vboxinfo,physicalp,Domain)
end


function TurbulenceParameters(TI::Float64, Vboxinfo::VirtualBox, physicalp::PhysicalParameters, Domain::TurbulenceDomain; seed = 1234)
    @assert TI>0.0 "Turbulence Intensity Value must be >0.0" 
    Random.seed!(seed)
    Re_stress, Eddies = initialize_eddies(physicalp.u_in_mag, TI, Vboxinfo)
    @info "Eddies initialized - total eddies: $(length(Eddies))"

    TurbulenceParameters(TI=TI, Re_stress=Re_stress, Eddies=Eddies,Vboxinfo=Vboxinfo,Domain=Domain)
end


@with_kw struct SolverParameters <: UserParameters
    θ::Float64=0.5 #theta for velocity ODE solver
    petsc_options::String=petsc_options_default() #solver options for velocity and pressure
    matrix_freq_update::Int64=20 #update matrices and stabilization parameters every fixed number of time steps
    a_err_threshold::Int64=200 ### threshold to be satisfied:: norm(accelation_initial)/norm(acceleration_final) > a_err_threshold
    M::Int64=20 ## maximum number of internal iterations
    Number_Skip_Expansion::Int64=100 #number of initial time steps where is not used the Taylor-Expansion to compute the velocity field at the next time step, but simply the velcity field at the previous one.
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


@with_kw struct InitialParameters <: UserParameters
    u0::Vector = [0.0, 0.0, 0.0]
    restartfile::String = ""
    restart::Bool= (isempty(restartfile)) ? false : true
end

function InitialParameters(restartfile::String)
    InitialParameters(restartfile=restartfile, restart=true)
end

