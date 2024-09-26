abstract type SimulationCase end
abstract type VelocityBoundaryCase <: SimulationCase end


abstract type TGVBoundaryConditionsType end


struct SimulationParameters{T<:TurbulenceDomain}
    timep::TimeParameters
    physicalp::PhysicalParameters
    turbulencep::TurbulenceParameters{T}
    solverp::SolverParameters
    exportp::ExportParameters
    intialp::InitialParameters
end

function SimulationParameters(timep::TimeParameters, physicalp::PhysicalParameters,
    solverp::SolverParameters, exportp::ExportParameters)
    turbulencep = TurbulenceParameters()
    intialp = InitialParameters()
    SimulationParameters(timep, physicalp, turbulencep, solverp, exportp, intialp)
end

function SimulationParameters(timep::TimeParameters, physicalp::PhysicalParameters, turbulencep::TurbulenceParameters,
    solverp::SolverParameters, exportp::ExportParameters)
    intialp = InitialParameters()
    SimulationParameters(timep, physicalp, turbulencep, solverp, exportp, intialp)
end

function SimulationParameters(timep::TimeParameters, physicalp::PhysicalParameters,
    solverp::SolverParameters, exportp::ExportParameters, intialp::InitialParameters)
    turbulencep = TurbulenceParameters()
    SimulationParameters(timep, physicalp, turbulencep, solverp, exportp, intialp)
end

"""
    create_new_case(case::Symbol)

It creates a new case, named `case`. It is useful to add new simulation cases in a high level way.
"""
function create_new_case(case::Symbol)
    @eval begin
        struct $case <: VelocityBoundaryCase
            meshp::MeshParameters
            simulationp::SimulationParameters
            sprob::StabilizedProblem
        end
        export $case
    end
end

for case in (:Airfoil, :WindTunnel, :Cylinder, :LidDriven)
    create_new_case(case::Symbol)
end



@with_kw struct TaylorGreenParameters
    Vs::Float64 = 1.0
    # diameter::Float64 = 0.5 # It is defined in the PhysicalParameter c
    Ua::Float64 = 0.3
    Va::Float64 = 0.2
    # ν::Float64 = 0.001  # It is now defined in the PhysicalParameter
end


function TaylorGreen_Periodic_Parameters()
    Vs::Float64 = 1.0
    # diameter::Float64 = 0.5 # It is defined in the PhysicalParameter c
    Ua::Float64 = 0.3
    Va::Float64 = 0.2
    # ν::Float64 = 0.001  # It is now defined in the PhysicalParameter
    return TaylorGreenParameters(Vs, Ua, Va)
end


function TaylorGreen_Natural_Parameters()
    Vs::Float64 = 1.0
    # diameter::Float64 = 0.5 # It is defined in the PhysicalParameter c
    Ua::Float64 = 0.0
    Va::Float64 = 0.0
    # ν::Float64 = 0.001  # It is now defined in the PhysicalParameter
    return  TaylorGreenParameters(Vs, Ua, Va)

end


@with_kw mutable struct Periodic <: TGVBoundaryConditionsType
    params::TaylorGreenParameters = TaylorGreen_Periodic_Parameters()
    bc::Tuple = (true, true)
    a_solution::Dict 
end

@with_kw mutable struct Natural <: TGVBoundaryConditionsType
    params::TaylorGreenParameters = TaylorGreen_Natural_Parameters()
    bc::Tuple = (false, false)
    a_solution::Dict 
end

function Periodic(mp::MeshParameters, fp::PhysicalParameters) 
    params = TaylorGreen_Periodic_Parameters()
    return Periodic(mp, fp, params)
end

function Periodic(mp::MeshParameters, fp::PhysicalParameters, params::TaylorGreenParameters) 
    @unpack D = mp
    @unpack Vs, Ua, Va = params
    @unpack ν,c = fp

    @assert D == 2 "TGV Periodic 3D not supported yet"

    bc = ntuple(i -> true, D)

    velocity, pressure, _ = analytical_solution(c*0.5, Vs, Ua, Va, ν)
    a_solution = Dict(:velocity => velocity, :pressure => pressure)
    return Periodic(params, bc, a_solution )
end

function Natural(mp::MeshParameters, fp::PhysicalParameters) 
    params = TaylorGreen_Natural_Parameters()
    return Natural(mp, fp, params)
end

function Natural(mp::MeshParameters, fp::PhysicalParameters, params::TaylorGreenParameters) 
    @unpack D = mp
    @unpack L = mp.meshinfo
    @unpack Vs, Ua, Va = params
    @unpack ν,c = fp
    
    bc = ntuple(i -> false, D)

    velocity, pressure  = TGV_initial(Vs, D, L)
    a_solution = Dict(:velocity => velocity, :pressure => pressure)

    return Natural(params, bc, a_solution)
end


struct TaylorGreen{T<:TGVBoundaryConditionsType} <: SimulationCase
    bc_type::T
    meshp::MeshParameters
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

MyStructurePrint = Union{SimulationCase,StabilizedProblem,SimulationParameters,
    UserParameters,MeshInfo,StabilizationMethod,StabilizationFormulation}


function printstructure(s::MyStructurePrint)
    fnames = fieldnames(typeof(s))
    for fn in fnames
        fnfield = getfield(s, fn)
        fnfield_type = typeof(fnfield) <: MyStructurePrint
        if fnfield_type
            println(typeof(fnfield))
            printstructure(fnfield)
        elseif typeof(fnfield) <: Vector{SemEddy}
            println("$fn = $(fnfield[1]) - total Eddies $(length(fnfield))")
        else
            println("$fn = $fnfield")
        end
    end
    println()
end

Base.show(io::IO, s::MyStructurePrint) = printstructure(s)


function search_field(s::MyStructurePrint, ::Val{f}, flag::Bool, a) where {f}
    fnames = fieldnames(typeof(s))
    i = 1
    while i <= length(fnames) && !flag
        fn = fnames[i]

        fnfield = getfield(s, fn)
        fnfield_type = typeof(fnfield) <: MyStructurePrint

        if f == fn
            flag = true
            a = fnfield
        elseif fnfield_type
            flag, a = search_field(fnfield, Val{f}(), flag, a)
        end

        i = i + 1
    end

    return flag, a
end


### Macro inspired from @unpack from UnPack.jl package
macro sunpack(args)
    args.head != :(=) && error("Expression needs to be of form `a, b = c`")
    items, suitecase = args.args
    items = isa(items, Symbol) ? [items] : items.args
    kd = Vector{Expr}(undef, length(items))

    for (i, key) in enumerate(items)
        kd[i] = quote
            flag, val = search_field($suitecase, Val{$(Expr(:quote, key))}(), false, nothing)
            $key = val
        end

    end
    kdblock = Expr(:block, kd...)
    esc(kdblock)
end


"""
    compute_fluctuation(x,t, simcase::SimulationCase)

For the point x, at the time t it computes the velocity fluctuations in all the direction. Each time the time t is increased the Eddy are convected.
The time informations are coded in the VirtualBox.
"""
function compute_fluctuation(x, t, (TurbulenceInlet, Eddies, u_in_mag, Vboxinfo, Re_stress, D, R))
   u_fluct = zeros(D)
   
   R_tol = 0.01



   if x[1]< 0.0 &&  R^2 - R_tol < x[1].^2 + x[2].^2 < R^2 + R_tol

        xplane = [x...]
        if D == 2
            xplane = [xplane..., 0.0]
        end

        xplane[1] = 0.0 #As the inlet was a vertical plane, reduce the dimensions of virtual box
        u_fluct = compute_fluct(xplane, t, Eddies, u_in_mag, Vboxinfo, Re_stress)
        u_fluct[1] = u_fluct[1] - u_in_mag #remove the inlet velocity component
        end 

    return VectorValue(u_fluct[1:D]...)

end
