abstract type SimulationCase end
abstract type MeshFileCase <:SimulationCase end
abstract type CartesianCase <:SimulationCase end


struct SimulationParameters
    timep::TimeParameters
    physicalp::PhysicalParameters
    solverp::SolverParameters
    exportp::ExportParameters
    clusterp::ClusterParameters
    restartp::RestartParameters
end

function SimulationParameters(timep::TimeParameters,physicalp::PhysicalParameters, 
                            solverp::SolverParameters,exportp::ExportParameters, clusterp::ClusterParameters)
    restartp=RestartParameters()
    SimulationParameters(timep,physicalp,solverp,exportp,clusterp,restartp)
end

struct Airfoil <: MeshFileCase 
    meshfile::String
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

struct WindTunnel <: MeshFileCase 
    meshfile::String
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

struct Cylinder <: MeshFileCase 
    meshfile::String
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

struct TaylorGreen <: CartesianCase
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

struct LidDriven <: CartesianCase 
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

VelocityBoundaryCase = Union{Airfoil,WindTunnel,Cylinder, LidDriven}

MyStructurePrint = Union{SimulationCase,StabilizedProblem,SimulationParameters,UserParameters,StabilizationMethod,StabilizationFormulation}

function printstructure(s::MyStructurePrint)
    fnames = fieldnames(typeof(s))
    for fn in fnames
        fnfield = getfield(s,fn)
        fnfield_type = typeof(fnfield) <: MyStructurePrint
        if fnfield_type
            println(typeof(fnfield))
            printstructure(fnfield)
        else
            println("$fn = $fnfield")
        end
    end 
    println()
end

Base.show(io::IO,s::MyStructurePrint) = printstructure(s)



function get_field(s::MyStructurePrint,f::Symbol)
    flag, val = search_field(s::MyStructurePrint,f::Symbol, false, nothing)
    if flag==false
        @error "field $f not found in $s"
    else
        return val
    end
end

function get_field(s::MyStructurePrint,fv::Vector{Symbol})
    return map(f->get_field(s,f),fv)
end

function search_field(s::MyStructurePrint,f::Symbol, flag::Bool, a)
    fnames = fieldnames(typeof(s))

    i = 1 
    while i<= length(fnames) && !flag
        fn = fnames[i]

        fnfield = getfield(s,fn)
        fnfield_type = typeof(fnfield) <: MyStructurePrint
        if f == fn
            flag= true
            a = fnfield
        elseif fnfield_type
            flag, a = search_field(fnfield,f, flag, a)
        end

        i = i+1
    end 


    return flag, a
end


