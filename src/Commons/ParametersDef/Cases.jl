abstract type SimulationCase end

struct SimulationParameters
    timep::TimeParameters
    physicalp::PhysicalParameters
    solverp::SolverParameters
    exportp::ExportParameters
    restartp::RestartParameters
end

function SimulationParameters(timep::TimeParameters,physicalp::PhysicalParameters, 
                            solverp::SolverParameters,exportp::ExportParameters)
    restartp=RestartParameters()
    SimulationParameters(timep,physicalp,solverp,exportp,restartp)
end


for case in (:Airfoil,:WindTunnel,:Cylinder,:LidDriven)
    @eval begin
        struct $case <: SimulationCase
            meshp::MeshParameters
            simulationp::SimulationParameters
            sprob::StabilizedProblem
        end
    end
end
  
struct TaylorGreen <: SimulationCase
    analyticalsol
    meshp::MeshParameters
    simulationp::SimulationParameters
    sprob::StabilizedProblem
end

function TaylorGreen(meshp::MeshParameters, simulationp::SimulationParameters, sprob::StabilizedProblem)
    diameter = 0.5 #0.5 [m] vortex dimension
    Vs = 1 #1[m/s]swirling speed
    Ua = 0.3 #0.3 [m/s]convective velocity in x
    Va = 0.2 #0.2 [m/s]convective velocity in y
    ν = 0.001 #0.001 m2/s 
  
  
    velocity, pressure, ωa = analytical_solution(diameter, Vs, Ua, Va, ν)
    analyticalsol = Dict(:velocity=>velocity, :pressure=>pressure)
    TaylorGreen(analyticalsol, meshp, simulationp, sprob)
end

VelocityBoundaryCase = Union{Airfoil,WindTunnel,Cylinder, LidDriven}


MyStructurePrint = Union{SimulationCase,StabilizedProblem,SimulationParameters,
                        UserParameters,MeshInfo, StabilizationMethod,StabilizationFormulation}


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



function search_field(s::MyStructurePrint, ::Val{f}, flag::Bool, a) where {f}
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
            flag, a = search_field(fnfield, Val{f}(), flag, a)
        end

        i = i+1
    end 

    return flag, a
end



### Macro inspired from @unpack from UnPack.jl package
macro sunpack(args)
    args.head!=:(=) && error("Expression needs to be of form `a, b = c`")
    items, suitecase = args.args
    items = isa(items, Symbol) ? [items] : items.args
    kd = Vector{Expr}(undef, length(items))
    
    for (i,key) in enumerate(items)
        kd[i] = quote
            flag, val = search_field($suitecase, Val{$(Expr(:quote, key))}(), false, nothing)
            $key = val
        end
    
    end
    kdblock = Expr(:block, kd...)
    esc(kdblock)
end