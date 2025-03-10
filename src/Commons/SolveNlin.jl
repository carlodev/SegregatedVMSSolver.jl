module SolveNlin
using Gridap
using GridapDistributed
using GridapPETSc
using SparseArrays
using PartitionedArrays
using Parameters
using Gridap.FESpaces

using SegregatedVMSSolver.ParametersDef
using SegregatedVMSSolver.SolverOptions
using SegregatedVMSSolver.CreateProblem
using SegregatedVMSSolver.MatrixCreation
using SegregatedVMSSolver.VectorsOperations
using SegregatedVMSSolver.ExportUtility
using SegregatedVMSSolver.Interfaces
using SegregatedVMSSolver.Equations

export solve_nlin

function intialize_fields(U,P,V,Q)
    Y = MultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])
    return X,Y
end

function initialize_solve(simcase::SimulationCase,params::Dict{Symbol,Any})
        @unpack trials,tests = params
   
        @sunpack t0,dt,save_sim_dir = simcase
        
        U,P = trials
        V,Q = tests

        X,Y = intialize_fields(U,P,V,Q)

        Ut0 = U(t0)
        Pt0 = P(t0)
      
        Ut0_1 = U(t0+dt)
        Pt0_1 = P(t0+dt)
      
        merge!(params, Dict(:Utn => Ut0, :Ptn => Pt0, :Utn1 => Ut0_1, :Ptn1 => Pt0_1))
      
      
        uh0, ph0 = create_initial_conditions(simcase,params)
        xh0 = interpolate_everywhere([uh0, ph0], X(t0))


    return xh0, X,Y

end

function solve_nlin(simcase::SimulationCase,params::Dict{Symbol,Any})
    xh0, X,Y = initialize_solve(simcase,params)
    
    res, jac, jac_t = coupled_equations(params,simcase)

    op = TransientFEOperator(res, jac, jac_t,X,Y)
    @sunpack t0,dt,tF,petsc_options,θ = simcase

        GridapPETSc.with(args=split(petsc_options)) do

        nls = PETScNonlinearSolver()


        ode_solver = ThetaMethod(nls,dt,θ)

        nl_solution = solve(ode_solver, op,xh0, t0, tF)
        
        for ((uh_tn,ph_tn), tn) in nl_solution
            println("Solution at time $(tn) computed")
    
        end

    end #end GridapPETSc

end

end #end module