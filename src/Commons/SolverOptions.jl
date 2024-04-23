module SolverOptions
using Gridap
using GridapDistributed
using GridapPETSc

using Gridap.Algebra
using MPI
import GridapPETSc:PETScLinearSolverNS


export vel_kspsetup
export pres_kspsetup
export petsc_options_default
export petsc_options
export create_PETSc_setup
export petsc_options_airfoil


function vel_kspsetup(ksp)
  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"vel_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])

end

function pres_kspsetup(ksp)
  pc = Ref{GridapPETSc.PETSC.PC}()
  @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"pres_")
  @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])
end

function petsc_options_default()
  return petsc_options()
end

"""
    petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp = "cg", pres_pc = "gamg")
It provides the command-line for `GridapPETSc` to solve the segregated linear systems
"""
function petsc_options(; vel_ksp="gmres", vel_pc="gamg", pres_ksp = "cg", pres_pc = "gamg")
  return " -vel_ksp_type $(vel_ksp) -vel_pc_type $(vel_pc) -vel_ksp_rtol 1.e-10 -vel_ksp_converged_reason \
  -pres_ksp_type $(pres_ksp) -pres_pc_type $(pres_pc)  -pres_ksp_rtol 1.e-6 -pres_ksp_converged_reason \
  -ksp_atol 0.0"
end

function petsc_options_airfoil()
  return "-vel_ksp_type gmres -vel_ksp_gmres_restart 300  -vel_ksp_rtol 1.e-6 -vel_pc_type hypre -vel_pc_hypre_type euclid -vel_ksp_converged_reason \
        -pres_ksp_type cg -pres_pc_type gamg -pres_ksp_rtol 1.e-2 -pres_ksp_converged_reason -ksp_atol 0.0"
end

"""
    create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)

Wrapper for creating PETSc symbolic and numeric setup for `GridapPETSc`  
"""
function create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)
      solver = PETScLinearSolver(ksp_setup)
      ss = symbolic_setup(solver, M)
      ns = numerical_setup(ss, M)
      @check_error_code GridapPETSc.PETSC.KSPView(ns.ksp[],C_NULL)

      return ns
end


function Algebra.solve!(x::PETScVector,ns::PETScLinearSolverNS,b::AbstractVector)
  if (x.comm != MPI.COMM_SELF)
    # gridap_petsc_gc() # Do garbage collection of PETSc objects
  end

  B = convert(PETScVector,b)
  solve!(x,ns,B)
  x
end


end