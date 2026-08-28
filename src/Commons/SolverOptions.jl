module SolverOptions
using Gridap
using GridapDistributed
using GridapPETSc
using GridapPETSc.PETSC
using PartitionedArrays
using SparseArrays

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


# Wrap the VMS NumericalSetup to allow specialized methods that omit
# garbage collection at every time step.
struct VMSPETScNS{T} <: NumericalSetup
  ns::PETScLinearSolverNS{T}
  Y::Ref{PETScVector}
  B::Ref{PETScVector}
  function VMSPETScNS{T}(ns::PETScLinearSolverNS{T}) where T
    return new(ns, Ref{PETScVector}(), Ref{PETScVector}())
  end
end
VMSPETScNS(ns::PETScLinearSolverNS{T}) where {T} = VMSPETScNS{T}(ns)


"""
    create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)

Wrapper for creating PETSc symbolic and numeric setup for `GridapPETSc`  
"""
function create_PETSc_setup(M::AbstractMatrix,ksp_setup::Function)
      solver = PETScLinearSolver(ksp_setup)
      ss = symbolic_setup(solver, M)
      ns = numerical_setup(ss, M)
      # @check_error_code GridapPETSc.PETSC.KSPView(ns.ksp[],C_NULL)
      return VMSPETScNS(ns)
end

function Algebra.numerical_setup!(vmsns::VMSPETScNS,A::AbstractMatrix)
  ns = vmsns.ns
  ns.A = A
  nnz_a = sum(map(x -> count(!iszero, nonzeros(x)), partition(ns.A)))
  nnz_b = nnz(ns.B)
  if nnz_a != nnz_b
    @info "Updating PETSc Matrix to increase nonzeros from $(nnz_b) to $(nnz_a)"
    ns.B = convert(PETScMatrix,A)
    # @check_error_code PETSC.KSPSetOperators(ns.ksp[],ns.B.mat[],ns.B.mat[])
  else
    GridapPETSc._copy!(ns.B.mat[], ns.A)
  end
  return ns
end

function Algebra.solve!(x::PartitionedArrays.PVector,vmsns::VMSPETScNS,b::PartitionedArrays.PVector)
  ns = vmsns.ns
  X = similar(b,(axes(ns.A)[2],))
  B = similar(b,(axes(ns.A)[2],))
  copy!(X,x)
  copy!(B,b)
  if !isassigned(vmsns.Y)
    vmsns.Y[] = convert(PETScVector,X)
  else
    GridapPETSc._copy!(vmsns.Y[].vec[], X)
  end
  
  solve!(vmsns.Y[],vmsns,B)
  copy!(x,vmsns.Y[])
  return x
end

function Algebra.solve!(x::PETScVector,vmsns::VMSPETScNS,b::AbstractVector)
  ns = vmsns.ns

  if !isassigned(vmsns.B)
    vmsns.B[] = convert(PETScVector,b)
  else
    GridapPETSc._copy!(vmsns.B[].vec[], b)
  end

  solve!(x,ns,vmsns.B[])
  return x
end

end