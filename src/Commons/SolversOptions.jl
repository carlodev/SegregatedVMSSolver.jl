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
