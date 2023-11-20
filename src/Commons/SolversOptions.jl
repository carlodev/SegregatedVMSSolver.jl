function vel_kspsetup(ksp)
    pc = Ref{GridapPETSc.PETSC.PC}()
    @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"vel_")
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])

    # @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[], GridapPETSc.PETSC.KSPPREONLY)
    # @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[], pc)
    # @check_error_code GridapPETSc.PETSC.PCSetType(pc[], GridapPETSc.PETSC.PCJACOBI)  
    # @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[], 1e-10,0.0,1e5,1e4)  

  end

  function pres_kspsetup(ksp)
    pc = Ref{GridapPETSc.PETSC.PC}()
    @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"pres_")
    @check_error_code GridapPETSc.PETSC.KSPSetFromOptions(ksp[])

    # @check_error_code GridapPETSc.PETSC.KSPSetType(ksp[], GridapPETSc.PETSC.KSPGMRES)
    # @check_error_code GridapPETSc.PETSC.KSPGetPC(ksp[], pc)
    # @check_error_code GridapPETSc.PETSC.PCSetType(pc[], GridapPETSc.PETSC.PCGAMG)
    # @check_error_code GridapPETSc.PETSC.KSPSetTolerances(ksp[],1e-10, 1e-8,1e4,1000)

  end


  # function vel_kspsetup(ksp)
  #   @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"vel_")
  #   @check_error_code GridapPETSc.PETSC.KSPSetReusePreconditioner(ksp[], GridapPETSc.PETSC.PETSC_TRUE)
    

  # end
  
  # function pres_kspsetup(ksp)
  #   @check_error_code GridapPETSc.PETSC.KSPSetOptionsPrefix(ksp[],"pres_")
  #   @check_error_code GridapPETSc.PETSC.KSPSetReusePreconditioner(ksp[], GridapPETSc.PETSC.PETSC_TRUE)

  # end