
using Test
using SegregatedVMSSolver

function create_simulation_test(case, D; meshfile="", rank_partition=(2,2))
      t0 =0.0
      dt = 0.1
      tF = 0.3
      Re = 1000
     
      sprob = StabilizedProblem()
      timep = TimeParameters(t0,dt,tF)
      physicalp = PhysicalParameters(Re=Re)
      solverp = SolverParameters(M=2)
      exportp = ExportParameters(printmodel=false)

      if isempty(meshfile)
            meshp = MeshParameters(rank_partition,D; N = 10, L=0.5)
      else
            meshp= MeshParameters(rank_partition,D,meshfile)
      end    
      
      simparams = SimulationParameters(timep,physicalp,solverp,exportp)
        
      new_case = case(meshp,simparams,sprob)
      return new_case

end

function iterate_test_cases(D::Int64)
      if D==2
            cases = [Airfoil,Cylinder,LidDriven,TaylorGreen]
            airfoil_mesh_file = joinpath(@__DIR__, "..", "models", "DU89_2D_A1_M.msh")
            cylinder_mesh_file = joinpath(@__DIR__, "..", "models", "Cylinder_2D.msh")
            mesh_files = [airfoil_mesh_file,cylinder_mesh_file,"",""]
      elseif D ==3
            cases = [Airfoil]
            airfoil_mesh_file = joinpath(@__DIR__, "..", "models", "sd7003s_3D_simple.msh")
            mesh_files = [airfoil_mesh_file]

      end

      return zip(cases,mesh_files)
end



"""
    create_turbulent_simulation_test()

Testing an `Airfoil` simulation where the inlet is set as turbulent
"""
function create_turbulent_simulation_test()
      TI = 0.001    
      D = 3
      meshfile = joinpath(@__DIR__, "..", "models", "sd7003s_3D_simple.msh")
      rank_partition=(2,2,1)
  
      #VirtualBox creation
      Vboxinfo = VirtualBox((-6,6), (-0.1,0.3); σ=0.02)
      physicalp = PhysicalParameters(Re=1000)
  
      turbulencep= TurbulenceParameters(TI, Vboxinfo, physicalp)
      sprob = StabilizedProblem()
  
      timep = TimeParameters(0.0,0.1, 1.0)
  
      solverp = SolverParameters(M=2)
      exportp = ExportParameters(printmodel=false)
      meshp= MeshParameters(rank_partition,D,meshfile)
  
      simparams = SimulationParameters(timep,physicalp,turbulencep,solverp,exportp)
      tcase = Airfoil(meshp,simparams,sprob)
  
      return tcase

end


#mpiexecjl --project=../. -n 4 julia case_test.jl

