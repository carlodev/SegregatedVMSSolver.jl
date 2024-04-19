
using Test
using SegregatedVMSSolver



function create_simulation_test(case, D; meshfile="", rank_partition=(2,2))
      t0 =0.0
      dt = 0.1
      tF = 1.0
      Re = 1000
     
      sprob = StabilizedProblem()
      timep = TimeParameters(t0,dt,tF)
      physicalp = PhysicalParameters(Re=Re)
      solverp = SolverParameters()
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

#mpiexecjl --project=../. -n 4 julia case_test.jl

