macro publish(mod,name)
    quote
      using SegregatedVMSSolver.$mod: $name; export $name
    end
end

@publish ParametersDef StabilizationMethod
@publish ParametersDef StabilizationFormulation
@publish ParametersDef StabilizationParameters
@publish ParametersDef StabilizedProblem
@publish ParametersDef ScalarStabilization
@publish ParametersDef TensorStabilization
@publish ParametersDef ScalarFormulation
@publish ParametersDef TensorFormulation
@publish ParametersDef VMS
@publish ParametersDef SUPG

@publish ParametersDef SimulationCase
@publish ParametersDef SimulationParameters
@publish ParametersDef Airfoil
@publish ParametersDef WindTunnel
@publish ParametersDef Cylinder
@publish ParametersDef TaylorGreen
@publish ParametersDef LidDriven
@publish ParametersDef VelocityBoundaryCase

@publish ParametersDef TimeParameters
@publish ParametersDef PhysicalParameters
@publish ParametersDef SolverParameters
@publish ParametersDef MeshParameters
@publish ParametersDef CartesianMeshParams
@publish ParametersDef GmshMeshParams
@publish ParametersDef ExportParameters
@publish ParametersDef RestartParameters

@publish Equations segregated_equations

