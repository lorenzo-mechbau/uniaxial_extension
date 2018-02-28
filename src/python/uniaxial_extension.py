#!/usr/bin/env python

# Based on the OpenCMISS-Iron uniaxial extension example, trying to replicate example 521 from classic-cm.

#> Main script
# Add Python bindings directory to PATH
import sys, os

# Intialise OpenCMISS
from opencmiss.iron import iron

def solveModel(compressible, useGeneratedMesh, zeroLoad, usePressureBasis):

    # Set problem parameters - Unit cube
    height = 1.0
    width = 1.0
    length = 1.0

    NumberOfGaussXi = 2

    coordinateSystemUserNumber = 1
    regionUserNumber = 1
    basisUserNumber = 1
    pressureBasisUserNumber = 2
    generatedMeshUserNumber = 1
    meshUserNumber = 1
    decompositionUserNumber = 1
    geometricFieldUserNumber = 1
    fibreFieldUserNumber = 2
    materialFieldUserNumber = 3
    dependentFieldUserNumber = 4
    equationsSetFieldUserNumber = 5
    deformedFieldUserNumber = 6
    equationsSetUserNumber = 1
    problemUserNumber = 1

    # Set all diganostic levels on for testing
    #iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

    numberGlobalXElements = 1
    numberGlobalYElements = 1
    numberGlobalZElements = 1
    totalNumberOfNodes=8
    totalNumberOfElements=1
    InterpolationType = 1
    if(usePressureBasis):
      numberOfMeshComponents = 2
    else:
      numberOfMeshComponents = 1
    if(numberGlobalZElements==0):
        numberOfXi = 2
    else:
        numberOfXi = 3

    # Get the number of computational nodes and this computational node number
    numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
    computationalNodeNumber = iron.ComputationalNodeNumberGet()

    # Create a 3D rectangular cartesian coordinate system
    coordinateSystem = iron.CoordinateSystem()
    coordinateSystem.CreateStart(coordinateSystemUserNumber)
    coordinateSystem.DimensionSet(3)
    coordinateSystem.CreateFinish()

    # Create a region and assign the coordinate system to the region
    region = iron.Region()
    region.CreateStart(regionUserNumber,iron.WorldRegion)
    region.LabelSet("Region")
    region.coordinateSystem = coordinateSystem
    region.CreateFinish()

    # Define basis
    basis = iron.Basis()
    basis.CreateStart(basisUserNumber)
    if InterpolationType in (1,2,3,4):
        basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    elif InterpolationType in (7,8,9):
        basis.type = iron.BasisTypes.SIMPLEX
    basis.numberOfXi = numberOfXi
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
    if(NumberOfGaussXi>0):
        basis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
    basis.CreateFinish()

    if(usePressureBasis):
        # Define pressure basis
        pressureBasis = iron.Basis()
        pressureBasis.CreateStart(pressureBasisUserNumber)
        if InterpolationType in (1,2,3,4):
            pressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
        elif InterpolationType in (7,8,9):
            pressureBasis.type = iron.BasisTypes.SIMPLEX
        pressureBasis.numberOfXi = numberOfXi
        pressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi
        if(NumberOfGaussXi>0):
            pressureBasis.quadratureNumberOfGaussXi = [NumberOfGaussXi]*numberOfXi
        pressureBasis.CreateFinish()

    mesh = iron.Mesh()
    if useGeneratedMesh:
        # Start the creation of a generated mesh in the region
        generatedMesh = iron.GeneratedMesh()
        generatedMesh.CreateStart(generatedMeshUserNumber,region)
        generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
        if(usePressureBasis):
            generatedMesh.basis = [basis,pressureBasis]
        else:
            generatedMesh.basis = [basis]
            generatedMesh.extent = [width,length,height]
            generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]
        # Finish the creation of a generated mesh in the region
        generatedMesh.CreateFinish(meshUserNumber,mesh)

    else:

        # Start the creation of a manually generated mesh in the region
        mesh = iron.Mesh()
        mesh.CreateStart(meshUserNumber,region,numberOfXi)
        mesh.NumberOfComponentsSet(numberOfMeshComponents)
        mesh.NumberOfElementsSet(totalNumberOfElements)

        #Define nodes for the mesh
        nodes = iron.Nodes()
        nodes.CreateStart(region,totalNumberOfNodes)
        nodes.CreateFinish()

        elements = iron.MeshElements()
        meshComponentNumber=1
        elements.CreateStart(mesh,meshComponentNumber,basis)
        elements.NodesSet(1,[1,2,3,4,5,6,7,8])
        elements.CreateFinish()

        mesh.CreateFinish()

    # Create a decomposition for the mesh
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decompositionUserNumber,mesh)
    decomposition.type = iron.DecompositionTypes.CALCULATED
    decomposition.numberOfDomains = numberOfComputationalNodes
    decomposition.CreateFinish()

    # Create a field for the geometry
    geometricField = iron.Field()
    geometricField.CreateStart(geometricFieldUserNumber,region)
    geometricField.MeshDecompositionSet(decomposition)
    geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
    geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
    if InterpolationType == 4:
        geometricField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    geometricField.CreateFinish()

    if useGeneratedMesh:
        # Update the geometric field parameters from generated mesh
        generatedMesh.GeometricParametersCalculate(geometricField)
    else:
        # Update the geometric field parameters manually
        geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
        # node 1
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,1,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,2,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,1,3,0.0)
        # node 2
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,1,height)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,2,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,2,3,0.0)
        # node 3
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,1,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,2,width)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,3,3,0.0)
        # node 4
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,1,height)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,2,width)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,4,3,0.0)
        # node 5
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,1,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,2,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,5,3,length)
        # node 6
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,1,height)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,2,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,6,3,length)
        # node 7
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,1,0.0)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,2,width)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,7,3,length)
        # node 8
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,1,height)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,2,width)
        geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,8,3,length)
        geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

    # Create a fibre field and attach it to the geometric field
    fibreField = iron.Field()
    fibreField.CreateStart(fibreFieldUserNumber,region)
    fibreField.TypeSet(iron.FieldTypes.FIBRE)
    fibreField.MeshDecompositionSet(decomposition)
    fibreField.GeometricFieldSet(geometricField)
    fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
    if InterpolationType == 4:
        fibreField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    fibreField.CreateFinish()

    # Create the material field
    if compressible:
        numberOfMaterialComponents = 3
    else:
        numberOfMaterialComponents = 2
    materialField = iron.Field()
    materialField.CreateStart(materialFieldUserNumber,region)
    materialField.TypeSet(iron.FieldTypes.MATERIAL)
    materialField.MeshDecompositionSet(decomposition)
    materialField.GeometricFieldSet(geometricField)
    materialField.NumberOfVariablesSet(1)
    materialField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfMaterialComponents)
    materialField.VariableLabelSet(iron.FieldVariableTypes.U,"Material")
    materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    if compressible:
        materialField.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 3, 1)
    if InterpolationType == 4:
        materialField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    materialField.CreateFinish()

    # Set Mooney-Rivlin constants c10 and c01 respectively.
    iron.Field.ComponentValuesInitialiseDP(
        materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,2.0)
    iron.Field.ComponentValuesInitialiseDP(
        materialField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,6.0)
    if compressible:
        iron.Field.ComponentValuesInitialiseDP(
            materialField, iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, 3, 1.0e9)

    # Create the dependent field
    if compressible:
        numberOfMaterialComponents = 3
    else:
        numberOfMaterialComponents = 4
    dependentField = iron.Field()
    dependentField.CreateStart(dependentFieldUserNumber,region)
    dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
    dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
    dependentField.MeshDecompositionSet(decomposition)
    dependentField.GeometricFieldSet(geometricField)
    dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
    dependentField.NumberOfVariablesSet(2)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,numberOfMaterialComponents)
    dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,numberOfMaterialComponents)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,1)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,1)
    if not compressible:
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
        dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.ELEMENT_BASED)
        if(usePressureBasis):
            # Set the pressure to be nodally based and use the second mesh component
            if InterpolationType == 4:
                dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,4,iron.FieldInterpolationTypes.NODE_BASED)
                dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,4,iron.FieldInterpolationTypes.NODE_BASED)
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
            dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)
    if InterpolationType == 4:
        dependentField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
    dependentField.CreateFinish()

    # Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
    iron.Field.ParametersToFieldParametersComponentCopy(
        geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
        dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
    if not compressible:
        iron.Field.ComponentValuesInitialiseDP(
            dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,-8.0)

    # Create a deformed geometry field, as cmgui doesn't like displaying
    # deformed fibres from the dependent field because it isn't a geometric field.
    deformedField = iron.Field()
    deformedField.CreateStart(deformedFieldUserNumber, region)
    deformedField.MeshDecompositionSet(decomposition)
    deformedField.TypeSet(iron.FieldTypes.GEOMETRIC)
    deformedField.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
    for component in [1, 2, 3]:
        deformedField.ComponentMeshComponentSet(
                iron.FieldVariableTypes.U, component, 1)
    if InterpolationType == 4:
        deformedField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
    deformedField.CreateFinish()

    # Create the equations_set
    equationsSetField = iron.Field()
    equationsSet = iron.EquationsSet()
    if compressible:
        equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
            iron.EquationsSetTypes.FINITE_ELASTICITY,
            iron.EquationsSetSubtypes.COMPRESSIBLE_FINITE_ELASTICITY]
    else:
        equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
            iron.EquationsSetTypes.FINITE_ELASTICITY,
            iron.EquationsSetSubtypes.MOONEY_RIVLIN]
    equationsSet.CreateStart(equationsSetUserNumber,region,fibreField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
    equationsSet.CreateFinish()

    equationsSet.MaterialsCreateStart(materialFieldUserNumber,materialField)
    equationsSet.MaterialsCreateFinish()

    equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
    equationsSet.DependentCreateFinish()

    # Create equations
    equations = iron.Equations()
    equationsSet.EquationsCreateStart(equations)
    equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
    equations.outputType = iron.EquationsOutputTypes.NONE
    equationsSet.EquationsCreateFinish()

    # Define the problem
    problem = iron.Problem()
    problemSpecification = [iron.ProblemClasses.ELASTICITY,
            iron.ProblemTypes.FINITE_ELASTICITY,
            iron.ProblemSubtypes.NONE]
    problem.CreateStart(problemUserNumber, problemSpecification)
    problem.CreateFinish()

    # Create control loops
    problem.ControlLoopCreateStart()
    problem.ControlLoopCreateFinish()

    # Create problem solver
    nonLinearSolver = iron.Solver()
    linearSolver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
    nonLinearSolver.outputType = iron.SolverOutputTypes.PROGRESS
    nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
    nonLinearSolver.NewtonLinearSolverGet(linearSolver)
    nonLinearSolver.NewtonAbsoluteToleranceSet(1e-14)
    nonLinearSolver.NewtonSolutionToleranceSet(1e-14)
    nonLinearSolver.NewtonRelativeToleranceSet(1e-14)
    linearSolver.linearType = iron.LinearSolverTypes.DIRECT
    #linearSolver.libraryType = iron.SolverLibraries.LAPACK
    problem.SolversCreateFinish()

    # Create solver equations and add equations set to solver equations
    solver = iron.Solver()
    solverEquations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
    solver.SolverEquationsGet(solverEquations)
    solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
    equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
    problem.SolverEquationsCreateFinish()

    # Prescribe boundary conditions (absolute nodal parameters)
    boundaryConditions = iron.BoundaryConditions()
    solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

    #Set x=0 nodes to no x displacment in x. Set x=width nodes to 10% x displacement
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,3,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,1,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,7,1,iron.BoundaryConditionsTypes.FIXED,0.0)

    if zeroLoad:
        load = 0.
    else:
        load = 0.1*width
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,1,iron.BoundaryConditionsTypes.FIXED,load)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,4,1,iron.BoundaryConditionsTypes.FIXED,load)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,6,1,iron.BoundaryConditionsTypes.FIXED,load)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,8,1,iron.BoundaryConditionsTypes.FIXED,load)

    # Set y=0 nodes to no y displacement
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,5,2,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,6,2,iron.BoundaryConditionsTypes.FIXED,0.0)

    # Set z=0 nodes to no y displacement
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,1,3,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,2,3,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,3,3,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundaryConditions.AddNode(dependentField,iron.FieldVariableTypes.U,1,1,4,3,iron.BoundaryConditionsTypes.FIXED,0.0)

    solverEquations.BoundaryConditionsCreateFinish()

    # Solve the problem
    problem.Solve()

    # Copy deformed geometry into deformed field
    for component in [1, 2, 3]:
        dependentField.ParametersToFieldParametersComponentCopy(
            iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, component,
            deformedField, iron.FieldVariableTypes.U,
            iron.FieldParameterSetTypes.VALUES, component)

    if useGeneratedMesh:
        output_file = "./results/unit_cube_generated_mesh"
    else:
        output_file = "./results/unit_cube_manual_mesh"
    if not os.path.exists("./results"):
        os.makedirs("./results")

    prefix = ''
    if compressible:
        prefix = "_compressible"

    if zeroLoad:
        prefix += '_zero_load'

    # Export results
    fields = iron.Fields()
    fields.CreateRegion(region)
    fields.NodesExport(output_file+prefix,"FORTRAN")
    fields.ElementsExport(output_file+prefix,"FORTRAN")
    fields.Finalise()

    problem.Destroy()
    if useGeneratedMesh:
      generatedMesh.Destroy()
    basis.Destroy()
    region.Destroy()
    coordinateSystem.Destroy()

if __name__ == "__main__":
    compressible = False
    useGeneratedMesh = True
    zeroLoad = True
    usePressureBasis = False
    # Arguments: compressible, useGeneratedMesh, zeroLoad, usePressureBasis
    solveModel(False, False, False, False)
    #solveModel(True, False, False, False)
    solveModel(False, True, False, False)
    solveModel(False, False, True, False)
    solveModel(False, True, False, True)
