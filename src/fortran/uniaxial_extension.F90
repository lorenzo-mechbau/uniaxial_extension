PROGRAM UNIAXIAL_EXTENSION

#ifndef NOMPIMOD
  USE MPI
#endif
  USE OpenCMISS
  USE OpenCMISS_Iron

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  TYPE(cmfe_RegionType)           :: worldRegion
  TYPE(cmfe_CoordinateSystemType) :: worldCoordinateSystem

  INTEGER(CMISSIntg)              :: Err

  ! Intialise cmiss
  CALL cmfe_Initialise(worldCoordinateSystem,worldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)

  ! Input arguments: compressible, useGeneratedMesh, zeroLoad, useSimplex, usePressureBasis

  CALL SOLVE_MODEL(.FALSE., .FALSE., .FALSE., .FALSE., .FALSE.)
  ! CALL SOLVE_MODEL(.TRUE., .FALSE. , .FALSE., .FALSE., .FALSE.)
  CALL SOLVE_MODEL(.FALSE., .TRUE. , .FALSE., .FALSE., .FALSE.)
  CALL SOLVE_MODEL(.FALSE., .FALSE., .TRUE., .FALSE., .FALSE.)
  CALL SOLVE_MODEL(.FALSE., .FALSE., .FALSE., .TRUE., .FALSE.)
  CALL SOLVE_MODEL(.FALSE., .TRUE., .FALSE., .FALSE., .TRUE.)

  CALL cmfe_Finalise(Err)

  STOP

CONTAINS

  SUBROUTINE SOLVE_MODEL(compressible, useGeneratedMesh, zeroLoad, useSimplex, usePressureBasis)

    LOGICAL, INTENT(IN)             :: compressible
    LOGICAL, INTENT(IN)             :: useGeneratedMesh
    LOGICAL, INTENT(IN)             :: zeroLoad
    LOGICAL, INTENT(IN)             :: useSimplex
    LOGICAL, INTENT(IN)             :: usePressureBasis

    REAL(CMISSRP)                   :: height = 1.0_CMISSRP
    REAL(CMISSRP)                   :: width = 1.0_CMISSRP
    REAL(CMISSRP)                   :: length = 1.0_CMISSRP
    REAL(CMISSRP)                   :: load = 0.0_CMISSRP

    LOGICAL                         :: directory_exists = .FALSE.

    INTEGER(CMISSIntg)              :: numberOfGaussXi = 2

    INTEGER(CMISSIntg)              :: coordinateSystemUserNumber = 1
    INTEGER(CMISSIntg)              :: regionUserNumber = 1
    INTEGER(CMISSIntg)              :: basisUserNumber = 1
    INTEGER(CMISSIntg)              :: pressureBasisUserNumber = 2
    INTEGER(CMISSIntg)              :: generatedMeshUserNumber = 1
    INTEGER(CMISSIntg)              :: meshUserNumber = 1
    INTEGER(CMISSIntg)              :: decompositionUserNumber = 1
    INTEGER(CMISSIntg)              :: geometricFieldUserNumber = 1
    INTEGER(CMISSIntg)              :: fibreFieldUserNumber = 2
    INTEGER(CMISSIntg)              :: materialFieldUserNumber = 3
    INTEGER(CMISSIntg)              :: dependentFieldUserNumber = 4
    INTEGER(CMISSIntg)              :: equationsSetFieldUserNumber = 5
    INTEGER(CMISSIntg)              :: deformedFieldUserNumber = 6
    INTEGER(CMISSIntg)              :: equationsSetUserNumber = 1
    INTEGER(CMISSIntg)              :: problemUserNumber = 1
    INTEGER(CMISSIntg)              :: equationsSetIndex = 1

    INTEGER(CMISSIntg)              :: numberGlobalXElements = 1
    INTEGER(CMISSIntg)              :: numberGlobalYElements = 1
    INTEGER(CMISSIntg)              :: numberGlobalZElements = 1
    INTEGER(CMISSIntg)              :: totalNumberOfNodes = 8
    INTEGER(CMISSIntg)              :: totalNumberOfElements = 1
    INTEGER(CMISSIntg)              :: InterpolationType
    INTEGER(CMISSIntg)              :: numberOfMeshComponents = 1
    INTEGER(CMISSIntg)              :: meshComponentNumber = 1

    INTEGER(CMISSIntg)              :: numberOfComputationalNodes,computationalNodeNumber
    INTEGER(CMISSIntg)              :: componentIdx,Err,numberOfMaterialComponents
    INTEGER(CMISSIntg)              :: numberOfXi,quadratureOrder

    CHARACTER(LEN=255)              :: output_file,prefix

    !CMISS variables

    TYPE(cmfe_BasisType)                  :: basis,pressureBasis
    TYPE(cmfe_BoundaryConditionsType)     :: boundaryConditions
    TYPE(cmfe_ComputationEnvironmentType) :: computationEnvironment
    TYPE(cmfe_CoordinateSystemType)       :: coordinateSystem
    TYPE(cmfe_DecompositionType)          :: decomposition
    TYPE(cmfe_EquationsType)              :: equations
    TYPE(cmfe_EquationsSetType)           :: equationsSet
    TYPE(cmfe_FieldType)                  :: geometricField,equationsSetField,fibreField
    TYPE(cmfe_FieldType)                  :: dependentField,materialField,deformedField
    TYPE(cmfe_FieldsType)                 :: fields
    TYPE(cmfe_MeshType)                   :: mesh
    TYPE(cmfe_GeneratedMeshType)          :: generatedMesh
    TYPE(cmfe_MeshElementsType)           :: elements
    TYPE(cmfe_NodesType)                  :: nodes
    TYPE(cmfe_ProblemType)                :: problem
    TYPE(cmfe_RegionType)                 :: region
    TYPE(cmfe_SolverType)                 :: solver,nonlinearSolver,linearSolver
    TYPE(cmfe_SolverEquationsType)        :: solverEquations
    TYPE(cmfe_ControlLoopType)            :: controlLoop

    WRITE(*,'(A)') "Program starting."

    ! Set all diganostic levels on for testing
    CALL cmfe_DiagnosticsSetOn(CMFE_FROM_DIAG_TYPE,[1,2,3,4,5],"Diagnostics", &
      & ["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"],Err)

    IF (usePressureBasis) THEN
      numberOfMeshComponents = 2
    ELSE
      numberOfMeshComponents = 1
    END IF
    IF (numberGlobalZElements==0) THEN
      numberOfXi = 2
    ELSE
      numberOfXi = 3
    END IF

    IF (useSimplex) THEN
      interpolationType = 7
      quadratureOrder   = 3
      IF (useGeneratedMesh) CALL HANDLE_ERROR("Generated simplex mesh not set up.")
    ELSE
      interpolationType = 1
    END IF

    ! Get the number of computational nodes and this computational node number
    CALL cmfe_ComputationEnvironment_Initialise(computationEnvironment,err)
    CALL cmfe_ComputationEnvironment_NumberOfWorldNodesGet(computationEnvironment,numberOfComputationalNodes,err)
    CALL cmfe_ComputationEnvironment_WorldNodeNumberGet(computationEnvironment,computationalNodeNumber,err)

    CALL cmfe_CoordinateSystem_Initialise(coordinateSystem,Err)
    CALL cmfe_CoordinateSystem_CreateStart(coordinateSystemUserNumber,coordinateSystem,Err)
    CALL cmfe_CoordinateSystem_DimensionSet(coordinateSystem,3,Err)
    CALL cmfe_CoordinateSystem_CreateFinish(coordinateSystem,Err)

    ! Create a region and assign the coordinate system to the region
    CALL cmfe_Region_Initialise(region,Err)
    CALL cmfe_Region_CreateStart(regionUserNumber,worldRegion,region,Err)
    CALL cmfe_Region_LabelSet(region,"Region",Err)
    CALL cmfe_Region_CoordinateSystemSet(region,coordinateSystem,Err)
    CALL cmfe_Region_CreateFinish(region,Err)

    ! Define basis
    CALL cmfe_Basis_Initialise(basis,Err)
    CALL cmfe_Basis_CreateStart(basisUserNumber,basis,Err)
    SELECT CASE (interpolationType)
    CASE(1,2,3,4)
      CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
      CALL cmfe_Basis_NumberOfXiSet(basis,numberOfXi,Err)
      IF(numberGlobalZElements==0) THEN
        CALL cmfe_Basis_InterpolationXiSet(basis, &
          & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
          &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      ELSE
        CALL cmfe_Basis_InterpolationXiSet(basis, &
           & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
      END IF
      IF (numberOfGaussXi>0) THEN
        IF(numberGlobalZElements==0) THEN
          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi],Err)
        ELSE
          CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(basis,[numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],Err)
        END IF
      END IF
    CASE(7,8,9)
      CALL cmfe_Basis_TypeSet(basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
      CALL cmfe_Basis_NumberOfXiSet(basis,numberOfXi,Err)
      IF(numberGlobalZElements==0) THEN
        CALL cmfe_Basis_InterpolationXiSet(basis, &
          & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
          &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
      ELSE
        CALL cmfe_Basis_InterpolationXiSet(basis, &
           & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
           &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
      END IF
      CALL cmfe_Basis_QuadratureOrderSet(basis,quadratureOrder,Err)
    CASE DEFAULT
      CALL HANDLE_ERROR("Invalid interpolation type.")
    END SELECT
    CALL cmfe_Basis_CreateFinish(basis,Err)

    IF (usePressureBasis) THEN
      ! Define pressure basis
      CALL cmfe_Basis_Initialise(pressureBasis,Err)
      CALL cmfe_Basis_CreateStart(pressureBasisUserNumber,pressureBasis,Err)
      SELECT CASE (interpolationType)
      CASE(1,2,3,4)
        CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_LAGRANGE_HERMITE_TP_TYPE,Err)
        CALL cmfe_Basis_NumberOfXiSet(pressureBasis,numberOfXi,Err)
        IF(numberGlobalZElements==0) THEN
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
            & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
            &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
        ELSE
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
             & [CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION],Err)
        ENDIF
        IF (numberOfGaussXi>0) THEN
          IF(numberGlobalZElements==0) THEN
            CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi],Err)
          ELSE
            CALL cmfe_Basis_QuadratureNumberOfGaussXiSet(pressureBasis, &
              & [numberOfGaussXi,numberOfGaussXi,numberOfGaussXi],Err)
          END IF
        END IF
      CASE(7,8,9)
        CALL cmfe_Basis_TypeSet(pressureBasis,CMFE_BASIS_SIMPLEX_TYPE,Err)
        CALL cmfe_Basis_NumberOfXiSet(pressureBasis,numberOfXi,Err)
        IF(numberGlobalZElements==0) THEN
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
            & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
            &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
        ELSE
          CALL cmfe_Basis_InterpolationXiSet(pressureBasis, &
             & [CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
             &  CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION],Err)
        ENDIF
      CALL cmfe_Basis_QuadratureOrderSet(pressureBasis,quadratureOrder,Err)
      CASE DEFAULT
        CALL HANDLE_ERROR("Invalid interpolation type.")
      END SELECT
      CALL cmfe_Basis_CreateFinish(pressureBasis,Err)
    END IF

    CALL cmfe_Mesh_Initialise(Mesh,Err)
    IF (useGeneratedMesh) THEN
      ! Start the creation of a generated mesh in the region
      CALL cmfe_GeneratedMesh_Initialise(generatedMesh,Err)
      CALL cmfe_GeneratedMesh_CreateStart(generatedMeshUserNumber,region,generatedMesh,Err)
      CALL cmfe_GeneratedMesh_TypeSet(generatedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
      IF (usePressureBasis) THEN
        CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,[basis,pressureBasis],Err)
      ELSE
        CALL cmfe_GeneratedMesh_BasisSet(generatedMesh,basis,Err)
        ! TODO those should not be in else, right?
        CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[width,length,height],Err)
        CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh, &
         & [numberGlobalXElements,numberGlobalYElements, &
         & numberGlobalZElements],Err)
        ! end TODO
      END IF
      CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
    ELSE
      ! Start the creation of a manually generated mesh in the region
      CALL cmfe_Mesh_CreateStart(meshUserNumber,region,numberOfXi,mesh,Err)
      CALL cmfe_Mesh_NumberOfComponentsSet(mesh,numberOfMeshComponents,Err)
      IF (useSimplex) THEN
        CALL cmfe_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements*5,Err)
      ELSE
        CALL cmfe_Mesh_NumberOfElementsSet(mesh,totalNumberOfElements,Err)
      END IF

      ! Define nodes for the mesh
      CALL cmfe_Nodes_Initialise(nodes,Err)
      CALL cmfe_Nodes_CreateStart(region,totalNumberOfNodes,nodes,Err)
      CALL cmfe_Nodes_CreateFinish(nodes,Err)

      CALL cmfe_MeshElements_Initialise(elements,Err)
      CALL cmfe_MeshElements_CreateStart(mesh,meshComponentNumber,basis,elements,Err)
      IF (useSimplex) THEN
        CALL cmfe_MeshElements_NodesSet(elements,1,[1,2,4,6],Err)
        CALL cmfe_MeshElements_NodesSet(elements,2,[1,4,3,7],Err)
        CALL cmfe_MeshElements_NodesSet(elements,3,[1,6,7,5],Err)
        CALL cmfe_MeshElements_NodesSet(elements,4,[6,4,7,8],Err)
        CALL cmfe_MeshElements_NodesSet(elements,5,[1,6,4,7],Err)
      ELSE
        CALL cmfe_MeshElements_NodesSet(elements,1,[1,2,3,4,5,6,7,8],Err)
      END IF
      CALL cmfe_MeshElements_CreateFinish(elements,Err)

      CALL cmfe_Mesh_CreateFinish(mesh,Err)
    END IF

    ! Create a decomposition for the mesh
    CALL cmfe_Decomposition_Initialise(decomposition,Err)
    CALL cmfe_Decomposition_CreateStart(decompositionUserNumber,mesh,decomposition,Err)
    CALL cmfe_Decomposition_TypeSet(decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
    CALL cmfe_Decomposition_NumberOfDomainsSet(decomposition,numberOfComputationalNodes,Err)
    CALL cmfe_Decomposition_CreateFinish(decomposition,Err)

    ! Create a field for the geometry
    CALL cmfe_Field_Initialise(geometricField,Err)
    CALL cmfe_Field_CreateStart(geometricFieldUserNumber,region,geometricField,Err)
    CALL cmfe_Field_MeshDecompositionSet(geometricField,Decomposition,Err)
    CALL cmfe_Field_TypeSet(geometricField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
    CALL cmfe_Field_VariableLabelSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,"Geometry",Err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(geometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
    IF (interpolationType==4) THEN
      CALL cmfe_Field_ScalingTypeSet(geometricField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    END IF
    CALL cmfe_Field_CreateFinish(geometricField,Err)

    IF (useGeneratedMesh) THEN
      ! Update the geometric field parameters from generated mesh
      CALL cmfe_GeneratedMesh_GeometricParametersCalculate(generatedMesh,geometricField,Err)
    ELSE
      ! Update the geometric field parameters manually
      CALL cmfe_Field_ParameterSetUpdateStart(geometricField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
      ! node 1
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,1,3,0.0_CMISSRP,Err)
      ! node 2
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,1,height,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,2,3,0.0_CMISSRP,Err)
      ! node 3
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,2,width,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,3,3,0.0_CMISSRP,Err)
      ! node 4
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,1,height,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,2,width,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,4,3,0.0_CMISSRP,Err)
      ! node 5
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,5,3,length,Err)
      ! node 6
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,1,height,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,2,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,6,3,length,Err)
      ! node 7
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,1,0.0_CMISSRP,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,2,width,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,7,3,length,Err)
      ! node 8
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,1,height,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,2,width,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(geometricField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & CMFE_FIELD_VALUES_SET_TYPE,1,1,8,3,length,Err)
      CALL cmfe_Field_ParameterSetUpdateFinish(geometricField, &
        & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
    END IF

    ! Create a fibre field and attach it to the geometric field
    CALL cmfe_Field_Initialise(fibreField,Err)
    CALL cmfe_Field_CreateStart(fibreFieldUserNumber,region,fibreField,Err)
    CALL cmfe_Field_TypeSet(fibreField,CMFE_FIELD_FIBRE_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(fibreField,decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(fibreField,geometricField,Err)
    CALL cmfe_Field_VariableLabelSet(fibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
    IF (interpolationType==4) THEN
      CALL cmfe_Field_ScalingTypeSet(fibreField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    END IF
    CALL cmfe_Field_CreateFinish(fibreField,Err)

    ! Create the material field
    IF (compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 2
    END IF
    CALL cmfe_Field_Initialise(materialField,Err)
    CALL cmfe_Field_CreateStart(materialFieldUserNumber,Region,materialField,Err)
    CALL cmfe_Field_TypeSet(materialField,CMFE_FIELD_MATERIAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(materialField,decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(materialField,geometricField,Err)
    CALL cmfe_Field_NumberOfVariablesSet(materialField,1,Err)
    CALL cmfe_Field_NumberOfComponentsSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,Err)
    CALL cmfe_Field_VariableLabelSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,"Material",Err)
    CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
    IF (compressible) THEN
      CALL cmfe_Field_ComponentMeshComponentSet(materialField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
    END IF
    IF (interpolationType==4) THEN
      CALL cmfe_Field_ScalingTypeSet(materialField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    END IF
    CALL cmfe_Field_CreateFinish(materialField,Err)

    ! Set Mooney-Rivlin constants c10 and c01 respectively.
    CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 1,2.0_CMISSRP,Err)
    CALL cmfe_Field_ComponentValuesInitialise(materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 2,6.0_CMISSRP,Err)
    IF (compressible) THEN
      CALL cmfe_Field_ComponentValuesInitialise( &
        & materialField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 3,1.0e9_CMISSRP,Err)
    END IF

    ! Create the dependent field
    IF (compressible) THEN
      numberOfMaterialComponents = 3
    ELSE
      numberOfMaterialComponents = 4
    END IF

    CALL cmfe_Field_Initialise(dependentField,Err)
    CALL cmfe_Field_CreateStart(dependentFieldUserNumber,region,dependentField,Err)
    CALL cmfe_Field_VariableLabelSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent",Err)
    CALL cmfe_Field_TypeSet(dependentField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
    CALL cmfe_Field_MeshDecompositionSet(dependentField,decomposition,Err)
    CALL cmfe_Field_GeometricFieldSet(dependentField,geometricField,Err)
    CALL cmfe_Field_DependentTypeSet(dependentField,CMFE_FIELD_DEPENDENT_TYPE,Err)
    CALL cmfe_Field_NumberOfVariablesSet(dependentField,2,Err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,numberOfMaterialComponents,Err)
    CALL cmfe_Field_NumberOfComponentsSet(dependentField, &
      & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,numberOfMaterialComponents,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,1,Err)
    CALL cmfe_Field_ComponentMeshComponentSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,1,Err)
    IF (.NOT.compressible) THEN
      ! TODO we always had node based interpolation, right? --> check!
      CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
        & 4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
      CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE, &
        & 4,CMFE_FIELD_ELEMENT_BASED_INTERPOLATION,Err)
      ! TODO end
      IF (usePressureBasis) THEN
        ! Set the pressure to be nodally based and use the second mesh component
        IF (interpolationType==4) THEN
          CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4, &
            & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
          CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4, &
            & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
        END IF
        CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_U_VARIABLE_TYPE,4,2,Err)
        CALL cmfe_Field_ComponentInterpolationSet(dependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,4,2,Err)
      END IF
    END IF
    IF (interpolationType==4) CALL cmfe_Field_ScalingTypeSet(dependentField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    CALL cmfe_Field_CreateFinish(dependentField,Err)

    ! Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
    CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
      & geometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3, &
      & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,Err)
    IF (.NOT.compressible) THEN
      CALL cmfe_Field_ComponentValuesInitialise( &
        & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
        & 4,-8.0_CMISSRP,Err)
    END IF

    ! Create a deformed geometry field, as cmgui doesn't like displaying
    ! deformed fibres from the dependent field because it isn't a geometric field.
    CALL cmfe_Field_Initialise(deformedField,Err)
    CALL cmfe_Field_CreateStart(deformedFieldUserNumber,region,deformedField,Err)
    CALL cmfe_Field_MeshDecompositionSet(deformedField,decomposition,Err)
    CALL cmfe_Field_TypeSet(deformedField,CMFE_FIELD_GEOMETRIC_TYPE,Err)
    CALL cmfe_Field_VariableLabelSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,"DeformedGeometry",Err)
    DO componentIdx=1,3
      CALL cmfe_Field_ComponentMeshComponentSet(deformedField,CMFE_FIELD_U_VARIABLE_TYPE,componentIdx,1,Err)
    END DO
    IF (interpolationType==4) CALL cmfe_Field_ScalingTypeSet(deformedField,CMFE_FIELD_ARITHMETIC_MEAN_SCALING,Err)
    CALL cmfe_Field_CreateFinish(deformedField,Err)

    ! Create the equations_set
    CALL cmfe_Field_Initialise(equationsSetField,Err)
    CALL cmfe_EquationsSet_Initialise(equationsSet,Err)
    IF (compressible) THEN
      CALL cmfe_EquationsSet_CreateStart(equationsSetUserNumber,region,fibreField, &
        & [CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  CMFE_EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE], &
        & equationsSetFieldUserNumber,equationsSetField,equationsSet,Err)
    ELSE
      CALL cmfe_EquationsSet_CreateStart(equationsSetUserNumber,region,fibreField, &
        & [CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
        &  CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE, &
        &  CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE], &
        & equationsSetFieldUserNumber,equationsSetField,equationsSet,Err)
    END IF
    CALL cmfe_EquationsSet_CreateFinish(equationsSet,Err)
    CALL cmfe_EquationsSet_MaterialsCreateStart(equationsSet,materialFieldUserNumber,materialField,Err)
    CALL cmfe_EquationsSet_MaterialsCreateFinish(equationsSet,Err)
    CALL cmfe_EquationsSet_DependentCreateStart(equationsSet,dependentFieldUserNumber,dependentField,Err)
    CALL cmfe_EquationsSet_DependentCreateFinish(equationsSet,Err)

    ! Create equations
    CALL cmfe_Equations_Initialise(equations,Err)
    CALL cmfe_EquationsSet_EquationsCreateStart(equationsSet,equations,Err)
    CALL cmfe_Equations_SparsityTypeSet(equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
    CALL cmfe_Equations_OutputTypeSet(equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
    CALL cmfe_EquationsSet_EquationsCreateFinish(equationsSet,Err)

    ! Define the problem
    CALL cmfe_Problem_Initialise(problem,Err)
    CALL cmfe_Problem_CreateStart(problemUserNumber, &
      & [CMFE_PROBLEM_ELASTICITY_CLASS, &
      &  CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
      &  CMFE_PROBLEM_NO_SUBTYPE],problem,Err)
    CALL cmfe_Problem_CreateFinish(problem,Err)

    ! Create control loops
    CALL cmfe_Problem_ControlLoopCreateStart(problem,Err)
    CALL cmfe_Problem_ControlLoopCreateFinish(problem,Err)

    ! Create problem solver
    CALL cmfe_Solver_Initialise(nonlinearSolver,Err)
    CALL cmfe_Solver_Initialise(linearSolver,Err)
    CALL cmfe_Problem_SolversCreateStart(problem,Err)
    CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,nonLinearSolver,Err)
    CALL cmfe_Solver_OutputTypeSet(nonlinearSolver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(nonlinearSolver, &
      & CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
    CALL cmfe_Solver_NewtonLinearSolverGet(nonlinearSolver,linearSolver,Err)
    CALL cmfe_Solver_NewtonAbsoluteToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,Err)
    CALL cmfe_Solver_NewtonSolutionToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,Err)
    CALL cmfe_Solver_NewtonRelativeToleranceSet(nonlinearSolver,1.0E-14_CMISSRP,Err)
    CALL cmfe_Solver_LinearTypeSet(linearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
    CALL cmfe_Problem_SolversCreateFinish(problem,Err)

    ! Create solver equations and add equations set to solver equations
    CALL cmfe_Solver_Initialise(solver,Err)
    CALL cmfe_SolverEquations_Initialise(solverEquations,Err)
    CALL cmfe_Problem_SolverEquationsCreateStart(problem,Err)
    CALL cmfe_Problem_SolverGet(problem,CMFE_CONTROL_LOOP_NODE,1,solver,Err)
    CALL cmfe_Solver_SolverEquationsGet(solver,solverEquations,Err)
    CALL cmfe_SolverEquations_SparsityTypeSet(solverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
    CALL cmfe_SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,Err)
    CALL cmfe_Problem_SolverEquationsCreateFinish(problem,Err)


    ! Prescribe boundary conditions (absolute nodal parameters)
    CALL cmfe_BoundaryConditions_Initialise(boundaryConditions,Err)
    CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,Err)

    ! Set x=0 nodes to no x displacment in x. Set x=width nodes to 10% x displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,7,1,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    IF (zeroLoad) THEN
      load = 0.0_CMISSRP
    ELSE
      load = 0.1_CMISSRP*width
    END IF
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,1,CMFE_BOUNDARY_CONDITION_FIXED,load,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,1,CMFE_BOUNDARY_CONDITION_FIXED,load,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,1,CMFE_BOUNDARY_CONDITION_FIXED,load,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,8,1,CMFE_BOUNDARY_CONDITION_FIXED,load,Err)

    ! Set y=0 nodes to no y displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,5,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,6,2,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)

    ! Set z=0 nodes to no y displacement
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,1,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,2,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,3,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(boundaryConditions,dependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & 1,1,4,3,CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)

    CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(solverEquations,Err)

    ! Solve the problem
    CALL cmfe_Problem_Solve(problem,Err)

    ! Copy deformed geometry into deformed field
    DO componentIdx=1,3
      CALL cmfe_Field_ParametersToFieldParametersComponentCopy( &
        & dependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx, &
        & deformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,componentIdx,Err)
    END DO

    IF (useGeneratedMesh) THEN
      output_file = "./results/unit_cube_generated_mesh"
    ELSE
      IF (useSimplex) THEN
        output_file = "./results/unit_cube_manual_mesh_simplex"
      ELSE
        output_file = "./results/unit_cube_manual_mesh"
      END IF
    END IF

    INQUIRE(file="./results", exist=directory_exists)
    IF (.NOT.directory_exists) THEN
      CALL execute_command_line ("mkdir ./results")
    END IF

    prefix = ''
    IF (compressible) THEN
      IF (zeroLoad) THEN
        prefix = "_compressible_zero_load"
      ELSE
        prefix = "_compressible"
      END IF
    ELSE
      IF (zeroLoad) THEN
        prefix = "_zero_load"
      END IF
    END IF

    ! Export results
    CALL cmfe_Fields_Initialise(fields,Err)
    CALL cmfe_Fields_Create(region,fields,Err)
    CALL cmfe_Fields_NodesExport(fields,trim(output_file)//trim(prefix),"FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(fields,trim(output_file)//trim(prefix),"FORTRAN",Err)
    CALL cmfe_Fields_Finalise(fields,Err)

    CALL cmfe_Problem_Destroy(ProblemUserNumber,Err)
    IF (useGeneratedMesh) THEN
      CALL cmfe_GeneratedMesh_Destroy(RegionUserNumber,GeneratedMeshUserNumber,Err)
    END IF
    CALL cmfe_Basis_Destroy(BasisUserNumber,Err)
    CALL cmfe_Region_Destroy(RegionUserNumber,Err)
    CALL cmfe_CoordinateSystem_Destroy(CoordinateSystemUserNumber,Err)

    WRITE(*,'(A)') "Program successfully completed."

  END SUBROUTINE SOLVE_MODEL

  SUBROUTINE HANDLE_ERROR(ERROR_STRING)
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING

    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP
  END SUBROUTINE HANDLE_ERROR

END PROGRAM UNIAXIAL_EXTENSION
