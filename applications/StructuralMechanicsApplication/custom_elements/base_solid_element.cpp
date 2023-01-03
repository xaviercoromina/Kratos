// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_elements/base_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

void BaseSolidElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if(this->UseGeometryIntegrationMethod()) {
            if( GetProperties().Has(INTEGRATION_ORDER) ) {
                const SizeType integration_order = GetProperties()[INTEGRATION_ORDER];
                switch ( integration_order )
                {
                case 1:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
                    break;
                case 2:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
                    break;
                case 3:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
                    break;
                case 4:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_4;
                    break;
                case 5:
                    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
                    break;
                default:
                    KRATOS_WARNING("BaseSolidElement") << "Integration order " << integration_order << " is not available, using default integration order for the geometry" << std::endl;
                    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
                }
            } else {
                mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
            }
        }

        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
            mConstitutiveLawVector.resize(integration_points.size());

        InitializeMaterial();

    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const bool is_rotated = IsElementRotated();
        const auto& r_geom = GetGeometry();
        const SizeType number_of_nodes = r_geom.size();
        const SizeType dimension = r_geom.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
        const Properties& r_properties = GetProperties();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geom,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points

        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints();

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

            // Rotate to local axes strain/F
            if (is_rotated)
                RotateToLocalAxes(Values, this_kinematic_variables);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->InitializeMaterialResponse(Values, GetStressMeasure());

            // TODO: Deprecated, remove this
            mConstitutiveLawVector[point_number]->InitializeSolutionStep( r_properties, r_geom, row( N_values, point_number ), rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->InitializeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeNonLinearIteration( const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        // TODO: Deprecated, remove this
        mConstitutiveLawVector[point_number]->FinalizeNonLinearIteration( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::FinalizeSolutionStep( const ProcessInfo& rCurrentProcessInfo )
{
    bool required = false;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const bool is_rotated = IsElementRotated();
        const auto &r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geometry,r_properties,rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions=Values.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        Values.SetStrainVector(this_constitutive_variables.StrainVector);
        Values.SetStressVector(this_constitutive_variables.StressVector);
        Values.SetConstitutiveMatrix(this_constitutive_variables.D);

        // Reading integration points
        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);

        // Reading integration points
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

            // Rotate to local axes strain/F
            if (is_rotated)
                RotateToLocalAxes(Values, this_kinematic_variables);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, GetStressMeasure());

            // TODO: Deprecated, remove this
            mConstitutiveLawVector[point_number]->FinalizeSolutionStep( r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure BaseSolidElement::GetStressMeasure() const
{
    return ConstitutiveLaw::StressMeasure_PK2;
}

/***********************************************************************************/
/***********************************************************************************/

bool BaseSolidElement::UseElementProvidedStrain() const
{
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& N_values = this->ShapeFunctionsValues(mThisIntegrationMethod);
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            mConstitutiveLawVector[point_number]->ResetMaterial( r_properties,  r_geometry, row( N_values, point_number ) );
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer BaseSolidElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    KRATOS_WARNING("BaseSolidElement") << " Call BaseSolidElement (base class) Clone " << std::endl;

    BaseSolidElement::Pointer p_new_elem = Kratos::make_intrusive<BaseSolidElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = r_geom[0].GetDofPosition(DISPLACEMENT_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = r_geom[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geom[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = r_geom[i].GetDof(DISPLACEMENT_X,pos).EquationId();
            rResult[index + 1] = r_geom[i].GetDof(DISPLACEMENT_Y,pos+1).EquationId();
            rResult[index + 2] = r_geom[i].GetDof(DISPLACEMENT_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension * number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_geom[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("")
};

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValuesVector(
    Vector& rValues,
    int Step
    ) const
{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& displacement = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k) {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = r_geom[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = velocity[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    ) const
{
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = r_geom[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = acceleration[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector, rCurrentProcessInfo);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * dimension;

            AtomicAdd(r_geom[i].GetValue(NODAL_MASS), element_mass_vector[index]);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector damping_residual_contribution = ZeroVector(element_size);

    // Calculate damping contribution to residual -->
    if (r_prop.Has(RAYLEIGH_ALPHA) || r_prop.Has(RAYLEIGH_BETA)) {
        Vector current_nodal_velocities = ZeroVector(element_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);

        Matrix damping_matrix(element_size, element_size);
        this->CalculateDampingMatrixWithLumpedMass(damping_matrix, rCurrentProcessInfo);

        // Current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < dimension; ++j) {
                AtomicAdd(r_force_residual[j], (rRHSVector[index + j] - damping_residual_contribution[index + j]));
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    //calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                             const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = false;
    VectorType RHS;

    CalculateAll( rLeftHandSideMatrix, RHS, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Calculation flags
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    noalias(rMassMatrix) = ZeroMatrix( mat_size, mat_size );

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Checking if computing lumped mass matrix
    const bool compute_lumped_mass_matrix = StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(r_prop, rCurrentProcessInfo);

    // LUMPED MASS MATRIX
    if (compute_lumped_mass_matrix) {
        VectorType temp_vector(mat_size);
        this->CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);
        for (IndexType i = 0; i < mat_size; ++i)
            rMassMatrix(i, i) = temp_vector[i];
    } else { // CONSISTENT MASS
        const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
        const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

        Matrix J0(dimension, dimension);

        IntegrationMethod integration_method = UseGeometryIntegrationMethod() ? IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom) : mThisIntegrationMethod ;
        const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints( integration_method );
        const Matrix& Ncontainer = this->ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
            GeometryUtils::JacobianOnInitialConfiguration(
                r_geom, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::Det(J0);
            const double integration_weight =
                GetIntegrationWeight(integration_points, point_number, detJ0) * thickness;
            const Vector& rN = row(Ncontainer,point_number);

            for ( IndexType i = 0; i < number_of_nodes; ++i ) {
                const SizeType index_i = i * dimension;

                for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                    const SizeType index_j = j * dimension;
                    const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                    for ( IndexType k = 0; k < dimension; ++k )
                        rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const unsigned int mat_size = GetGeometry().PointsNumber() * GetGeometry().WorkingSpaceDimension();

    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        mat_size);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points, false);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number <number_of_integration_points; ++point_number ) {
            bool value;
            mConstitutiveLawVector[point_number]->GetValue( rVariable, value);
            rOutput[point_number] = value;
        }
    } else {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ++ii ) {
            bool solution;
            solution = mConstitutiveLawVector[ii]->CalculateValue( Values, rVariable, solution);
            rOutput[ii] = solution;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points, false);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    const std::size_t number_of_integration_points = integration_points.size();
    const auto& r_geometry = GetGeometry();
    const auto& r_properties = GetProperties();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType number_of_nodes = r_geometry.size();

    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points, false );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_WEIGHT) {

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            double integration_weight;

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                this_kinematic_variables.detJ0 = CalculateDerivativesOnReferenceConfiguration(this_kinematic_variables.J0,
                                                                                    this_kinematic_variables.InvJ0,
                                                                                    this_kinematic_variables.DN_DX,
                                                                                    point_number,
                                                                                    this->GetIntegrationMethod());

                integration_weight = GetIntegrationWeight(integration_points,
                                                                    point_number,
                                                                    this_kinematic_variables.detJ0);

                if (dimension == 2 && r_properties.Has(THICKNESS))
                    integration_weight *= r_properties[THICKNESS];

                rOutput[point_number] = integration_weight;
            }
        } else if ( rVariable == STRAIN_ENERGY ) {

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,r_properties,rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                double StrainEnergy = 0.0;

                mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);

                rOutput[point_number] = StrainEnergy;
            }
        } else if ( rVariable == ERROR_INTEGRATION_POINT ) {

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,r_properties,rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            //Calculate Cauchy Stresses from the FE solution
            std::vector<Vector> sigma_FE_solution(number_of_nodes);
            const Variable<Vector>& r_variable_stress = CAUCHY_STRESS_VECTOR;
            CalculateOnIntegrationPoints(r_variable_stress, sigma_FE_solution, rCurrentProcessInfo);

            // calculate the determinatn of the Jacobian in the current configuration
            Vector detJ(number_of_integration_points);
            if (UseGeometryIntegrationMethod()){
                detJ = r_geometry.DeterminantOfJacobian(detJ);
            } else {
                for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                   detJ[point_number] = r_geometry.DeterminantOfJacobian(integration_points[point_number]);
                }
            }

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            if (r_geometry[0].Has(RECOVERED_STRESS)) {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    // Compute element kinematics B, F, DN_DX ...
                    CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), false);

                    double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ[point_number]);

                    if (dimension == 2 && r_properties.Has(THICKNESS))
                        integration_weight *= r_properties[THICKNESS];

                    // Calculate recovered stresses at integration points
                    Vector sigma_recovered = ZeroVector(strain_size);

                    // sigma_recovered = sum(N_i * sigma_recovered_i)
                    for (IndexType node_number=0; node_number<number_of_nodes; node_number++) {
                        const auto& r_sigma_recovered_node = r_geometry[node_number].GetValue(RECOVERED_STRESS);
                        for (IndexType stress_component = 0; stress_component<strain_size; stress_component++) {
                            sigma_recovered[stress_component] += this_kinematic_variables.N[node_number] * r_sigma_recovered_node[stress_component];
                        }
                    }

                    // Calculate error_sigma
                    Vector error_sigma(strain_size);
                    error_sigma = sigma_recovered - sigma_FE_solution[point_number];

                    // For debug
                    KRATOS_TRACE("ERROR_INTEGRATION_POINT")
                    <<"sigma recovered: " << sigma_recovered << std::endl
                    <<"sigma FE: " << sigma_FE_solution[point_number] << std::endl;

                    // Calculate inverse of material matrix
                    Matrix invD(strain_size,strain_size);
                    double detD;
                    MathUtils<double>::InvertMatrix(this_constitutive_variables.D, invD,detD);

                    // Calculate error_energy
                    rOutput[point_number] = integration_weight * inner_prod(error_sigma, prod(invD, error_sigma));
                }
            } else {
                for (IndexType point_number = 0; point_number < number_of_integration_points; point_number++) {
                    rOutput[point_number] = 0.0;
                }
            }
        } else if (rVariable == VON_MISES_STRESS) {

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,r_properties,rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());
                // Compute material reponse, not encessary to rotate since it's an invariant
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), false);

                // Compute VM stress
                if (dimension == 2) {
                    if (strain_size == 3) {
                        rOutput[point_number] = ConstitutiveLawUtilities<3>::CalculateVonMisesEquivalentStress(this_constitutive_variables.StressVector);
                    } else { // Axysimmetric 4
                        Vector aux_stress(6);
                        noalias(aux_stress) = ZeroVector(6);
                        for (IndexType i = 0; i < 4; ++i)
                            aux_stress(i) = this_constitutive_variables.StressVector(i);
                        rOutput[point_number] = ConstitutiveLawUtilities<6>::CalculateVonMisesEquivalentStress(aux_stress);
                    }
                } else { // 3D
                    rOutput[point_number] = ConstitutiveLawUtilities<6>::CalculateVonMisesEquivalentStress(this_constitutive_variables.StressVector);
                }
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_COORDINATES) {

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                Point global_point;
                r_geom.GlobalCoordinates(global_point, integration_points[point_number]);

                noalias(rOutput[point_number]) = global_point.Coordinates();
            }
        } else if (rVariable == LOCAL_AXIS_1 || rVariable == LOCAL_AXIS_2 || rVariable == LOCAL_AXIS_3) {
            if (this->Has(rVariable)) {
                for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
                    noalias(rOutput[point_number]) = this->GetValue(rVariable);
            } else if (rVariable == LOCAL_AXIS_3) {
                const array_1d<double, 3> r_local_axis_1 = this->GetValue(LOCAL_AXIS_1);
                const array_1d<double, 3> local_axis_2 = this->GetValue(LOCAL_AXIS_2);
                for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number)
                    noalias(rOutput[point_number]) = MathUtils<double>::CrossProduct(r_local_axis_1, local_axis_2);
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }  else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints( this->GetIntegrationMethod() );
    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == INSITU_STRESS ) {
            Vector strain_vector( strain_size );

            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size() != strain_vector.size() )
                    rOutput[point_number].resize( strain_vector.size(), false );

                rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
            }
        } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
            const bool is_rotated = IsElementRotated();
            // Create and initialize element variables:

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geom,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());
                //call the constitutive law to update material variables
                if( rVariable == CAUCHY_STRESS_VECTOR) {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy, is_rotated);
                } else {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points,ConstitutiveLaw::StressMeasure_PK2, is_rotated);
                }

                if (strain_size == 4) { // Axysimmetric
                    if (rOutput[point_number].size() != 6)
                        rOutput[point_number].resize(6, false);
                    noalias(rOutput[point_number]) = ZeroVector(6);
                    for (IndexType i = 0; i < 4; ++i)
                        rOutput[point_number](i) = this_constitutive_variables.StressVector(i);
                } else {
                    if (rOutput[point_number].size() != strain_size)
                        rOutput[point_number].resize(strain_size, false);

                    noalias(rOutput[point_number]) = this_constitutive_variables.StressVector;
                }

            }
        } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR ) {
            // Create and initialize element variables:

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geom,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags &ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

            //reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());
                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this_stress_measure, false);

                if (strain_size == 4) { // Axysimmetric
                    if (rOutput[point_number].size() != 6)
                        rOutput[point_number].resize(6, false);
                    noalias(rOutput[point_number]) = ZeroVector(6);
                    for (IndexType i = 0; i < 4; ++i)
                        rOutput[point_number](i) = this_constitutive_variables.StrainVector(i);
                } else {
                    if (rOutput[point_number].size() != strain_size)
                        rOutput[point_number].resize(strain_size, false);

                    noalias(rOutput[point_number]) = this_constitutive_variables.StrainVector;
                }
            }
        } else if (rVariable == INITIAL_STRESS_VECTOR) {
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                if (mConstitutiveLawVector[point_number]->HasInitialState()) {
                    const Vector& r_initial_stress = mConstitutiveLawVector[point_number]->GetInitialState().GetInitialStressVector();

                    if ( rOutput[point_number].size() != strain_size)
                        rOutput[point_number].resize( strain_size, false );

                    noalias(rOutput[point_number]) = r_initial_stress;
                } else {
                    noalias(rOutput[point_number]) = ZeroVector(strain_size);
                }
            }
        } else if (rVariable == INITIAL_STRAIN_VECTOR) {
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                if (mConstitutiveLawVector[point_number]->HasInitialState()) {
                    const Vector& r_initial_strain = mConstitutiveLawVector[point_number]->GetInitialState().GetInitialStrainVector();

                    if ( rOutput[point_number].size() != strain_size)
                        rOutput[point_number].resize( strain_size, false );

                    noalias(rOutput[point_number]) = r_initial_strain;
                } else {
                    noalias(rOutput[point_number]) = ZeroVector(strain_size);
                }
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
        

    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints( this->GetIntegrationMethod() );
    const auto& r_geom = GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
            std::vector<Vector> stress_vector;

            if( rVariable == CAUCHY_STRESS_TENSOR )
                this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
            else
                this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                noalias(rOutput[point_number]) = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
            }
        }
        else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR) {
            std::vector<Vector> strain_vector;
            if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
                CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );
            else
                CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != dimension )
                    rOutput[point_number].resize( dimension, dimension, false );

                noalias(rOutput[point_number]) = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
            }
        } else if ( rVariable == CONSTITUTIVE_MATRIX ) {
            const bool is_rotated = IsElementRotated();
            // Create and initialize element variables:
            const SizeType number_of_nodes = r_geom.size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geom,GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);
            Values.SetConstitutiveMatrix(this_constitutive_variables.D); //this is the output parameter

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), is_rotated);

                if( rOutput[point_number].size2() != this_constitutive_variables.D.size2() )
                    rOutput[point_number].resize( this_constitutive_variables.D.size1() , this_constitutive_variables.D.size2() , false );

                noalias(rOutput[point_number]) = this_constitutive_variables.D;
            }
        } else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES
            // Create and initialize element variables:
            const SizeType number_of_nodes = r_geom.size();
            const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

            KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
            ConstitutiveVariables this_constitutive_variables(strain_size);

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geom,GetProperties(),rCurrentProcessInfo);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

                if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                    rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

                noalias(rOutput[point_number]) = this_kinematic_variables.F;
            }
        }  else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<bool>& rVariable,
    const std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    const std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    const std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        for ( IndexType point_number = 0; point_number < integration_points_number; ++point_number ) {
            mConstitutiveLawVector[point_number] = rValues[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    const std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    const std::vector<array_1d<double, 6>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetValuesOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    const std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("BaseSolidElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

int  BaseSolidElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    // Basic check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_ERROR << "You have called to the CalculateAll from the base class for solid elements" << std::endl;
}

//***********************************************************************
//***********************************************************************

double BaseSolidElement::GetIntegrationWeight(
    const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
    const IndexType PointNumber,
    const double detJ
    ) const
{
    return rThisIntegrationPoints[PointNumber].Weight() * detJ;
}

//***********************************************************************
//***********************************************************************

void BaseSolidElement::CalculateShapeGradientOfMassMatrix(MatrixType& rMassMatrix, ShapeParameter Deriv) const
{
    KRATOS_TRY;

    // Properties
    const auto& r_prop = GetProperties();

    // Geometry information
    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    // Getting density
    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    const IntegrationMethod integration_method = this->UseGeometryIntegrationMethod() ? IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geom) : mThisIntegrationMethod;
    const Matrix& Ncontainer = this->ShapeFunctionsValues(integration_method);
    Matrix J0(dimension, dimension), DN_DX0_deriv;
    const auto& integration_points = this->IntegrationPoints(integration_method);
    for (unsigned point_number = 0; point_number < integration_points.size(); ++point_number) {
        Matrix DN_De;
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        if(UseGeometryIntegrationMethod()) {
            DN_De = r_geom.ShapeFunctionsLocalGradients(integration_method)[point_number];
        } else {
            r_geom.ShapeFunctionsLocalGradients(DN_De, integration_points[point_number]);
        }
        GeometricalSensitivityUtility geometrical_sensitivity(J0, DN_De);
        double detJ0_deriv;
        geometrical_sensitivity.CalculateSensitivity(Deriv, detJ0_deriv, DN_DX0_deriv);
        const double integration_weight =
            GetIntegrationWeight(integration_points, point_number, detJ0_deriv) * thickness;
        const Vector& rN = row(Ncontainer, point_number);

        for (unsigned i = 0; i < r_geom.size(); ++i) {
            const unsigned index_i = i * dimension;

            for (unsigned j = 0; j < r_geom.size(); ++j) {
                const unsigned index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight * density;

                for (unsigned k = 0; k < dimension; ++k)
                    rMassMatrix(index_i + k, index_j + k) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    KRATOS_ERROR << "You have called to the CalculateKinematicVariables from the base class for solid elements" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure,
    const bool IsElementRotated
    )
{
    // Setting the variables for the CL
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // rotate to local axes strain/F
    if (IsElementRotated)
        RotateToLocalAxes(rValues, rThisKinematicVariables);

    // Actually do the computations in the ConstitutiveLaw in local axes
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done

    // We undo the rotation of strain/F, C, stress
    if (IsElementRotated)
        RotateToGlobalAxes(rValues, rThisKinematicVariables);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::BuildRotationSystem(
    BoundedMatrix<double, 3, 3>& rRotationMatrix,
    const SizeType StrainSize
    )
{
    const array_1d<double, 3>& r_local_axis_1 = this->GetValue(LOCAL_AXIS_1);
    array_1d<double, 3> local_axis_2;
    array_1d<double, 3> local_axis_3;

    if (StrainSize == 6) {
        noalias(local_axis_2) = this->GetValue(LOCAL_AXIS_2);
        noalias(local_axis_3) = MathUtils<double>::CrossProduct(r_local_axis_1, local_axis_2);
    } else if (StrainSize == 3) { // we assume xy plane
        local_axis_2[0] = r_local_axis_1[1];
        local_axis_2[1] = -r_local_axis_1[0];
        local_axis_2[2] = 0.0;
        local_axis_3[0] = 0.0;
        local_axis_3[1] = 0.0;
        local_axis_3[2] = 1.0;
    }
    StructuralMechanicsElementUtilities::InitialCheckLocalAxes(r_local_axis_1, local_axis_2, local_axis_3);
    StructuralMechanicsElementUtilities::BuildRotationMatrix(rRotationMatrix, r_local_axis_1, local_axis_2, local_axis_3);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::RotateToLocalAxes(
    ConstitutiveLaw::Parameters& rValues,
    KinematicVariables& rThisKinematicVariables
    )
{
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    BoundedMatrix<double, 3, 3> rotation_matrix;

    BuildRotationSystem(rotation_matrix, strain_size);

    if (UseElementProvidedStrain()) { // we rotate strain
        if (strain_size == 6) {
            BoundedMatrix<double, 6, 6> voigt_rotation_matrix;
            ConstitutiveLawUtilities<6>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
            rValues.GetStrainVector() = prod(voigt_rotation_matrix, rValues.GetStrainVector());
        } else if (strain_size == 3) {
            BoundedMatrix<double, 3, 3> voigt_rotation_matrix;
            ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
            rValues.GetStrainVector() = prod(voigt_rotation_matrix, rValues.GetStrainVector());
        }
    } else { // Rotate F
        BoundedMatrix<double, 3, 3> inv_rotation_matrix;
        double aux_det;
        MathUtils<double>::InvertMatrix3(rotation_matrix, inv_rotation_matrix, aux_det);
        rThisKinematicVariables.F = prod(rotation_matrix, rThisKinematicVariables.F);
        rThisKinematicVariables.F = prod(rThisKinematicVariables.F, inv_rotation_matrix);
        rValues.SetDeformationGradientF(rThisKinematicVariables.F);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::RotateToGlobalAxes(
    ConstitutiveLaw::Parameters& rValues,
    KinematicVariables& rThisKinematicVariables
    )
{
    const auto& r_options = rValues.GetOptions();
    const bool stress_option = r_options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool constitutive_matrix_option = r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    BoundedMatrix<double, 3, 3> rotation_matrix;

    BuildRotationSystem(rotation_matrix, strain_size);

    // Undo the rotation in strain, stress and C
    if (strain_size == 6) {
        BoundedMatrix<double, 6, 6> voigt_rotation_matrix;
        ConstitutiveLawUtilities<6>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
        rValues.GetStrainVector() = prod(trans(voigt_rotation_matrix), rValues.GetStrainVector());
        if (stress_option)
            rValues.GetStressVector() = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
        if (constitutive_matrix_option) {
            BoundedMatrix<double, 6, 6> aux;
            noalias(aux) = prod(trans(voigt_rotation_matrix), rValues.GetConstitutiveMatrix());
            noalias(rValues.GetConstitutiveMatrix()) = prod(aux, voigt_rotation_matrix);
        }
    } else if (strain_size == 3) {
        BoundedMatrix<double, 3, 3> voigt_rotation_matrix;
        ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(rotation_matrix, voigt_rotation_matrix);
        rValues.GetStrainVector() = prod(trans(voigt_rotation_matrix), rValues.GetStrainVector());
        if (stress_option)
            rValues.GetStressVector() = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
        if (constitutive_matrix_option) {
            BoundedMatrix<double, 3, 3> aux;
            noalias(aux) = prod(trans(voigt_rotation_matrix), rValues.GetConstitutiveMatrix());
            noalias(rValues.GetConstitutiveMatrix()) = prod(aux, voigt_rotation_matrix);
        }
    }
    // Now undo the rotation in F if required
    if (!UseElementProvidedStrain()) {
        BoundedMatrix<double, 3, 3> inv_rotation_matrix;
        double aux_det;
        MathUtils<double>::InvertMatrix3(rotation_matrix, inv_rotation_matrix, aux_det);
        rThisKinematicVariables.F = prod(inv_rotation_matrix, rThisKinematicVariables.F);
        rThisKinematicVariables.F = prod(rThisKinematicVariables.F, rotation_matrix);
        rValues.SetDeformationGradientF(rThisKinematicVariables.F);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& BaseSolidElement::CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    DeltaDisplacement.resize(number_of_nodes , dimension, false);

    for ( IndexType i_node = 0; i_node < number_of_nodes; i_node++ ) {
        const array_1d<double, 3 >& current_displacement  = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        const array_1d<double, 3 >& previous_displacement = GetGeometry()[i_node].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( IndexType j_dim = 0; j_dim < dimension; ++j_dim )
            DeltaDisplacement(i_node, j_dim) = current_displacement[j_dim] - previous_displacement[j_dim];
    }

    return DeltaDisplacement;

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    const GeometryType& r_geom = GetGeometry();
    double detJ0;
    if (UseGeometryIntegrationMethod()) {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom,
            this->IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
        MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
        const Matrix& rDN_De =
            GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
        GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    } else {
        const auto& integration_points =  this->IntegrationPoints();
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[PointNumber],rJ0);
        MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
        Matrix DN_De;
        GetGeometry().ShapeFunctionsLocalGradients(DN_De, integration_points[PointNumber]);
        GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ0, rDN_DX);
    }
    return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

double BaseSolidElement::CalculateDerivativesOnCurrentConfiguration(
    Matrix& rJ,
    Matrix& rInvJ,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    double detJ;
    if (UseGeometryIntegrationMethod()) {
        rJ = GetGeometry().Jacobian( rJ, PointNumber, ThisIntegrationMethod );
        const Matrix& DN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
        MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
        GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    } else{
        const auto& integration_points =  this->IntegrationPoints();
        rJ = GetGeometry().Jacobian( rJ, integration_points[PointNumber] );
        Matrix DN_De;
        GetGeometry().ShapeFunctionsLocalGradients(DN_De, integration_points[PointNumber]);
        MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
        GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    }
    return detJ;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> BaseSolidElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, 6, 24>& rB,
    const ConstitutiveLaw::VoigtSizeMatrixType& rD,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY
    // option 1
    BoundedMatrix<double, 6, 24> aux;
    BoundedMatrix<double, 24, 6> aux2;
    noalias(aux2) = trans(rB);
    noalias(aux) = prod(rD, rB); // now IntegrationWeight done outside
    noalias( rLeftHandSideMatrix ) += prod(aux2, aux);

    // option 3
    // BoundedMatrix<double, 24, 24> aux;
    // aux.clear();
    // noalias( rLeftHandSideMatrix ) += aux;

    // option 2
    // noalias( rLeftHandSideMatrix ) += prod( trans( rB ), Matrix(prod(rD, rB)));



//     BoundedMatrix<double, 6, 6> aux_D;
//     noalias(aux_D) = IntegrationWeight * rD;


// const double crLeftHandSideMatrix0 = rB(0,0)*rD(0,0) + rB(1,0)*rD(0,1) + rB(2,0)*rD(0,2) + rB(3,0)*rD(0,3) + rB(4,0)*rD(0,4) + rB(5,0)*rD(0,5);
// const double crLeftHandSideMatrix1 = rB(0,0)*rD(0,1) + rB(1,0)*rD(1,1) + rB(2,0)*rD(1,2) + rB(3,0)*rD(1,3) + rB(4,0)*rD(1,4) + rB(5,0)*rD(1,5);
// const double crLeftHandSideMatrix2 = rB(0,0)*rD(0,2) + rB(1,0)*rD(1,2) + rB(2,0)*rD(2,2) + rB(3,0)*rD(2,3) + rB(4,0)*rD(2,4) + rB(5,0)*rD(2,5);
// const double crLeftHandSideMatrix3 = rB(0,0)*rD(0,3) + rB(1,0)*rD(1,3) + rB(2,0)*rD(2,3) + rB(3,0)*rD(3,3) + rB(4,0)*rD(3,4) + rB(5,0)*rD(3,5);
// const double crLeftHandSideMatrix4 = rB(0,0)*rD(0,4) + rB(1,0)*rD(1,4) + rB(2,0)*rD(2,4) + rB(3,0)*rD(3,4) + rB(4,0)*rD(4,4) + rB(5,0)*rD(4,5);
// const double crLeftHandSideMatrix5 = rB(0,0)*rD(0,5) + rB(1,0)*rD(1,5) + rB(2,0)*rD(2,5) + rB(3,0)*rD(3,5) + rB(4,0)*rD(4,5) + rB(5,0)*rD(5,5);
// const double crLeftHandSideMatrix6 = rB(0,1)*rD(0,0) + rB(1,1)*rD(0,1) + rB(2,1)*rD(0,2) + rB(3,1)*rD(0,3) + rB(4,1)*rD(0,4) + rB(5,1)*rD(0,5);
// const double crLeftHandSideMatrix7 = rB(0,1)*rD(0,1) + rB(1,1)*rD(1,1) + rB(2,1)*rD(1,2) + rB(3,1)*rD(1,3) + rB(4,1)*rD(1,4) + rB(5,1)*rD(1,5);
// const double crLeftHandSideMatrix8 = rB(0,1)*rD(0,2) + rB(1,1)*rD(1,2) + rB(2,1)*rD(2,2) + rB(3,1)*rD(2,3) + rB(4,1)*rD(2,4) + rB(5,1)*rD(2,5);
// const double crLeftHandSideMatrix9 = rB(0,1)*rD(0,3) + rB(1,1)*rD(1,3) + rB(2,1)*rD(2,3) + rB(3,1)*rD(3,3) + rB(4,1)*rD(3,4) + rB(5,1)*rD(3,5);
// const double crLeftHandSideMatrix10 = rB(0,1)*rD(0,4) + rB(1,1)*rD(1,4) + rB(2,1)*rD(2,4) + rB(3,1)*rD(3,4) + rB(4,1)*rD(4,4) + rB(5,1)*rD(4,5);
// const double crLeftHandSideMatrix11 = rB(0,1)*rD(0,5) + rB(1,1)*rD(1,5) + rB(2,1)*rD(2,5) + rB(3,1)*rD(3,5) + rB(4,1)*rD(4,5) + rB(5,1)*rD(5,5);
// const double crLeftHandSideMatrix12 = rB(0,2)*rD(0,0) + rB(1,2)*rD(0,1) + rB(2,2)*rD(0,2) + rB(3,2)*rD(0,3) + rB(4,2)*rD(0,4) + rB(5,2)*rD(0,5);
// const double crLeftHandSideMatrix13 = rB(0,2)*rD(0,1) + rB(1,2)*rD(1,1) + rB(2,2)*rD(1,2) + rB(3,2)*rD(1,3) + rB(4,2)*rD(1,4) + rB(5,2)*rD(1,5);
// const double crLeftHandSideMatrix14 = rB(0,2)*rD(0,2) + rB(1,2)*rD(1,2) + rB(2,2)*rD(2,2) + rB(3,2)*rD(2,3) + rB(4,2)*rD(2,4) + rB(5,2)*rD(2,5);
// const double crLeftHandSideMatrix15 = rB(0,2)*rD(0,3) + rB(1,2)*rD(1,3) + rB(2,2)*rD(2,3) + rB(3,2)*rD(3,3) + rB(4,2)*rD(3,4) + rB(5,2)*rD(3,5);
// const double crLeftHandSideMatrix16 = rB(0,2)*rD(0,4) + rB(1,2)*rD(1,4) + rB(2,2)*rD(2,4) + rB(3,2)*rD(3,4) + rB(4,2)*rD(4,4) + rB(5,2)*rD(4,5);
// const double crLeftHandSideMatrix17 = rB(0,2)*rD(0,5) + rB(1,2)*rD(1,5) + rB(2,2)*rD(2,5) + rB(3,2)*rD(3,5) + rB(4,2)*rD(4,5) + rB(5,2)*rD(5,5);
// const double crLeftHandSideMatrix18 = rB(0,3)*rD(0,0) + rB(1,3)*rD(0,1) + rB(2,3)*rD(0,2) + rB(3,3)*rD(0,3) + rB(4,3)*rD(0,4) + rB(5,3)*rD(0,5);
// const double crLeftHandSideMatrix19 = rB(0,3)*rD(0,1) + rB(1,3)*rD(1,1) + rB(2,3)*rD(1,2) + rB(3,3)*rD(1,3) + rB(4,3)*rD(1,4) + rB(5,3)*rD(1,5);
// const double crLeftHandSideMatrix20 = rB(0,3)*rD(0,2) + rB(1,3)*rD(1,2) + rB(2,3)*rD(2,2) + rB(3,3)*rD(2,3) + rB(4,3)*rD(2,4) + rB(5,3)*rD(2,5);
// const double crLeftHandSideMatrix21 = rB(0,3)*rD(0,3) + rB(1,3)*rD(1,3) + rB(2,3)*rD(2,3) + rB(3,3)*rD(3,3) + rB(4,3)*rD(3,4) + rB(5,3)*rD(3,5);
// const double crLeftHandSideMatrix22 = rB(0,3)*rD(0,4) + rB(1,3)*rD(1,4) + rB(2,3)*rD(2,4) + rB(3,3)*rD(3,4) + rB(4,3)*rD(4,4) + rB(5,3)*rD(4,5);
// const double crLeftHandSideMatrix23 = rB(0,3)*rD(0,5) + rB(1,3)*rD(1,5) + rB(2,3)*rD(2,5) + rB(3,3)*rD(3,5) + rB(4,3)*rD(4,5) + rB(5,3)*rD(5,5);
// const double crLeftHandSideMatrix24 = rB(0,4)*rD(0,0) + rB(1,4)*rD(0,1) + rB(2,4)*rD(0,2) + rB(3,4)*rD(0,3) + rB(4,4)*rD(0,4) + rB(5,4)*rD(0,5);
// const double crLeftHandSideMatrix25 = rB(0,4)*rD(0,1) + rB(1,4)*rD(1,1) + rB(2,4)*rD(1,2) + rB(3,4)*rD(1,3) + rB(4,4)*rD(1,4) + rB(5,4)*rD(1,5);
// const double crLeftHandSideMatrix26 = rB(0,4)*rD(0,2) + rB(1,4)*rD(1,2) + rB(2,4)*rD(2,2) + rB(3,4)*rD(2,3) + rB(4,4)*rD(2,4) + rB(5,4)*rD(2,5);
// const double crLeftHandSideMatrix27 = rB(0,4)*rD(0,3) + rB(1,4)*rD(1,3) + rB(2,4)*rD(2,3) + rB(3,4)*rD(3,3) + rB(4,4)*rD(3,4) + rB(5,4)*rD(3,5);
// const double crLeftHandSideMatrix28 = rB(0,4)*rD(0,4) + rB(1,4)*rD(1,4) + rB(2,4)*rD(2,4) + rB(3,4)*rD(3,4) + rB(4,4)*rD(4,4) + rB(5,4)*rD(4,5);
// const double crLeftHandSideMatrix29 = rB(0,4)*rD(0,5) + rB(1,4)*rD(1,5) + rB(2,4)*rD(2,5) + rB(3,4)*rD(3,5) + rB(4,4)*rD(4,5) + rB(5,4)*rD(5,5);
// const double crLeftHandSideMatrix30 = rB(0,5)*rD(0,0) + rB(1,5)*rD(0,1) + rB(2,5)*rD(0,2) + rB(3,5)*rD(0,3) + rB(4,5)*rD(0,4) + rB(5,5)*rD(0,5);
// const double crLeftHandSideMatrix31 = rB(0,5)*rD(0,1) + rB(1,5)*rD(1,1) + rB(2,5)*rD(1,2) + rB(3,5)*rD(1,3) + rB(4,5)*rD(1,4) + rB(5,5)*rD(1,5);
// const double crLeftHandSideMatrix32 = rB(0,5)*rD(0,2) + rB(1,5)*rD(1,2) + rB(2,5)*rD(2,2) + rB(3,5)*rD(2,3) + rB(4,5)*rD(2,4) + rB(5,5)*rD(2,5);
// const double crLeftHandSideMatrix33 = rB(0,5)*rD(0,3) + rB(1,5)*rD(1,3) + rB(2,5)*rD(2,3) + rB(3,5)*rD(3,3) + rB(4,5)*rD(3,4) + rB(5,5)*rD(3,5);
// const double crLeftHandSideMatrix34 = rB(0,5)*rD(0,4) + rB(1,5)*rD(1,4) + rB(2,5)*rD(2,4) + rB(3,5)*rD(3,4) + rB(4,5)*rD(4,4) + rB(5,5)*rD(4,5);
// const double crLeftHandSideMatrix35 = rB(0,5)*rD(0,5) + rB(1,5)*rD(1,5) + rB(2,5)*rD(2,5) + rB(3,5)*rD(3,5) + rB(4,5)*rD(4,5) + rB(5,5)*rD(5,5);
// const double crLeftHandSideMatrix36 = rB(0,6)*rD(0,0) + rB(1,6)*rD(0,1) + rB(2,6)*rD(0,2) + rB(3,6)*rD(0,3) + rB(4,6)*rD(0,4) + rB(5,6)*rD(0,5);
// const double crLeftHandSideMatrix37 = rB(0,6)*rD(0,1) + rB(1,6)*rD(1,1) + rB(2,6)*rD(1,2) + rB(3,6)*rD(1,3) + rB(4,6)*rD(1,4) + rB(5,6)*rD(1,5);
// const double crLeftHandSideMatrix38 = rB(0,6)*rD(0,2) + rB(1,6)*rD(1,2) + rB(2,6)*rD(2,2) + rB(3,6)*rD(2,3) + rB(4,6)*rD(2,4) + rB(5,6)*rD(2,5);
// const double crLeftHandSideMatrix39 = rB(0,6)*rD(0,3) + rB(1,6)*rD(1,3) + rB(2,6)*rD(2,3) + rB(3,6)*rD(3,3) + rB(4,6)*rD(3,4) + rB(5,6)*rD(3,5);
// const double crLeftHandSideMatrix40 = rB(0,6)*rD(0,4) + rB(1,6)*rD(1,4) + rB(2,6)*rD(2,4) + rB(3,6)*rD(3,4) + rB(4,6)*rD(4,4) + rB(5,6)*rD(4,5);
// const double crLeftHandSideMatrix41 = rB(0,6)*rD(0,5) + rB(1,6)*rD(1,5) + rB(2,6)*rD(2,5) + rB(3,6)*rD(3,5) + rB(4,6)*rD(4,5) + rB(5,6)*rD(5,5);
// const double crLeftHandSideMatrix42 = rB(0,7)*rD(0,0) + rB(1,7)*rD(0,1) + rB(2,7)*rD(0,2) + rB(3,7)*rD(0,3) + rB(4,7)*rD(0,4) + rB(5,7)*rD(0,5);
// const double crLeftHandSideMatrix43 = rB(0,7)*rD(0,1) + rB(1,7)*rD(1,1) + rB(2,7)*rD(1,2) + rB(3,7)*rD(1,3) + rB(4,7)*rD(1,4) + rB(5,7)*rD(1,5);
// const double crLeftHandSideMatrix44 = rB(0,7)*rD(0,2) + rB(1,7)*rD(1,2) + rB(2,7)*rD(2,2) + rB(3,7)*rD(2,3) + rB(4,7)*rD(2,4) + rB(5,7)*rD(2,5);
// const double crLeftHandSideMatrix45 = rB(0,7)*rD(0,3) + rB(1,7)*rD(1,3) + rB(2,7)*rD(2,3) + rB(3,7)*rD(3,3) + rB(4,7)*rD(3,4) + rB(5,7)*rD(3,5);
// const double crLeftHandSideMatrix46 = rB(0,7)*rD(0,4) + rB(1,7)*rD(1,4) + rB(2,7)*rD(2,4) + rB(3,7)*rD(3,4) + rB(4,7)*rD(4,4) + rB(5,7)*rD(4,5);
// const double crLeftHandSideMatrix47 = rB(0,7)*rD(0,5) + rB(1,7)*rD(1,5) + rB(2,7)*rD(2,5) + rB(3,7)*rD(3,5) + rB(4,7)*rD(4,5) + rB(5,7)*rD(5,5);
// const double crLeftHandSideMatrix48 = rB(0,8)*rD(0,0) + rB(1,8)*rD(0,1) + rB(2,8)*rD(0,2) + rB(3,8)*rD(0,3) + rB(4,8)*rD(0,4) + rB(5,8)*rD(0,5);
// const double crLeftHandSideMatrix49 = rB(0,8)*rD(0,1) + rB(1,8)*rD(1,1) + rB(2,8)*rD(1,2) + rB(3,8)*rD(1,3) + rB(4,8)*rD(1,4) + rB(5,8)*rD(1,5);
// const double crLeftHandSideMatrix50 = rB(0,8)*rD(0,2) + rB(1,8)*rD(1,2) + rB(2,8)*rD(2,2) + rB(3,8)*rD(2,3) + rB(4,8)*rD(2,4) + rB(5,8)*rD(2,5);
// const double crLeftHandSideMatrix51 = rB(0,8)*rD(0,3) + rB(1,8)*rD(1,3) + rB(2,8)*rD(2,3) + rB(3,8)*rD(3,3) + rB(4,8)*rD(3,4) + rB(5,8)*rD(3,5);
// const double crLeftHandSideMatrix52 = rB(0,8)*rD(0,4) + rB(1,8)*rD(1,4) + rB(2,8)*rD(2,4) + rB(3,8)*rD(3,4) + rB(4,8)*rD(4,4) + rB(5,8)*rD(4,5);
// const double crLeftHandSideMatrix53 = rB(0,8)*rD(0,5) + rB(1,8)*rD(1,5) + rB(2,8)*rD(2,5) + rB(3,8)*rD(3,5) + rB(4,8)*rD(4,5) + rB(5,8)*rD(5,5);
// const double crLeftHandSideMatrix54 = rB(0,9)*rD(0,0) + rB(1,9)*rD(0,1) + rB(2,9)*rD(0,2) + rB(3,9)*rD(0,3) + rB(4,9)*rD(0,4) + rB(5,9)*rD(0,5);
// const double crLeftHandSideMatrix55 = rB(0,9)*rD(0,1) + rB(1,9)*rD(1,1) + rB(2,9)*rD(1,2) + rB(3,9)*rD(1,3) + rB(4,9)*rD(1,4) + rB(5,9)*rD(1,5);
// const double crLeftHandSideMatrix56 = rB(0,9)*rD(0,2) + rB(1,9)*rD(1,2) + rB(2,9)*rD(2,2) + rB(3,9)*rD(2,3) + rB(4,9)*rD(2,4) + rB(5,9)*rD(2,5);
// const double crLeftHandSideMatrix57 = rB(0,9)*rD(0,3) + rB(1,9)*rD(1,3) + rB(2,9)*rD(2,3) + rB(3,9)*rD(3,3) + rB(4,9)*rD(3,4) + rB(5,9)*rD(3,5);
// const double crLeftHandSideMatrix58 = rB(0,9)*rD(0,4) + rB(1,9)*rD(1,4) + rB(2,9)*rD(2,4) + rB(3,9)*rD(3,4) + rB(4,9)*rD(4,4) + rB(5,9)*rD(4,5);
// const double crLeftHandSideMatrix59 = rB(0,9)*rD(0,5) + rB(1,9)*rD(1,5) + rB(2,9)*rD(2,5) + rB(3,9)*rD(3,5) + rB(4,9)*rD(4,5) + rB(5,9)*rD(5,5);
// const double crLeftHandSideMatrix60 = rB(0,10)*rD(0,0) + rB(1,10)*rD(0,1) + rB(2,10)*rD(0,2) + rB(3,10)*rD(0,3) + rB(4,10)*rD(0,4) + rB(5,10)*rD(0,5);
// const double crLeftHandSideMatrix61 = rB(0,10)*rD(0,1) + rB(1,10)*rD(1,1) + rB(2,10)*rD(1,2) + rB(3,10)*rD(1,3) + rB(4,10)*rD(1,4) + rB(5,10)*rD(1,5);
// const double crLeftHandSideMatrix62 = rB(0,10)*rD(0,2) + rB(1,10)*rD(1,2) + rB(2,10)*rD(2,2) + rB(3,10)*rD(2,3) + rB(4,10)*rD(2,4) + rB(5,10)*rD(2,5);
// const double crLeftHandSideMatrix63 = rB(0,10)*rD(0,3) + rB(1,10)*rD(1,3) + rB(2,10)*rD(2,3) + rB(3,10)*rD(3,3) + rB(4,10)*rD(3,4) + rB(5,10)*rD(3,5);
// const double crLeftHandSideMatrix64 = rB(0,10)*rD(0,4) + rB(1,10)*rD(1,4) + rB(2,10)*rD(2,4) + rB(3,10)*rD(3,4) + rB(4,10)*rD(4,4) + rB(5,10)*rD(4,5);
// const double crLeftHandSideMatrix65 = rB(0,10)*rD(0,5) + rB(1,10)*rD(1,5) + rB(2,10)*rD(2,5) + rB(3,10)*rD(3,5) + rB(4,10)*rD(4,5) + rB(5,10)*rD(5,5);
// const double crLeftHandSideMatrix66 = rB(0,11)*rD(0,0) + rB(1,11)*rD(0,1) + rB(2,11)*rD(0,2) + rB(3,11)*rD(0,3) + rB(4,11)*rD(0,4) + rB(5,11)*rD(0,5);
// const double crLeftHandSideMatrix67 = rB(0,11)*rD(0,1) + rB(1,11)*rD(1,1) + rB(2,11)*rD(1,2) + rB(3,11)*rD(1,3) + rB(4,11)*rD(1,4) + rB(5,11)*rD(1,5);
// const double crLeftHandSideMatrix68 = rB(0,11)*rD(0,2) + rB(1,11)*rD(1,2) + rB(2,11)*rD(2,2) + rB(3,11)*rD(2,3) + rB(4,11)*rD(2,4) + rB(5,11)*rD(2,5);
// const double crLeftHandSideMatrix69 = rB(0,11)*rD(0,3) + rB(1,11)*rD(1,3) + rB(2,11)*rD(2,3) + rB(3,11)*rD(3,3) + rB(4,11)*rD(3,4) + rB(5,11)*rD(3,5);
// const double crLeftHandSideMatrix70 = rB(0,11)*rD(0,4) + rB(1,11)*rD(1,4) + rB(2,11)*rD(2,4) + rB(3,11)*rD(3,4) + rB(4,11)*rD(4,4) + rB(5,11)*rD(4,5);
// const double crLeftHandSideMatrix71 = rB(0,11)*rD(0,5) + rB(1,11)*rD(1,5) + rB(2,11)*rD(2,5) + rB(3,11)*rD(3,5) + rB(4,11)*rD(4,5) + rB(5,11)*rD(5,5);
// const double crLeftHandSideMatrix72 = rB(0,12)*rD(0,0) + rB(1,12)*rD(0,1) + rB(2,12)*rD(0,2) + rB(3,12)*rD(0,3) + rB(4,12)*rD(0,4) + rB(5,12)*rD(0,5);
// const double crLeftHandSideMatrix73 = rB(0,12)*rD(0,1) + rB(1,12)*rD(1,1) + rB(2,12)*rD(1,2) + rB(3,12)*rD(1,3) + rB(4,12)*rD(1,4) + rB(5,12)*rD(1,5);
// const double crLeftHandSideMatrix74 = rB(0,12)*rD(0,2) + rB(1,12)*rD(1,2) + rB(2,12)*rD(2,2) + rB(3,12)*rD(2,3) + rB(4,12)*rD(2,4) + rB(5,12)*rD(2,5);
// const double crLeftHandSideMatrix75 = rB(0,12)*rD(0,3) + rB(1,12)*rD(1,3) + rB(2,12)*rD(2,3) + rB(3,12)*rD(3,3) + rB(4,12)*rD(3,4) + rB(5,12)*rD(3,5);
// const double crLeftHandSideMatrix76 = rB(0,12)*rD(0,4) + rB(1,12)*rD(1,4) + rB(2,12)*rD(2,4) + rB(3,12)*rD(3,4) + rB(4,12)*rD(4,4) + rB(5,12)*rD(4,5);
// const double crLeftHandSideMatrix77 = rB(0,12)*rD(0,5) + rB(1,12)*rD(1,5) + rB(2,12)*rD(2,5) + rB(3,12)*rD(3,5) + rB(4,12)*rD(4,5) + rB(5,12)*rD(5,5);
// const double crLeftHandSideMatrix78 = rB(0,13)*rD(0,0) + rB(1,13)*rD(0,1) + rB(2,13)*rD(0,2) + rB(3,13)*rD(0,3) + rB(4,13)*rD(0,4) + rB(5,13)*rD(0,5);
// const double crLeftHandSideMatrix79 = rB(0,13)*rD(0,1) + rB(1,13)*rD(1,1) + rB(2,13)*rD(1,2) + rB(3,13)*rD(1,3) + rB(4,13)*rD(1,4) + rB(5,13)*rD(1,5);
// const double crLeftHandSideMatrix80 = rB(0,13)*rD(0,2) + rB(1,13)*rD(1,2) + rB(2,13)*rD(2,2) + rB(3,13)*rD(2,3) + rB(4,13)*rD(2,4) + rB(5,13)*rD(2,5);
// const double crLeftHandSideMatrix81 = rB(0,13)*rD(0,3) + rB(1,13)*rD(1,3) + rB(2,13)*rD(2,3) + rB(3,13)*rD(3,3) + rB(4,13)*rD(3,4) + rB(5,13)*rD(3,5);
// const double crLeftHandSideMatrix82 = rB(0,13)*rD(0,4) + rB(1,13)*rD(1,4) + rB(2,13)*rD(2,4) + rB(3,13)*rD(3,4) + rB(4,13)*rD(4,4) + rB(5,13)*rD(4,5);
// const double crLeftHandSideMatrix83 = rB(0,13)*rD(0,5) + rB(1,13)*rD(1,5) + rB(2,13)*rD(2,5) + rB(3,13)*rD(3,5) + rB(4,13)*rD(4,5) + rB(5,13)*rD(5,5);
// const double crLeftHandSideMatrix84 = rB(0,14)*rD(0,0) + rB(1,14)*rD(0,1) + rB(2,14)*rD(0,2) + rB(3,14)*rD(0,3) + rB(4,14)*rD(0,4) + rB(5,14)*rD(0,5);
// const double crLeftHandSideMatrix85 = rB(0,14)*rD(0,1) + rB(1,14)*rD(1,1) + rB(2,14)*rD(1,2) + rB(3,14)*rD(1,3) + rB(4,14)*rD(1,4) + rB(5,14)*rD(1,5);
// const double crLeftHandSideMatrix86 = rB(0,14)*rD(0,2) + rB(1,14)*rD(1,2) + rB(2,14)*rD(2,2) + rB(3,14)*rD(2,3) + rB(4,14)*rD(2,4) + rB(5,14)*rD(2,5);
// const double crLeftHandSideMatrix87 = rB(0,14)*rD(0,3) + rB(1,14)*rD(1,3) + rB(2,14)*rD(2,3) + rB(3,14)*rD(3,3) + rB(4,14)*rD(3,4) + rB(5,14)*rD(3,5);
// const double crLeftHandSideMatrix88 = rB(0,14)*rD(0,4) + rB(1,14)*rD(1,4) + rB(2,14)*rD(2,4) + rB(3,14)*rD(3,4) + rB(4,14)*rD(4,4) + rB(5,14)*rD(4,5);
// const double crLeftHandSideMatrix89 = rB(0,14)*rD(0,5) + rB(1,14)*rD(1,5) + rB(2,14)*rD(2,5) + rB(3,14)*rD(3,5) + rB(4,14)*rD(4,5) + rB(5,14)*rD(5,5);
// const double crLeftHandSideMatrix90 = rB(0,15)*rD(0,0) + rB(1,15)*rD(0,1) + rB(2,15)*rD(0,2) + rB(3,15)*rD(0,3) + rB(4,15)*rD(0,4) + rB(5,15)*rD(0,5);
// const double crLeftHandSideMatrix91 = rB(0,15)*rD(0,1) + rB(1,15)*rD(1,1) + rB(2,15)*rD(1,2) + rB(3,15)*rD(1,3) + rB(4,15)*rD(1,4) + rB(5,15)*rD(1,5);
// const double crLeftHandSideMatrix92 = rB(0,15)*rD(0,2) + rB(1,15)*rD(1,2) + rB(2,15)*rD(2,2) + rB(3,15)*rD(2,3) + rB(4,15)*rD(2,4) + rB(5,15)*rD(2,5);
// const double crLeftHandSideMatrix93 = rB(0,15)*rD(0,3) + rB(1,15)*rD(1,3) + rB(2,15)*rD(2,3) + rB(3,15)*rD(3,3) + rB(4,15)*rD(3,4) + rB(5,15)*rD(3,5);
// const double crLeftHandSideMatrix94 = rB(0,15)*rD(0,4) + rB(1,15)*rD(1,4) + rB(2,15)*rD(2,4) + rB(3,15)*rD(3,4) + rB(4,15)*rD(4,4) + rB(5,15)*rD(4,5);
// const double crLeftHandSideMatrix95 = rB(0,15)*rD(0,5) + rB(1,15)*rD(1,5) + rB(2,15)*rD(2,5) + rB(3,15)*rD(3,5) + rB(4,15)*rD(4,5) + rB(5,15)*rD(5,5);
// const double crLeftHandSideMatrix96 = rB(0,16)*rD(0,0) + rB(1,16)*rD(0,1) + rB(2,16)*rD(0,2) + rB(3,16)*rD(0,3) + rB(4,16)*rD(0,4) + rB(5,16)*rD(0,5);
// const double crLeftHandSideMatrix97 = rB(0,16)*rD(0,1) + rB(1,16)*rD(1,1) + rB(2,16)*rD(1,2) + rB(3,16)*rD(1,3) + rB(4,16)*rD(1,4) + rB(5,16)*rD(1,5);
// const double crLeftHandSideMatrix98 = rB(0,16)*rD(0,2) + rB(1,16)*rD(1,2) + rB(2,16)*rD(2,2) + rB(3,16)*rD(2,3) + rB(4,16)*rD(2,4) + rB(5,16)*rD(2,5);
// const double crLeftHandSideMatrix99 = rB(0,16)*rD(0,3) + rB(1,16)*rD(1,3) + rB(2,16)*rD(2,3) + rB(3,16)*rD(3,3) + rB(4,16)*rD(3,4) + rB(5,16)*rD(3,5);
// const double crLeftHandSideMatrix100 = rB(0,16)*rD(0,4) + rB(1,16)*rD(1,4) + rB(2,16)*rD(2,4) + rB(3,16)*rD(3,4) + rB(4,16)*rD(4,4) + rB(5,16)*rD(4,5);
// const double crLeftHandSideMatrix101 = rB(0,16)*rD(0,5) + rB(1,16)*rD(1,5) + rB(2,16)*rD(2,5) + rB(3,16)*rD(3,5) + rB(4,16)*rD(4,5) + rB(5,16)*rD(5,5);
// const double crLeftHandSideMatrix102 = rB(0,17)*rD(0,0) + rB(1,17)*rD(0,1) + rB(2,17)*rD(0,2) + rB(3,17)*rD(0,3) + rB(4,17)*rD(0,4) + rB(5,17)*rD(0,5);
// const double crLeftHandSideMatrix103 = rB(0,17)*rD(0,1) + rB(1,17)*rD(1,1) + rB(2,17)*rD(1,2) + rB(3,17)*rD(1,3) + rB(4,17)*rD(1,4) + rB(5,17)*rD(1,5);
// const double crLeftHandSideMatrix104 = rB(0,17)*rD(0,2) + rB(1,17)*rD(1,2) + rB(2,17)*rD(2,2) + rB(3,17)*rD(2,3) + rB(4,17)*rD(2,4) + rB(5,17)*rD(2,5);
// const double crLeftHandSideMatrix105 = rB(0,17)*rD(0,3) + rB(1,17)*rD(1,3) + rB(2,17)*rD(2,3) + rB(3,17)*rD(3,3) + rB(4,17)*rD(3,4) + rB(5,17)*rD(3,5);
// const double crLeftHandSideMatrix106 = rB(0,17)*rD(0,4) + rB(1,17)*rD(1,4) + rB(2,17)*rD(2,4) + rB(3,17)*rD(3,4) + rB(4,17)*rD(4,4) + rB(5,17)*rD(4,5);
// const double crLeftHandSideMatrix107 = rB(0,17)*rD(0,5) + rB(1,17)*rD(1,5) + rB(2,17)*rD(2,5) + rB(3,17)*rD(3,5) + rB(4,17)*rD(4,5) + rB(5,17)*rD(5,5);
// const double crLeftHandSideMatrix108 = rB(0,18)*rD(0,0) + rB(1,18)*rD(0,1) + rB(2,18)*rD(0,2) + rB(3,18)*rD(0,3) + rB(4,18)*rD(0,4) + rB(5,18)*rD(0,5);
// const double crLeftHandSideMatrix109 = rB(0,18)*rD(0,1) + rB(1,18)*rD(1,1) + rB(2,18)*rD(1,2) + rB(3,18)*rD(1,3) + rB(4,18)*rD(1,4) + rB(5,18)*rD(1,5);
// const double crLeftHandSideMatrix110 = rB(0,18)*rD(0,2) + rB(1,18)*rD(1,2) + rB(2,18)*rD(2,2) + rB(3,18)*rD(2,3) + rB(4,18)*rD(2,4) + rB(5,18)*rD(2,5);
// const double crLeftHandSideMatrix111 = rB(0,18)*rD(0,3) + rB(1,18)*rD(1,3) + rB(2,18)*rD(2,3) + rB(3,18)*rD(3,3) + rB(4,18)*rD(3,4) + rB(5,18)*rD(3,5);
// const double crLeftHandSideMatrix112 = rB(0,18)*rD(0,4) + rB(1,18)*rD(1,4) + rB(2,18)*rD(2,4) + rB(3,18)*rD(3,4) + rB(4,18)*rD(4,4) + rB(5,18)*rD(4,5);
// const double crLeftHandSideMatrix113 = rB(0,18)*rD(0,5) + rB(1,18)*rD(1,5) + rB(2,18)*rD(2,5) + rB(3,18)*rD(3,5) + rB(4,18)*rD(4,5) + rB(5,18)*rD(5,5);
// const double crLeftHandSideMatrix114 = rB(0,19)*rD(0,0) + rB(1,19)*rD(0,1) + rB(2,19)*rD(0,2) + rB(3,19)*rD(0,3) + rB(4,19)*rD(0,4) + rB(5,19)*rD(0,5);
// const double crLeftHandSideMatrix115 = rB(0,19)*rD(0,1) + rB(1,19)*rD(1,1) + rB(2,19)*rD(1,2) + rB(3,19)*rD(1,3) + rB(4,19)*rD(1,4) + rB(5,19)*rD(1,5);
// const double crLeftHandSideMatrix116 = rB(0,19)*rD(0,2) + rB(1,19)*rD(1,2) + rB(2,19)*rD(2,2) + rB(3,19)*rD(2,3) + rB(4,19)*rD(2,4) + rB(5,19)*rD(2,5);
// const double crLeftHandSideMatrix117 = rB(0,19)*rD(0,3) + rB(1,19)*rD(1,3) + rB(2,19)*rD(2,3) + rB(3,19)*rD(3,3) + rB(4,19)*rD(3,4) + rB(5,19)*rD(3,5);
// const double crLeftHandSideMatrix118 = rB(0,19)*rD(0,4) + rB(1,19)*rD(1,4) + rB(2,19)*rD(2,4) + rB(3,19)*rD(3,4) + rB(4,19)*rD(4,4) + rB(5,19)*rD(4,5);
// const double crLeftHandSideMatrix119 = rB(0,19)*rD(0,5) + rB(1,19)*rD(1,5) + rB(2,19)*rD(2,5) + rB(3,19)*rD(3,5) + rB(4,19)*rD(4,5) + rB(5,19)*rD(5,5);
// const double crLeftHandSideMatrix120 = rB(0,20)*rD(0,0) + rB(1,20)*rD(0,1) + rB(2,20)*rD(0,2) + rB(3,20)*rD(0,3) + rB(4,20)*rD(0,4) + rB(5,20)*rD(0,5);
// const double crLeftHandSideMatrix121 = rB(0,20)*rD(0,1) + rB(1,20)*rD(1,1) + rB(2,20)*rD(1,2) + rB(3,20)*rD(1,3) + rB(4,20)*rD(1,4) + rB(5,20)*rD(1,5);
// const double crLeftHandSideMatrix122 = rB(0,20)*rD(0,2) + rB(1,20)*rD(1,2) + rB(2,20)*rD(2,2) + rB(3,20)*rD(2,3) + rB(4,20)*rD(2,4) + rB(5,20)*rD(2,5);
// const double crLeftHandSideMatrix123 = rB(0,20)*rD(0,3) + rB(1,20)*rD(1,3) + rB(2,20)*rD(2,3) + rB(3,20)*rD(3,3) + rB(4,20)*rD(3,4) + rB(5,20)*rD(3,5);
// const double crLeftHandSideMatrix124 = rB(0,20)*rD(0,4) + rB(1,20)*rD(1,4) + rB(2,20)*rD(2,4) + rB(3,20)*rD(3,4) + rB(4,20)*rD(4,4) + rB(5,20)*rD(4,5);
// const double crLeftHandSideMatrix125 = rB(0,20)*rD(0,5) + rB(1,20)*rD(1,5) + rB(2,20)*rD(2,5) + rB(3,20)*rD(3,5) + rB(4,20)*rD(4,5) + rB(5,20)*rD(5,5);
// const double crLeftHandSideMatrix126 = rB(0,21)*rD(0,0) + rB(1,21)*rD(0,1) + rB(2,21)*rD(0,2) + rB(3,21)*rD(0,3) + rB(4,21)*rD(0,4) + rB(5,21)*rD(0,5);
// const double crLeftHandSideMatrix127 = rB(0,21)*rD(0,1) + rB(1,21)*rD(1,1) + rB(2,21)*rD(1,2) + rB(3,21)*rD(1,3) + rB(4,21)*rD(1,4) + rB(5,21)*rD(1,5);
// const double crLeftHandSideMatrix128 = rB(0,21)*rD(0,2) + rB(1,21)*rD(1,2) + rB(2,21)*rD(2,2) + rB(3,21)*rD(2,3) + rB(4,21)*rD(2,4) + rB(5,21)*rD(2,5);
// const double crLeftHandSideMatrix129 = rB(0,21)*rD(0,3) + rB(1,21)*rD(1,3) + rB(2,21)*rD(2,3) + rB(3,21)*rD(3,3) + rB(4,21)*rD(3,4) + rB(5,21)*rD(3,5);
// const double crLeftHandSideMatrix130 = rB(0,21)*rD(0,4) + rB(1,21)*rD(1,4) + rB(2,21)*rD(2,4) + rB(3,21)*rD(3,4) + rB(4,21)*rD(4,4) + rB(5,21)*rD(4,5);
// const double crLeftHandSideMatrix131 = rB(0,21)*rD(0,5) + rB(1,21)*rD(1,5) + rB(2,21)*rD(2,5) + rB(3,21)*rD(3,5) + rB(4,21)*rD(4,5) + rB(5,21)*rD(5,5);
// const double crLeftHandSideMatrix132 = rB(0,22)*rD(0,0) + rB(1,22)*rD(0,1) + rB(2,22)*rD(0,2) + rB(3,22)*rD(0,3) + rB(4,22)*rD(0,4) + rB(5,22)*rD(0,5);
// const double crLeftHandSideMatrix133 = rB(0,22)*rD(0,1) + rB(1,22)*rD(1,1) + rB(2,22)*rD(1,2) + rB(3,22)*rD(1,3) + rB(4,22)*rD(1,4) + rB(5,22)*rD(1,5);
// const double crLeftHandSideMatrix134 = rB(0,22)*rD(0,2) + rB(1,22)*rD(1,2) + rB(2,22)*rD(2,2) + rB(3,22)*rD(2,3) + rB(4,22)*rD(2,4) + rB(5,22)*rD(2,5);
// const double crLeftHandSideMatrix135 = rB(0,22)*rD(0,3) + rB(1,22)*rD(1,3) + rB(2,22)*rD(2,3) + rB(3,22)*rD(3,3) + rB(4,22)*rD(3,4) + rB(5,22)*rD(3,5);
// const double crLeftHandSideMatrix136 = rB(0,22)*rD(0,4) + rB(1,22)*rD(1,4) + rB(2,22)*rD(2,4) + rB(3,22)*rD(3,4) + rB(4,22)*rD(4,4) + rB(5,22)*rD(4,5);
// const double crLeftHandSideMatrix137 = rB(0,22)*rD(0,5) + rB(1,22)*rD(1,5) + rB(2,22)*rD(2,5) + rB(3,22)*rD(3,5) + rB(4,22)*rD(4,5) + rB(5,22)*rD(5,5);
// const double crLeftHandSideMatrix138 = rB(0,23)*rD(0,0) + rB(1,23)*rD(0,1) + rB(2,23)*rD(0,2) + rB(3,23)*rD(0,3) + rB(4,23)*rD(0,4) + rB(5,23)*rD(0,5);
// const double crLeftHandSideMatrix139 = rB(0,23)*rD(0,1) + rB(1,23)*rD(1,1) + rB(2,23)*rD(1,2) + rB(3,23)*rD(1,3) + rB(4,23)*rD(1,4) + rB(5,23)*rD(1,5);
// const double crLeftHandSideMatrix140 = rB(0,23)*rD(0,2) + rB(1,23)*rD(1,2) + rB(2,23)*rD(2,2) + rB(3,23)*rD(2,3) + rB(4,23)*rD(2,4) + rB(5,23)*rD(2,5);
// const double crLeftHandSideMatrix141 = rB(0,23)*rD(0,3) + rB(1,23)*rD(1,3) + rB(2,23)*rD(2,3) + rB(3,23)*rD(3,3) + rB(4,23)*rD(3,4) + rB(5,23)*rD(3,5);
// const double crLeftHandSideMatrix142 = rB(0,23)*rD(0,4) + rB(1,23)*rD(1,4) + rB(2,23)*rD(2,4) + rB(3,23)*rD(3,4) + rB(4,23)*rD(4,4) + rB(5,23)*rD(4,5);
// const double crLeftHandSideMatrix143 = rB(0,23)*rD(0,5) + rB(1,23)*rD(1,5) + rB(2,23)*rD(2,5) + rB(3,23)*rD(3,5) + rB(4,23)*rD(4,5) + rB(5,23)*rD(5,5);
// rLeftHandSideMatrix(0,0)+=(crLeftHandSideMatrix0*rB(0,0) + crLeftHandSideMatrix1*rB(1,0) + crLeftHandSideMatrix2*rB(2,0) + crLeftHandSideMatrix3*rB(3,0) + crLeftHandSideMatrix4*rB(4,0) + crLeftHandSideMatrix5*rB(5,0));
// rLeftHandSideMatrix(0,1)+=(crLeftHandSideMatrix0*rB(0,1) + crLeftHandSideMatrix1*rB(1,1) + crLeftHandSideMatrix2*rB(2,1) + crLeftHandSideMatrix3*rB(3,1) + crLeftHandSideMatrix4*rB(4,1) + crLeftHandSideMatrix5*rB(5,1));
// rLeftHandSideMatrix(0,2)+=(crLeftHandSideMatrix0*rB(0,2) + crLeftHandSideMatrix1*rB(1,2) + crLeftHandSideMatrix2*rB(2,2) + crLeftHandSideMatrix3*rB(3,2) + crLeftHandSideMatrix4*rB(4,2) + crLeftHandSideMatrix5*rB(5,2));
// rLeftHandSideMatrix(0,3)+=(crLeftHandSideMatrix0*rB(0,3) + crLeftHandSideMatrix1*rB(1,3) + crLeftHandSideMatrix2*rB(2,3) + crLeftHandSideMatrix3*rB(3,3) + crLeftHandSideMatrix4*rB(4,3) + crLeftHandSideMatrix5*rB(5,3));
// rLeftHandSideMatrix(0,4)+=(crLeftHandSideMatrix0*rB(0,4) + crLeftHandSideMatrix1*rB(1,4) + crLeftHandSideMatrix2*rB(2,4) + crLeftHandSideMatrix3*rB(3,4) + crLeftHandSideMatrix4*rB(4,4) + crLeftHandSideMatrix5*rB(5,4));
// rLeftHandSideMatrix(0,5)+=(crLeftHandSideMatrix0*rB(0,5) + crLeftHandSideMatrix1*rB(1,5) + crLeftHandSideMatrix2*rB(2,5) + crLeftHandSideMatrix3*rB(3,5) + crLeftHandSideMatrix4*rB(4,5) + crLeftHandSideMatrix5*rB(5,5));
// rLeftHandSideMatrix(0,6)+=(crLeftHandSideMatrix0*rB(0,6) + crLeftHandSideMatrix1*rB(1,6) + crLeftHandSideMatrix2*rB(2,6) + crLeftHandSideMatrix3*rB(3,6) + crLeftHandSideMatrix4*rB(4,6) + crLeftHandSideMatrix5*rB(5,6));
// rLeftHandSideMatrix(0,7)+=(crLeftHandSideMatrix0*rB(0,7) + crLeftHandSideMatrix1*rB(1,7) + crLeftHandSideMatrix2*rB(2,7) + crLeftHandSideMatrix3*rB(3,7) + crLeftHandSideMatrix4*rB(4,7) + crLeftHandSideMatrix5*rB(5,7));
// rLeftHandSideMatrix(0,8)+=(crLeftHandSideMatrix0*rB(0,8) + crLeftHandSideMatrix1*rB(1,8) + crLeftHandSideMatrix2*rB(2,8) + crLeftHandSideMatrix3*rB(3,8) + crLeftHandSideMatrix4*rB(4,8) + crLeftHandSideMatrix5*rB(5,8));
// rLeftHandSideMatrix(0,9)+=(crLeftHandSideMatrix0*rB(0,9) + crLeftHandSideMatrix1*rB(1,9) + crLeftHandSideMatrix2*rB(2,9) + crLeftHandSideMatrix3*rB(3,9) + crLeftHandSideMatrix4*rB(4,9) + crLeftHandSideMatrix5*rB(5,9));
// rLeftHandSideMatrix(0,10)+=(crLeftHandSideMatrix0*rB(0,10) + crLeftHandSideMatrix1*rB(1,10) + crLeftHandSideMatrix2*rB(2,10) + crLeftHandSideMatrix3*rB(3,10) + crLeftHandSideMatrix4*rB(4,10) + crLeftHandSideMatrix5*rB(5,10));
// rLeftHandSideMatrix(0,11)+=(crLeftHandSideMatrix0*rB(0,11) + crLeftHandSideMatrix1*rB(1,11) + crLeftHandSideMatrix2*rB(2,11) + crLeftHandSideMatrix3*rB(3,11) + crLeftHandSideMatrix4*rB(4,11) + crLeftHandSideMatrix5*rB(5,11));
// rLeftHandSideMatrix(0,12)+=(crLeftHandSideMatrix0*rB(0,12) + crLeftHandSideMatrix1*rB(1,12) + crLeftHandSideMatrix2*rB(2,12) + crLeftHandSideMatrix3*rB(3,12) + crLeftHandSideMatrix4*rB(4,12) + crLeftHandSideMatrix5*rB(5,12));
// rLeftHandSideMatrix(0,13)+=(crLeftHandSideMatrix0*rB(0,13) + crLeftHandSideMatrix1*rB(1,13) + crLeftHandSideMatrix2*rB(2,13) + crLeftHandSideMatrix3*rB(3,13) + crLeftHandSideMatrix4*rB(4,13) + crLeftHandSideMatrix5*rB(5,13));
// rLeftHandSideMatrix(0,14)+=(crLeftHandSideMatrix0*rB(0,14) + crLeftHandSideMatrix1*rB(1,14) + crLeftHandSideMatrix2*rB(2,14) + crLeftHandSideMatrix3*rB(3,14) + crLeftHandSideMatrix4*rB(4,14) + crLeftHandSideMatrix5*rB(5,14));
// rLeftHandSideMatrix(0,15)+=(crLeftHandSideMatrix0*rB(0,15) + crLeftHandSideMatrix1*rB(1,15) + crLeftHandSideMatrix2*rB(2,15) + crLeftHandSideMatrix3*rB(3,15) + crLeftHandSideMatrix4*rB(4,15) + crLeftHandSideMatrix5*rB(5,15));
// rLeftHandSideMatrix(0,16)+=(crLeftHandSideMatrix0*rB(0,16) + crLeftHandSideMatrix1*rB(1,16) + crLeftHandSideMatrix2*rB(2,16) + crLeftHandSideMatrix3*rB(3,16) + crLeftHandSideMatrix4*rB(4,16) + crLeftHandSideMatrix5*rB(5,16));
// rLeftHandSideMatrix(0,17)+=(crLeftHandSideMatrix0*rB(0,17) + crLeftHandSideMatrix1*rB(1,17) + crLeftHandSideMatrix2*rB(2,17) + crLeftHandSideMatrix3*rB(3,17) + crLeftHandSideMatrix4*rB(4,17) + crLeftHandSideMatrix5*rB(5,17));
// rLeftHandSideMatrix(0,18)+=(crLeftHandSideMatrix0*rB(0,18) + crLeftHandSideMatrix1*rB(1,18) + crLeftHandSideMatrix2*rB(2,18) + crLeftHandSideMatrix3*rB(3,18) + crLeftHandSideMatrix4*rB(4,18) + crLeftHandSideMatrix5*rB(5,18));
// rLeftHandSideMatrix(0,19)+=(crLeftHandSideMatrix0*rB(0,19) + crLeftHandSideMatrix1*rB(1,19) + crLeftHandSideMatrix2*rB(2,19) + crLeftHandSideMatrix3*rB(3,19) + crLeftHandSideMatrix4*rB(4,19) + crLeftHandSideMatrix5*rB(5,19));
// rLeftHandSideMatrix(0,20)+=(crLeftHandSideMatrix0*rB(0,20) + crLeftHandSideMatrix1*rB(1,20) + crLeftHandSideMatrix2*rB(2,20) + crLeftHandSideMatrix3*rB(3,20) + crLeftHandSideMatrix4*rB(4,20) + crLeftHandSideMatrix5*rB(5,20));
// rLeftHandSideMatrix(0,21)+=(crLeftHandSideMatrix0*rB(0,21) + crLeftHandSideMatrix1*rB(1,21) + crLeftHandSideMatrix2*rB(2,21) + crLeftHandSideMatrix3*rB(3,21) + crLeftHandSideMatrix4*rB(4,21) + crLeftHandSideMatrix5*rB(5,21));
// rLeftHandSideMatrix(0,22)+=(crLeftHandSideMatrix0*rB(0,22) + crLeftHandSideMatrix1*rB(1,22) + crLeftHandSideMatrix2*rB(2,22) + crLeftHandSideMatrix3*rB(3,22) + crLeftHandSideMatrix4*rB(4,22) + crLeftHandSideMatrix5*rB(5,22));
// rLeftHandSideMatrix(0,23)+=(crLeftHandSideMatrix0*rB(0,23) + crLeftHandSideMatrix1*rB(1,23) + crLeftHandSideMatrix2*rB(2,23) + crLeftHandSideMatrix3*rB(3,23) + crLeftHandSideMatrix4*rB(4,23) + crLeftHandSideMatrix5*rB(5,23));

// // rLeftHandSideMatrix(1,0)+=(crLeftHandSideMatrix10*rB(4,0) + crLeftHandSideMatrix11*rB(5,0) + crLeftHandSideMatrix6*rB(0,0) + crLeftHandSideMatrix7*rB(1,0) + crLeftHandSideMatrix8*rB(2,0) + crLeftHandSideMatrix9*rB(3,0));
// rLeftHandSideMatrix(1,1)+=(crLeftHandSideMatrix10*rB(4,1) + crLeftHandSideMatrix11*rB(5,1) + crLeftHandSideMatrix6*rB(0,1) + crLeftHandSideMatrix7*rB(1,1) + crLeftHandSideMatrix8*rB(2,1) + crLeftHandSideMatrix9*rB(3,1));
// rLeftHandSideMatrix(1,2)+=(crLeftHandSideMatrix10*rB(4,2) + crLeftHandSideMatrix11*rB(5,2) + crLeftHandSideMatrix6*rB(0,2) + crLeftHandSideMatrix7*rB(1,2) + crLeftHandSideMatrix8*rB(2,2) + crLeftHandSideMatrix9*rB(3,2));
// rLeftHandSideMatrix(1,3)+=(crLeftHandSideMatrix10*rB(4,3) + crLeftHandSideMatrix11*rB(5,3) + crLeftHandSideMatrix6*rB(0,3) + crLeftHandSideMatrix7*rB(1,3) + crLeftHandSideMatrix8*rB(2,3) + crLeftHandSideMatrix9*rB(3,3));
// rLeftHandSideMatrix(1,4)+=(crLeftHandSideMatrix10*rB(4,4) + crLeftHandSideMatrix11*rB(5,4) + crLeftHandSideMatrix6*rB(0,4) + crLeftHandSideMatrix7*rB(1,4) + crLeftHandSideMatrix8*rB(2,4) + crLeftHandSideMatrix9*rB(3,4));
// rLeftHandSideMatrix(1,5)+=(crLeftHandSideMatrix10*rB(4,5) + crLeftHandSideMatrix11*rB(5,5) + crLeftHandSideMatrix6*rB(0,5) + crLeftHandSideMatrix7*rB(1,5) + crLeftHandSideMatrix8*rB(2,5) + crLeftHandSideMatrix9*rB(3,5));
// rLeftHandSideMatrix(1,6)+=(crLeftHandSideMatrix10*rB(4,6) + crLeftHandSideMatrix11*rB(5,6) + crLeftHandSideMatrix6*rB(0,6) + crLeftHandSideMatrix7*rB(1,6) + crLeftHandSideMatrix8*rB(2,6) + crLeftHandSideMatrix9*rB(3,6));
// rLeftHandSideMatrix(1,7)+=(crLeftHandSideMatrix10*rB(4,7) + crLeftHandSideMatrix11*rB(5,7) + crLeftHandSideMatrix6*rB(0,7) + crLeftHandSideMatrix7*rB(1,7) + crLeftHandSideMatrix8*rB(2,7) + crLeftHandSideMatrix9*rB(3,7));
// rLeftHandSideMatrix(1,8)+=(crLeftHandSideMatrix10*rB(4,8) + crLeftHandSideMatrix11*rB(5,8) + crLeftHandSideMatrix6*rB(0,8) + crLeftHandSideMatrix7*rB(1,8) + crLeftHandSideMatrix8*rB(2,8) + crLeftHandSideMatrix9*rB(3,8));
// rLeftHandSideMatrix(1,9)+=(crLeftHandSideMatrix10*rB(4,9) + crLeftHandSideMatrix11*rB(5,9) + crLeftHandSideMatrix6*rB(0,9) + crLeftHandSideMatrix7*rB(1,9) + crLeftHandSideMatrix8*rB(2,9) + crLeftHandSideMatrix9*rB(3,9));
// rLeftHandSideMatrix(1,10)+=(crLeftHandSideMatrix10*rB(4,10) + crLeftHandSideMatrix11*rB(5,10) + crLeftHandSideMatrix6*rB(0,10) + crLeftHandSideMatrix7*rB(1,10) + crLeftHandSideMatrix8*rB(2,10) + crLeftHandSideMatrix9*rB(3,10));
// rLeftHandSideMatrix(1,11)+=(crLeftHandSideMatrix10*rB(4,11) + crLeftHandSideMatrix11*rB(5,11) + crLeftHandSideMatrix6*rB(0,11) + crLeftHandSideMatrix7*rB(1,11) + crLeftHandSideMatrix8*rB(2,11) + crLeftHandSideMatrix9*rB(3,11));
// rLeftHandSideMatrix(1,12)+=(crLeftHandSideMatrix10*rB(4,12) + crLeftHandSideMatrix11*rB(5,12) + crLeftHandSideMatrix6*rB(0,12) + crLeftHandSideMatrix7*rB(1,12) + crLeftHandSideMatrix8*rB(2,12) + crLeftHandSideMatrix9*rB(3,12));
// rLeftHandSideMatrix(1,13)+=(crLeftHandSideMatrix10*rB(4,13) + crLeftHandSideMatrix11*rB(5,13) + crLeftHandSideMatrix6*rB(0,13) + crLeftHandSideMatrix7*rB(1,13) + crLeftHandSideMatrix8*rB(2,13) + crLeftHandSideMatrix9*rB(3,13));
// rLeftHandSideMatrix(1,14)+=(crLeftHandSideMatrix10*rB(4,14) + crLeftHandSideMatrix11*rB(5,14) + crLeftHandSideMatrix6*rB(0,14) + crLeftHandSideMatrix7*rB(1,14) + crLeftHandSideMatrix8*rB(2,14) + crLeftHandSideMatrix9*rB(3,14));
// rLeftHandSideMatrix(1,15)+=(crLeftHandSideMatrix10*rB(4,15) + crLeftHandSideMatrix11*rB(5,15) + crLeftHandSideMatrix6*rB(0,15) + crLeftHandSideMatrix7*rB(1,15) + crLeftHandSideMatrix8*rB(2,15) + crLeftHandSideMatrix9*rB(3,15));
// rLeftHandSideMatrix(1,16)+=(crLeftHandSideMatrix10*rB(4,16) + crLeftHandSideMatrix11*rB(5,16) + crLeftHandSideMatrix6*rB(0,16) + crLeftHandSideMatrix7*rB(1,16) + crLeftHandSideMatrix8*rB(2,16) + crLeftHandSideMatrix9*rB(3,16));
// rLeftHandSideMatrix(1,17)+=(crLeftHandSideMatrix10*rB(4,17) + crLeftHandSideMatrix11*rB(5,17) + crLeftHandSideMatrix6*rB(0,17) + crLeftHandSideMatrix7*rB(1,17) + crLeftHandSideMatrix8*rB(2,17) + crLeftHandSideMatrix9*rB(3,17));
// rLeftHandSideMatrix(1,18)+=(crLeftHandSideMatrix10*rB(4,18) + crLeftHandSideMatrix11*rB(5,18) + crLeftHandSideMatrix6*rB(0,18) + crLeftHandSideMatrix7*rB(1,18) + crLeftHandSideMatrix8*rB(2,18) + crLeftHandSideMatrix9*rB(3,18));
// rLeftHandSideMatrix(1,19)+=(crLeftHandSideMatrix10*rB(4,19) + crLeftHandSideMatrix11*rB(5,19) + crLeftHandSideMatrix6*rB(0,19) + crLeftHandSideMatrix7*rB(1,19) + crLeftHandSideMatrix8*rB(2,19) + crLeftHandSideMatrix9*rB(3,19));
// rLeftHandSideMatrix(1,20)+=(crLeftHandSideMatrix10*rB(4,20) + crLeftHandSideMatrix11*rB(5,20) + crLeftHandSideMatrix6*rB(0,20) + crLeftHandSideMatrix7*rB(1,20) + crLeftHandSideMatrix8*rB(2,20) + crLeftHandSideMatrix9*rB(3,20));
// rLeftHandSideMatrix(1,21)+=(crLeftHandSideMatrix10*rB(4,21) + crLeftHandSideMatrix11*rB(5,21) + crLeftHandSideMatrix6*rB(0,21) + crLeftHandSideMatrix7*rB(1,21) + crLeftHandSideMatrix8*rB(2,21) + crLeftHandSideMatrix9*rB(3,21));
// rLeftHandSideMatrix(1,22)+=(crLeftHandSideMatrix10*rB(4,22) + crLeftHandSideMatrix11*rB(5,22) + crLeftHandSideMatrix6*rB(0,22) + crLeftHandSideMatrix7*rB(1,22) + crLeftHandSideMatrix8*rB(2,22) + crLeftHandSideMatrix9*rB(3,22));
// rLeftHandSideMatrix(1,23)+=(crLeftHandSideMatrix10*rB(4,23) + crLeftHandSideMatrix11*rB(5,23) + crLeftHandSideMatrix6*rB(0,23) + crLeftHandSideMatrix7*rB(1,23) + crLeftHandSideMatrix8*rB(2,23) + crLeftHandSideMatrix9*rB(3,23));

// // rLeftHandSideMatrix(2,0)+=(crLeftHandSideMatrix12*rB(0,0) + crLeftHandSideMatrix13*rB(1,0) + crLeftHandSideMatrix14*rB(2,0) + crLeftHandSideMatrix15*rB(3,0) + crLeftHandSideMatrix16*rB(4,0) + crLeftHandSideMatrix17*rB(5,0));
// // rLeftHandSideMatrix(2,1)+=(crLeftHandSideMatrix12*rB(0,1) + crLeftHandSideMatrix13*rB(1,1) + crLeftHandSideMatrix14*rB(2,1) + crLeftHandSideMatrix15*rB(3,1) + crLeftHandSideMatrix16*rB(4,1) + crLeftHandSideMatrix17*rB(5,1));
// rLeftHandSideMatrix(2,2)+=(crLeftHandSideMatrix12*rB(0,2) + crLeftHandSideMatrix13*rB(1,2) + crLeftHandSideMatrix14*rB(2,2) + crLeftHandSideMatrix15*rB(3,2) + crLeftHandSideMatrix16*rB(4,2) + crLeftHandSideMatrix17*rB(5,2));
// rLeftHandSideMatrix(2,3)+=(crLeftHandSideMatrix12*rB(0,3) + crLeftHandSideMatrix13*rB(1,3) + crLeftHandSideMatrix14*rB(2,3) + crLeftHandSideMatrix15*rB(3,3) + crLeftHandSideMatrix16*rB(4,3) + crLeftHandSideMatrix17*rB(5,3));
// rLeftHandSideMatrix(2,4)+=(crLeftHandSideMatrix12*rB(0,4) + crLeftHandSideMatrix13*rB(1,4) + crLeftHandSideMatrix14*rB(2,4) + crLeftHandSideMatrix15*rB(3,4) + crLeftHandSideMatrix16*rB(4,4) + crLeftHandSideMatrix17*rB(5,4));
// rLeftHandSideMatrix(2,5)+=(crLeftHandSideMatrix12*rB(0,5) + crLeftHandSideMatrix13*rB(1,5) + crLeftHandSideMatrix14*rB(2,5) + crLeftHandSideMatrix15*rB(3,5) + crLeftHandSideMatrix16*rB(4,5) + crLeftHandSideMatrix17*rB(5,5));
// rLeftHandSideMatrix(2,6)+=(crLeftHandSideMatrix12*rB(0,6) + crLeftHandSideMatrix13*rB(1,6) + crLeftHandSideMatrix14*rB(2,6) + crLeftHandSideMatrix15*rB(3,6) + crLeftHandSideMatrix16*rB(4,6) + crLeftHandSideMatrix17*rB(5,6));
// rLeftHandSideMatrix(2,7)+=(crLeftHandSideMatrix12*rB(0,7) + crLeftHandSideMatrix13*rB(1,7) + crLeftHandSideMatrix14*rB(2,7) + crLeftHandSideMatrix15*rB(3,7) + crLeftHandSideMatrix16*rB(4,7) + crLeftHandSideMatrix17*rB(5,7));
// rLeftHandSideMatrix(2,8)+=(crLeftHandSideMatrix12*rB(0,8) + crLeftHandSideMatrix13*rB(1,8) + crLeftHandSideMatrix14*rB(2,8) + crLeftHandSideMatrix15*rB(3,8) + crLeftHandSideMatrix16*rB(4,8) + crLeftHandSideMatrix17*rB(5,8));
// rLeftHandSideMatrix(2,9)+=(crLeftHandSideMatrix12*rB(0,9) + crLeftHandSideMatrix13*rB(1,9) + crLeftHandSideMatrix14*rB(2,9) + crLeftHandSideMatrix15*rB(3,9) + crLeftHandSideMatrix16*rB(4,9) + crLeftHandSideMatrix17*rB(5,9));
// rLeftHandSideMatrix(2,10)+=(crLeftHandSideMatrix12*rB(0,10) + crLeftHandSideMatrix13*rB(1,10) + crLeftHandSideMatrix14*rB(2,10) + crLeftHandSideMatrix15*rB(3,10) + crLeftHandSideMatrix16*rB(4,10) + crLeftHandSideMatrix17*rB(5,10));
// rLeftHandSideMatrix(2,11)+=(crLeftHandSideMatrix12*rB(0,11) + crLeftHandSideMatrix13*rB(1,11) + crLeftHandSideMatrix14*rB(2,11) + crLeftHandSideMatrix15*rB(3,11) + crLeftHandSideMatrix16*rB(4,11) + crLeftHandSideMatrix17*rB(5,11));
// rLeftHandSideMatrix(2,12)+=(crLeftHandSideMatrix12*rB(0,12) + crLeftHandSideMatrix13*rB(1,12) + crLeftHandSideMatrix14*rB(2,12) + crLeftHandSideMatrix15*rB(3,12) + crLeftHandSideMatrix16*rB(4,12) + crLeftHandSideMatrix17*rB(5,12));
// rLeftHandSideMatrix(2,13)+=(crLeftHandSideMatrix12*rB(0,13) + crLeftHandSideMatrix13*rB(1,13) + crLeftHandSideMatrix14*rB(2,13) + crLeftHandSideMatrix15*rB(3,13) + crLeftHandSideMatrix16*rB(4,13) + crLeftHandSideMatrix17*rB(5,13));
// rLeftHandSideMatrix(2,14)+=(crLeftHandSideMatrix12*rB(0,14) + crLeftHandSideMatrix13*rB(1,14) + crLeftHandSideMatrix14*rB(2,14) + crLeftHandSideMatrix15*rB(3,14) + crLeftHandSideMatrix16*rB(4,14) + crLeftHandSideMatrix17*rB(5,14));
// rLeftHandSideMatrix(2,15)+=(crLeftHandSideMatrix12*rB(0,15) + crLeftHandSideMatrix13*rB(1,15) + crLeftHandSideMatrix14*rB(2,15) + crLeftHandSideMatrix15*rB(3,15) + crLeftHandSideMatrix16*rB(4,15) + crLeftHandSideMatrix17*rB(5,15));
// rLeftHandSideMatrix(2,16)+=(crLeftHandSideMatrix12*rB(0,16) + crLeftHandSideMatrix13*rB(1,16) + crLeftHandSideMatrix14*rB(2,16) + crLeftHandSideMatrix15*rB(3,16) + crLeftHandSideMatrix16*rB(4,16) + crLeftHandSideMatrix17*rB(5,16));
// rLeftHandSideMatrix(2,17)+=(crLeftHandSideMatrix12*rB(0,17) + crLeftHandSideMatrix13*rB(1,17) + crLeftHandSideMatrix14*rB(2,17) + crLeftHandSideMatrix15*rB(3,17) + crLeftHandSideMatrix16*rB(4,17) + crLeftHandSideMatrix17*rB(5,17));
// rLeftHandSideMatrix(2,18)+=(crLeftHandSideMatrix12*rB(0,18) + crLeftHandSideMatrix13*rB(1,18) + crLeftHandSideMatrix14*rB(2,18) + crLeftHandSideMatrix15*rB(3,18) + crLeftHandSideMatrix16*rB(4,18) + crLeftHandSideMatrix17*rB(5,18));
// rLeftHandSideMatrix(2,19)+=(crLeftHandSideMatrix12*rB(0,19) + crLeftHandSideMatrix13*rB(1,19) + crLeftHandSideMatrix14*rB(2,19) + crLeftHandSideMatrix15*rB(3,19) + crLeftHandSideMatrix16*rB(4,19) + crLeftHandSideMatrix17*rB(5,19));
// rLeftHandSideMatrix(2,20)+=(crLeftHandSideMatrix12*rB(0,20) + crLeftHandSideMatrix13*rB(1,20) + crLeftHandSideMatrix14*rB(2,20) + crLeftHandSideMatrix15*rB(3,20) + crLeftHandSideMatrix16*rB(4,20) + crLeftHandSideMatrix17*rB(5,20));
// rLeftHandSideMatrix(2,21)+=(crLeftHandSideMatrix12*rB(0,21) + crLeftHandSideMatrix13*rB(1,21) + crLeftHandSideMatrix14*rB(2,21) + crLeftHandSideMatrix15*rB(3,21) + crLeftHandSideMatrix16*rB(4,21) + crLeftHandSideMatrix17*rB(5,21));
// rLeftHandSideMatrix(2,22)+=(crLeftHandSideMatrix12*rB(0,22) + crLeftHandSideMatrix13*rB(1,22) + crLeftHandSideMatrix14*rB(2,22) + crLeftHandSideMatrix15*rB(3,22) + crLeftHandSideMatrix16*rB(4,22) + crLeftHandSideMatrix17*rB(5,22));
// rLeftHandSideMatrix(2,23)+=(crLeftHandSideMatrix12*rB(0,23) + crLeftHandSideMatrix13*rB(1,23) + crLeftHandSideMatrix14*rB(2,23) + crLeftHandSideMatrix15*rB(3,23) + crLeftHandSideMatrix16*rB(4,23) + crLeftHandSideMatrix17*rB(5,23));

// // rLeftHandSideMatrix(3,0)+=(crLeftHandSideMatrix18*rB(0,0) + crLeftHandSideMatrix19*rB(1,0) + crLeftHandSideMatrix20*rB(2,0) + crLeftHandSideMatrix21*rB(3,0) + crLeftHandSideMatrix22*rB(4,0) + crLeftHandSideMatrix23*rB(5,0));
// // rLeftHandSideMatrix(3,1)+=(crLeftHandSideMatrix18*rB(0,1) + crLeftHandSideMatrix19*rB(1,1) + crLeftHandSideMatrix20*rB(2,1) + crLeftHandSideMatrix21*rB(3,1) + crLeftHandSideMatrix22*rB(4,1) + crLeftHandSideMatrix23*rB(5,1));
// // rLeftHandSideMatrix(3,2)+=(crLeftHandSideMatrix18*rB(0,2) + crLeftHandSideMatrix19*rB(1,2) + crLeftHandSideMatrix20*rB(2,2) + crLeftHandSideMatrix21*rB(3,2) + crLeftHandSideMatrix22*rB(4,2) + crLeftHandSideMatrix23*rB(5,2));
// rLeftHandSideMatrix(3,3)+=(crLeftHandSideMatrix18*rB(0,3) + crLeftHandSideMatrix19*rB(1,3) + crLeftHandSideMatrix20*rB(2,3) + crLeftHandSideMatrix21*rB(3,3) + crLeftHandSideMatrix22*rB(4,3) + crLeftHandSideMatrix23*rB(5,3));
// rLeftHandSideMatrix(3,4)+=(crLeftHandSideMatrix18*rB(0,4) + crLeftHandSideMatrix19*rB(1,4) + crLeftHandSideMatrix20*rB(2,4) + crLeftHandSideMatrix21*rB(3,4) + crLeftHandSideMatrix22*rB(4,4) + crLeftHandSideMatrix23*rB(5,4));
// rLeftHandSideMatrix(3,5)+=(crLeftHandSideMatrix18*rB(0,5) + crLeftHandSideMatrix19*rB(1,5) + crLeftHandSideMatrix20*rB(2,5) + crLeftHandSideMatrix21*rB(3,5) + crLeftHandSideMatrix22*rB(4,5) + crLeftHandSideMatrix23*rB(5,5));
// rLeftHandSideMatrix(3,6)+=(crLeftHandSideMatrix18*rB(0,6) + crLeftHandSideMatrix19*rB(1,6) + crLeftHandSideMatrix20*rB(2,6) + crLeftHandSideMatrix21*rB(3,6) + crLeftHandSideMatrix22*rB(4,6) + crLeftHandSideMatrix23*rB(5,6));
// rLeftHandSideMatrix(3,7)+=(crLeftHandSideMatrix18*rB(0,7) + crLeftHandSideMatrix19*rB(1,7) + crLeftHandSideMatrix20*rB(2,7) + crLeftHandSideMatrix21*rB(3,7) + crLeftHandSideMatrix22*rB(4,7) + crLeftHandSideMatrix23*rB(5,7));
// rLeftHandSideMatrix(3,8)+=(crLeftHandSideMatrix18*rB(0,8) + crLeftHandSideMatrix19*rB(1,8) + crLeftHandSideMatrix20*rB(2,8) + crLeftHandSideMatrix21*rB(3,8) + crLeftHandSideMatrix22*rB(4,8) + crLeftHandSideMatrix23*rB(5,8));
// rLeftHandSideMatrix(3,9)+=(crLeftHandSideMatrix18*rB(0,9) + crLeftHandSideMatrix19*rB(1,9) + crLeftHandSideMatrix20*rB(2,9) + crLeftHandSideMatrix21*rB(3,9) + crLeftHandSideMatrix22*rB(4,9) + crLeftHandSideMatrix23*rB(5,9));
// rLeftHandSideMatrix(3,10)+=(crLeftHandSideMatrix18*rB(0,10) + crLeftHandSideMatrix19*rB(1,10) + crLeftHandSideMatrix20*rB(2,10) + crLeftHandSideMatrix21*rB(3,10) + crLeftHandSideMatrix22*rB(4,10) + crLeftHandSideMatrix23*rB(5,10));
// rLeftHandSideMatrix(3,11)+=(crLeftHandSideMatrix18*rB(0,11) + crLeftHandSideMatrix19*rB(1,11) + crLeftHandSideMatrix20*rB(2,11) + crLeftHandSideMatrix21*rB(3,11) + crLeftHandSideMatrix22*rB(4,11) + crLeftHandSideMatrix23*rB(5,11));
// rLeftHandSideMatrix(3,12)+=(crLeftHandSideMatrix18*rB(0,12) + crLeftHandSideMatrix19*rB(1,12) + crLeftHandSideMatrix20*rB(2,12) + crLeftHandSideMatrix21*rB(3,12) + crLeftHandSideMatrix22*rB(4,12) + crLeftHandSideMatrix23*rB(5,12));
// rLeftHandSideMatrix(3,13)+=(crLeftHandSideMatrix18*rB(0,13) + crLeftHandSideMatrix19*rB(1,13) + crLeftHandSideMatrix20*rB(2,13) + crLeftHandSideMatrix21*rB(3,13) + crLeftHandSideMatrix22*rB(4,13) + crLeftHandSideMatrix23*rB(5,13));
// rLeftHandSideMatrix(3,14)+=(crLeftHandSideMatrix18*rB(0,14) + crLeftHandSideMatrix19*rB(1,14) + crLeftHandSideMatrix20*rB(2,14) + crLeftHandSideMatrix21*rB(3,14) + crLeftHandSideMatrix22*rB(4,14) + crLeftHandSideMatrix23*rB(5,14));
// rLeftHandSideMatrix(3,15)+=(crLeftHandSideMatrix18*rB(0,15) + crLeftHandSideMatrix19*rB(1,15) + crLeftHandSideMatrix20*rB(2,15) + crLeftHandSideMatrix21*rB(3,15) + crLeftHandSideMatrix22*rB(4,15) + crLeftHandSideMatrix23*rB(5,15));
// rLeftHandSideMatrix(3,16)+=(crLeftHandSideMatrix18*rB(0,16) + crLeftHandSideMatrix19*rB(1,16) + crLeftHandSideMatrix20*rB(2,16) + crLeftHandSideMatrix21*rB(3,16) + crLeftHandSideMatrix22*rB(4,16) + crLeftHandSideMatrix23*rB(5,16));
// rLeftHandSideMatrix(3,17)+=(crLeftHandSideMatrix18*rB(0,17) + crLeftHandSideMatrix19*rB(1,17) + crLeftHandSideMatrix20*rB(2,17) + crLeftHandSideMatrix21*rB(3,17) + crLeftHandSideMatrix22*rB(4,17) + crLeftHandSideMatrix23*rB(5,17));
// rLeftHandSideMatrix(3,18)+=(crLeftHandSideMatrix18*rB(0,18) + crLeftHandSideMatrix19*rB(1,18) + crLeftHandSideMatrix20*rB(2,18) + crLeftHandSideMatrix21*rB(3,18) + crLeftHandSideMatrix22*rB(4,18) + crLeftHandSideMatrix23*rB(5,18));
// rLeftHandSideMatrix(3,19)+=(crLeftHandSideMatrix18*rB(0,19) + crLeftHandSideMatrix19*rB(1,19) + crLeftHandSideMatrix20*rB(2,19) + crLeftHandSideMatrix21*rB(3,19) + crLeftHandSideMatrix22*rB(4,19) + crLeftHandSideMatrix23*rB(5,19));
// rLeftHandSideMatrix(3,20)+=(crLeftHandSideMatrix18*rB(0,20) + crLeftHandSideMatrix19*rB(1,20) + crLeftHandSideMatrix20*rB(2,20) + crLeftHandSideMatrix21*rB(3,20) + crLeftHandSideMatrix22*rB(4,20) + crLeftHandSideMatrix23*rB(5,20));
// rLeftHandSideMatrix(3,21)+=(crLeftHandSideMatrix18*rB(0,21) + crLeftHandSideMatrix19*rB(1,21) + crLeftHandSideMatrix20*rB(2,21) + crLeftHandSideMatrix21*rB(3,21) + crLeftHandSideMatrix22*rB(4,21) + crLeftHandSideMatrix23*rB(5,21));
// rLeftHandSideMatrix(3,22)+=(crLeftHandSideMatrix18*rB(0,22) + crLeftHandSideMatrix19*rB(1,22) + crLeftHandSideMatrix20*rB(2,22) + crLeftHandSideMatrix21*rB(3,22) + crLeftHandSideMatrix22*rB(4,22) + crLeftHandSideMatrix23*rB(5,22));
// rLeftHandSideMatrix(3,23)+=(crLeftHandSideMatrix18*rB(0,23) + crLeftHandSideMatrix19*rB(1,23) + crLeftHandSideMatrix20*rB(2,23) + crLeftHandSideMatrix21*rB(3,23) + crLeftHandSideMatrix22*rB(4,23) + crLeftHandSideMatrix23*rB(5,23));

// // rLeftHandSideMatrix(4,0)+=(crLeftHandSideMatrix24*rB(0,0) + crLeftHandSideMatrix25*rB(1,0) + crLeftHandSideMatrix26*rB(2,0) + crLeftHandSideMatrix27*rB(3,0) + crLeftHandSideMatrix28*rB(4,0) + crLeftHandSideMatrix29*rB(5,0));
// // rLeftHandSideMatrix(4,1)+=(crLeftHandSideMatrix24*rB(0,1) + crLeftHandSideMatrix25*rB(1,1) + crLeftHandSideMatrix26*rB(2,1) + crLeftHandSideMatrix27*rB(3,1) + crLeftHandSideMatrix28*rB(4,1) + crLeftHandSideMatrix29*rB(5,1));
// // rLeftHandSideMatrix(4,2)+=(crLeftHandSideMatrix24*rB(0,2) + crLeftHandSideMatrix25*rB(1,2) + crLeftHandSideMatrix26*rB(2,2) + crLeftHandSideMatrix27*rB(3,2) + crLeftHandSideMatrix28*rB(4,2) + crLeftHandSideMatrix29*rB(5,2));
// // rLeftHandSideMatrix(4,3)+=(crLeftHandSideMatrix24*rB(0,3) + crLeftHandSideMatrix25*rB(1,3) + crLeftHandSideMatrix26*rB(2,3) + crLeftHandSideMatrix27*rB(3,3) + crLeftHandSideMatrix28*rB(4,3) + crLeftHandSideMatrix29*rB(5,3));
// rLeftHandSideMatrix(4,4)+=(crLeftHandSideMatrix24*rB(0,4) + crLeftHandSideMatrix25*rB(1,4) + crLeftHandSideMatrix26*rB(2,4) + crLeftHandSideMatrix27*rB(3,4) + crLeftHandSideMatrix28*rB(4,4) + crLeftHandSideMatrix29*rB(5,4));
// rLeftHandSideMatrix(4,5)+=(crLeftHandSideMatrix24*rB(0,5) + crLeftHandSideMatrix25*rB(1,5) + crLeftHandSideMatrix26*rB(2,5) + crLeftHandSideMatrix27*rB(3,5) + crLeftHandSideMatrix28*rB(4,5) + crLeftHandSideMatrix29*rB(5,5));
// rLeftHandSideMatrix(4,6)+=(crLeftHandSideMatrix24*rB(0,6) + crLeftHandSideMatrix25*rB(1,6) + crLeftHandSideMatrix26*rB(2,6) + crLeftHandSideMatrix27*rB(3,6) + crLeftHandSideMatrix28*rB(4,6) + crLeftHandSideMatrix29*rB(5,6));
// rLeftHandSideMatrix(4,7)+=(crLeftHandSideMatrix24*rB(0,7) + crLeftHandSideMatrix25*rB(1,7) + crLeftHandSideMatrix26*rB(2,7) + crLeftHandSideMatrix27*rB(3,7) + crLeftHandSideMatrix28*rB(4,7) + crLeftHandSideMatrix29*rB(5,7));
// rLeftHandSideMatrix(4,8)+=(crLeftHandSideMatrix24*rB(0,8) + crLeftHandSideMatrix25*rB(1,8) + crLeftHandSideMatrix26*rB(2,8) + crLeftHandSideMatrix27*rB(3,8) + crLeftHandSideMatrix28*rB(4,8) + crLeftHandSideMatrix29*rB(5,8));
// rLeftHandSideMatrix(4,9)+=(crLeftHandSideMatrix24*rB(0,9) + crLeftHandSideMatrix25*rB(1,9) + crLeftHandSideMatrix26*rB(2,9) + crLeftHandSideMatrix27*rB(3,9) + crLeftHandSideMatrix28*rB(4,9) + crLeftHandSideMatrix29*rB(5,9));
// rLeftHandSideMatrix(4,10)+=(crLeftHandSideMatrix24*rB(0,10) + crLeftHandSideMatrix25*rB(1,10) + crLeftHandSideMatrix26*rB(2,10) + crLeftHandSideMatrix27*rB(3,10) + crLeftHandSideMatrix28*rB(4,10) + crLeftHandSideMatrix29*rB(5,10));
// rLeftHandSideMatrix(4,11)+=(crLeftHandSideMatrix24*rB(0,11) + crLeftHandSideMatrix25*rB(1,11) + crLeftHandSideMatrix26*rB(2,11) + crLeftHandSideMatrix27*rB(3,11) + crLeftHandSideMatrix28*rB(4,11) + crLeftHandSideMatrix29*rB(5,11));
// rLeftHandSideMatrix(4,12)+=(crLeftHandSideMatrix24*rB(0,12) + crLeftHandSideMatrix25*rB(1,12) + crLeftHandSideMatrix26*rB(2,12) + crLeftHandSideMatrix27*rB(3,12) + crLeftHandSideMatrix28*rB(4,12) + crLeftHandSideMatrix29*rB(5,12));
// rLeftHandSideMatrix(4,13)+=(crLeftHandSideMatrix24*rB(0,13) + crLeftHandSideMatrix25*rB(1,13) + crLeftHandSideMatrix26*rB(2,13) + crLeftHandSideMatrix27*rB(3,13) + crLeftHandSideMatrix28*rB(4,13) + crLeftHandSideMatrix29*rB(5,13));
// rLeftHandSideMatrix(4,14)+=(crLeftHandSideMatrix24*rB(0,14) + crLeftHandSideMatrix25*rB(1,14) + crLeftHandSideMatrix26*rB(2,14) + crLeftHandSideMatrix27*rB(3,14) + crLeftHandSideMatrix28*rB(4,14) + crLeftHandSideMatrix29*rB(5,14));
// rLeftHandSideMatrix(4,15)+=(crLeftHandSideMatrix24*rB(0,15) + crLeftHandSideMatrix25*rB(1,15) + crLeftHandSideMatrix26*rB(2,15) + crLeftHandSideMatrix27*rB(3,15) + crLeftHandSideMatrix28*rB(4,15) + crLeftHandSideMatrix29*rB(5,15));
// rLeftHandSideMatrix(4,16)+=(crLeftHandSideMatrix24*rB(0,16) + crLeftHandSideMatrix25*rB(1,16) + crLeftHandSideMatrix26*rB(2,16) + crLeftHandSideMatrix27*rB(3,16) + crLeftHandSideMatrix28*rB(4,16) + crLeftHandSideMatrix29*rB(5,16));
// rLeftHandSideMatrix(4,17)+=(crLeftHandSideMatrix24*rB(0,17) + crLeftHandSideMatrix25*rB(1,17) + crLeftHandSideMatrix26*rB(2,17) + crLeftHandSideMatrix27*rB(3,17) + crLeftHandSideMatrix28*rB(4,17) + crLeftHandSideMatrix29*rB(5,17));
// rLeftHandSideMatrix(4,18)+=(crLeftHandSideMatrix24*rB(0,18) + crLeftHandSideMatrix25*rB(1,18) + crLeftHandSideMatrix26*rB(2,18) + crLeftHandSideMatrix27*rB(3,18) + crLeftHandSideMatrix28*rB(4,18) + crLeftHandSideMatrix29*rB(5,18));
// rLeftHandSideMatrix(4,19)+=(crLeftHandSideMatrix24*rB(0,19) + crLeftHandSideMatrix25*rB(1,19) + crLeftHandSideMatrix26*rB(2,19) + crLeftHandSideMatrix27*rB(3,19) + crLeftHandSideMatrix28*rB(4,19) + crLeftHandSideMatrix29*rB(5,19));
// rLeftHandSideMatrix(4,20)+=(crLeftHandSideMatrix24*rB(0,20) + crLeftHandSideMatrix25*rB(1,20) + crLeftHandSideMatrix26*rB(2,20) + crLeftHandSideMatrix27*rB(3,20) + crLeftHandSideMatrix28*rB(4,20) + crLeftHandSideMatrix29*rB(5,20));
// rLeftHandSideMatrix(4,21)+=(crLeftHandSideMatrix24*rB(0,21) + crLeftHandSideMatrix25*rB(1,21) + crLeftHandSideMatrix26*rB(2,21) + crLeftHandSideMatrix27*rB(3,21) + crLeftHandSideMatrix28*rB(4,21) + crLeftHandSideMatrix29*rB(5,21));
// rLeftHandSideMatrix(4,22)+=(crLeftHandSideMatrix24*rB(0,22) + crLeftHandSideMatrix25*rB(1,22) + crLeftHandSideMatrix26*rB(2,22) + crLeftHandSideMatrix27*rB(3,22) + crLeftHandSideMatrix28*rB(4,22) + crLeftHandSideMatrix29*rB(5,22));
// rLeftHandSideMatrix(4,23)+=(crLeftHandSideMatrix24*rB(0,23) + crLeftHandSideMatrix25*rB(1,23) + crLeftHandSideMatrix26*rB(2,23) + crLeftHandSideMatrix27*rB(3,23) + crLeftHandSideMatrix28*rB(4,23) + crLeftHandSideMatrix29*rB(5,23));

// // rLeftHandSideMatrix(5,0)+=(crLeftHandSideMatrix30*rB(0,0) + crLeftHandSideMatrix31*rB(1,0) + crLeftHandSideMatrix32*rB(2,0) + crLeftHandSideMatrix33*rB(3,0) + crLeftHandSideMatrix34*rB(4,0) + crLeftHandSideMatrix35*rB(5,0));
// // rLeftHandSideMatrix(5,1)+=(crLeftHandSideMatrix30*rB(0,1) + crLeftHandSideMatrix31*rB(1,1) + crLeftHandSideMatrix32*rB(2,1) + crLeftHandSideMatrix33*rB(3,1) + crLeftHandSideMatrix34*rB(4,1) + crLeftHandSideMatrix35*rB(5,1));
// // rLeftHandSideMatrix(5,2)+=(crLeftHandSideMatrix30*rB(0,2) + crLeftHandSideMatrix31*rB(1,2) + crLeftHandSideMatrix32*rB(2,2) + crLeftHandSideMatrix33*rB(3,2) + crLeftHandSideMatrix34*rB(4,2) + crLeftHandSideMatrix35*rB(5,2));
// // rLeftHandSideMatrix(5,3)+=(crLeftHandSideMatrix30*rB(0,3) + crLeftHandSideMatrix31*rB(1,3) + crLeftHandSideMatrix32*rB(2,3) + crLeftHandSideMatrix33*rB(3,3) + crLeftHandSideMatrix34*rB(4,3) + crLeftHandSideMatrix35*rB(5,3));
// // rLeftHandSideMatrix(5,4)+=(crLeftHandSideMatrix30*rB(0,4) + crLeftHandSideMatrix31*rB(1,4) + crLeftHandSideMatrix32*rB(2,4) + crLeftHandSideMatrix33*rB(3,4) + crLeftHandSideMatrix34*rB(4,4) + crLeftHandSideMatrix35*rB(5,4));
// rLeftHandSideMatrix(5,5)+=(crLeftHandSideMatrix30*rB(0,5) + crLeftHandSideMatrix31*rB(1,5) + crLeftHandSideMatrix32*rB(2,5) + crLeftHandSideMatrix33*rB(3,5) + crLeftHandSideMatrix34*rB(4,5) + crLeftHandSideMatrix35*rB(5,5));
// rLeftHandSideMatrix(5,6)+=(crLeftHandSideMatrix30*rB(0,6) + crLeftHandSideMatrix31*rB(1,6) + crLeftHandSideMatrix32*rB(2,6) + crLeftHandSideMatrix33*rB(3,6) + crLeftHandSideMatrix34*rB(4,6) + crLeftHandSideMatrix35*rB(5,6));
// rLeftHandSideMatrix(5,7)+=(crLeftHandSideMatrix30*rB(0,7) + crLeftHandSideMatrix31*rB(1,7) + crLeftHandSideMatrix32*rB(2,7) + crLeftHandSideMatrix33*rB(3,7) + crLeftHandSideMatrix34*rB(4,7) + crLeftHandSideMatrix35*rB(5,7));
// rLeftHandSideMatrix(5,8)+=(crLeftHandSideMatrix30*rB(0,8) + crLeftHandSideMatrix31*rB(1,8) + crLeftHandSideMatrix32*rB(2,8) + crLeftHandSideMatrix33*rB(3,8) + crLeftHandSideMatrix34*rB(4,8) + crLeftHandSideMatrix35*rB(5,8));
// rLeftHandSideMatrix(5,9)+=(crLeftHandSideMatrix30*rB(0,9) + crLeftHandSideMatrix31*rB(1,9) + crLeftHandSideMatrix32*rB(2,9) + crLeftHandSideMatrix33*rB(3,9) + crLeftHandSideMatrix34*rB(4,9) + crLeftHandSideMatrix35*rB(5,9));
// rLeftHandSideMatrix(5,10)+=(crLeftHandSideMatrix30*rB(0,10) + crLeftHandSideMatrix31*rB(1,10) + crLeftHandSideMatrix32*rB(2,10) + crLeftHandSideMatrix33*rB(3,10) + crLeftHandSideMatrix34*rB(4,10) + crLeftHandSideMatrix35*rB(5,10));
// rLeftHandSideMatrix(5,11)+=(crLeftHandSideMatrix30*rB(0,11) + crLeftHandSideMatrix31*rB(1,11) + crLeftHandSideMatrix32*rB(2,11) + crLeftHandSideMatrix33*rB(3,11) + crLeftHandSideMatrix34*rB(4,11) + crLeftHandSideMatrix35*rB(5,11));
// rLeftHandSideMatrix(5,12)+=(crLeftHandSideMatrix30*rB(0,12) + crLeftHandSideMatrix31*rB(1,12) + crLeftHandSideMatrix32*rB(2,12) + crLeftHandSideMatrix33*rB(3,12) + crLeftHandSideMatrix34*rB(4,12) + crLeftHandSideMatrix35*rB(5,12));
// rLeftHandSideMatrix(5,13)+=(crLeftHandSideMatrix30*rB(0,13) + crLeftHandSideMatrix31*rB(1,13) + crLeftHandSideMatrix32*rB(2,13) + crLeftHandSideMatrix33*rB(3,13) + crLeftHandSideMatrix34*rB(4,13) + crLeftHandSideMatrix35*rB(5,13));
// rLeftHandSideMatrix(5,14)+=(crLeftHandSideMatrix30*rB(0,14) + crLeftHandSideMatrix31*rB(1,14) + crLeftHandSideMatrix32*rB(2,14) + crLeftHandSideMatrix33*rB(3,14) + crLeftHandSideMatrix34*rB(4,14) + crLeftHandSideMatrix35*rB(5,14));
// rLeftHandSideMatrix(5,15)+=(crLeftHandSideMatrix30*rB(0,15) + crLeftHandSideMatrix31*rB(1,15) + crLeftHandSideMatrix32*rB(2,15) + crLeftHandSideMatrix33*rB(3,15) + crLeftHandSideMatrix34*rB(4,15) + crLeftHandSideMatrix35*rB(5,15));
// rLeftHandSideMatrix(5,16)+=(crLeftHandSideMatrix30*rB(0,16) + crLeftHandSideMatrix31*rB(1,16) + crLeftHandSideMatrix32*rB(2,16) + crLeftHandSideMatrix33*rB(3,16) + crLeftHandSideMatrix34*rB(4,16) + crLeftHandSideMatrix35*rB(5,16));
// rLeftHandSideMatrix(5,17)+=(crLeftHandSideMatrix30*rB(0,17) + crLeftHandSideMatrix31*rB(1,17) + crLeftHandSideMatrix32*rB(2,17) + crLeftHandSideMatrix33*rB(3,17) + crLeftHandSideMatrix34*rB(4,17) + crLeftHandSideMatrix35*rB(5,17));
// rLeftHandSideMatrix(5,18)+=(crLeftHandSideMatrix30*rB(0,18) + crLeftHandSideMatrix31*rB(1,18) + crLeftHandSideMatrix32*rB(2,18) + crLeftHandSideMatrix33*rB(3,18) + crLeftHandSideMatrix34*rB(4,18) + crLeftHandSideMatrix35*rB(5,18));
// rLeftHandSideMatrix(5,19)+=(crLeftHandSideMatrix30*rB(0,19) + crLeftHandSideMatrix31*rB(1,19) + crLeftHandSideMatrix32*rB(2,19) + crLeftHandSideMatrix33*rB(3,19) + crLeftHandSideMatrix34*rB(4,19) + crLeftHandSideMatrix35*rB(5,19));
// rLeftHandSideMatrix(5,20)+=(crLeftHandSideMatrix30*rB(0,20) + crLeftHandSideMatrix31*rB(1,20) + crLeftHandSideMatrix32*rB(2,20) + crLeftHandSideMatrix33*rB(3,20) + crLeftHandSideMatrix34*rB(4,20) + crLeftHandSideMatrix35*rB(5,20));
// rLeftHandSideMatrix(5,21)+=(crLeftHandSideMatrix30*rB(0,21) + crLeftHandSideMatrix31*rB(1,21) + crLeftHandSideMatrix32*rB(2,21) + crLeftHandSideMatrix33*rB(3,21) + crLeftHandSideMatrix34*rB(4,21) + crLeftHandSideMatrix35*rB(5,21));
// rLeftHandSideMatrix(5,22)+=(crLeftHandSideMatrix30*rB(0,22) + crLeftHandSideMatrix31*rB(1,22) + crLeftHandSideMatrix32*rB(2,22) + crLeftHandSideMatrix33*rB(3,22) + crLeftHandSideMatrix34*rB(4,22) + crLeftHandSideMatrix35*rB(5,22));
// rLeftHandSideMatrix(5,23)+=(crLeftHandSideMatrix30*rB(0,23) + crLeftHandSideMatrix31*rB(1,23) + crLeftHandSideMatrix32*rB(2,23) + crLeftHandSideMatrix33*rB(3,23) + crLeftHandSideMatrix34*rB(4,23) + crLeftHandSideMatrix35*rB(5,23));

// // rLeftHandSideMatrix(6,0)+=(crLeftHandSideMatrix36*rB(0,0) + crLeftHandSideMatrix37*rB(1,0) + crLeftHandSideMatrix38*rB(2,0) + crLeftHandSideMatrix39*rB(3,0) + crLeftHandSideMatrix40*rB(4,0) + crLeftHandSideMatrix41*rB(5,0));
// // rLeftHandSideMatrix(6,1)+=(crLeftHandSideMatrix36*rB(0,1) + crLeftHandSideMatrix37*rB(1,1) + crLeftHandSideMatrix38*rB(2,1) + crLeftHandSideMatrix39*rB(3,1) + crLeftHandSideMatrix40*rB(4,1) + crLeftHandSideMatrix41*rB(5,1));
// // rLeftHandSideMatrix(6,2)+=(crLeftHandSideMatrix36*rB(0,2) + crLeftHandSideMatrix37*rB(1,2) + crLeftHandSideMatrix38*rB(2,2) + crLeftHandSideMatrix39*rB(3,2) + crLeftHandSideMatrix40*rB(4,2) + crLeftHandSideMatrix41*rB(5,2));
// // rLeftHandSideMatrix(6,3)+=(crLeftHandSideMatrix36*rB(0,3) + crLeftHandSideMatrix37*rB(1,3) + crLeftHandSideMatrix38*rB(2,3) + crLeftHandSideMatrix39*rB(3,3) + crLeftHandSideMatrix40*rB(4,3) + crLeftHandSideMatrix41*rB(5,3));
// // rLeftHandSideMatrix(6,4)+=(crLeftHandSideMatrix36*rB(0,4) + crLeftHandSideMatrix37*rB(1,4) + crLeftHandSideMatrix38*rB(2,4) + crLeftHandSideMatrix39*rB(3,4) + crLeftHandSideMatrix40*rB(4,4) + crLeftHandSideMatrix41*rB(5,4));
// // rLeftHandSideMatrix(6,5)+=(crLeftHandSideMatrix36*rB(0,5) + crLeftHandSideMatrix37*rB(1,5) + crLeftHandSideMatrix38*rB(2,5) + crLeftHandSideMatrix39*rB(3,5) + crLeftHandSideMatrix40*rB(4,5) + crLeftHandSideMatrix41*rB(5,5));
// rLeftHandSideMatrix(6,6)+=(crLeftHandSideMatrix36*rB(0,6) + crLeftHandSideMatrix37*rB(1,6) + crLeftHandSideMatrix38*rB(2,6) + crLeftHandSideMatrix39*rB(3,6) + crLeftHandSideMatrix40*rB(4,6) + crLeftHandSideMatrix41*rB(5,6));
// rLeftHandSideMatrix(6,7)+=(crLeftHandSideMatrix36*rB(0,7) + crLeftHandSideMatrix37*rB(1,7) + crLeftHandSideMatrix38*rB(2,7) + crLeftHandSideMatrix39*rB(3,7) + crLeftHandSideMatrix40*rB(4,7) + crLeftHandSideMatrix41*rB(5,7));
// rLeftHandSideMatrix(6,8)+=(crLeftHandSideMatrix36*rB(0,8) + crLeftHandSideMatrix37*rB(1,8) + crLeftHandSideMatrix38*rB(2,8) + crLeftHandSideMatrix39*rB(3,8) + crLeftHandSideMatrix40*rB(4,8) + crLeftHandSideMatrix41*rB(5,8));
// rLeftHandSideMatrix(6,9)+=(crLeftHandSideMatrix36*rB(0,9) + crLeftHandSideMatrix37*rB(1,9) + crLeftHandSideMatrix38*rB(2,9) + crLeftHandSideMatrix39*rB(3,9) + crLeftHandSideMatrix40*rB(4,9) + crLeftHandSideMatrix41*rB(5,9));
// rLeftHandSideMatrix(6,10)+=(crLeftHandSideMatrix36*rB(0,10) + crLeftHandSideMatrix37*rB(1,10) + crLeftHandSideMatrix38*rB(2,10) + crLeftHandSideMatrix39*rB(3,10) + crLeftHandSideMatrix40*rB(4,10) + crLeftHandSideMatrix41*rB(5,10));
// rLeftHandSideMatrix(6,11)+=(crLeftHandSideMatrix36*rB(0,11) + crLeftHandSideMatrix37*rB(1,11) + crLeftHandSideMatrix38*rB(2,11) + crLeftHandSideMatrix39*rB(3,11) + crLeftHandSideMatrix40*rB(4,11) + crLeftHandSideMatrix41*rB(5,11));
// rLeftHandSideMatrix(6,12)+=(crLeftHandSideMatrix36*rB(0,12) + crLeftHandSideMatrix37*rB(1,12) + crLeftHandSideMatrix38*rB(2,12) + crLeftHandSideMatrix39*rB(3,12) + crLeftHandSideMatrix40*rB(4,12) + crLeftHandSideMatrix41*rB(5,12));
// rLeftHandSideMatrix(6,13)+=(crLeftHandSideMatrix36*rB(0,13) + crLeftHandSideMatrix37*rB(1,13) + crLeftHandSideMatrix38*rB(2,13) + crLeftHandSideMatrix39*rB(3,13) + crLeftHandSideMatrix40*rB(4,13) + crLeftHandSideMatrix41*rB(5,13));
// rLeftHandSideMatrix(6,14)+=(crLeftHandSideMatrix36*rB(0,14) + crLeftHandSideMatrix37*rB(1,14) + crLeftHandSideMatrix38*rB(2,14) + crLeftHandSideMatrix39*rB(3,14) + crLeftHandSideMatrix40*rB(4,14) + crLeftHandSideMatrix41*rB(5,14));
// rLeftHandSideMatrix(6,15)+=(crLeftHandSideMatrix36*rB(0,15) + crLeftHandSideMatrix37*rB(1,15) + crLeftHandSideMatrix38*rB(2,15) + crLeftHandSideMatrix39*rB(3,15) + crLeftHandSideMatrix40*rB(4,15) + crLeftHandSideMatrix41*rB(5,15));
// rLeftHandSideMatrix(6,16)+=(crLeftHandSideMatrix36*rB(0,16) + crLeftHandSideMatrix37*rB(1,16) + crLeftHandSideMatrix38*rB(2,16) + crLeftHandSideMatrix39*rB(3,16) + crLeftHandSideMatrix40*rB(4,16) + crLeftHandSideMatrix41*rB(5,16));
// rLeftHandSideMatrix(6,17)+=(crLeftHandSideMatrix36*rB(0,17) + crLeftHandSideMatrix37*rB(1,17) + crLeftHandSideMatrix38*rB(2,17) + crLeftHandSideMatrix39*rB(3,17) + crLeftHandSideMatrix40*rB(4,17) + crLeftHandSideMatrix41*rB(5,17));
// rLeftHandSideMatrix(6,18)+=(crLeftHandSideMatrix36*rB(0,18) + crLeftHandSideMatrix37*rB(1,18) + crLeftHandSideMatrix38*rB(2,18) + crLeftHandSideMatrix39*rB(3,18) + crLeftHandSideMatrix40*rB(4,18) + crLeftHandSideMatrix41*rB(5,18));
// rLeftHandSideMatrix(6,19)+=(crLeftHandSideMatrix36*rB(0,19) + crLeftHandSideMatrix37*rB(1,19) + crLeftHandSideMatrix38*rB(2,19) + crLeftHandSideMatrix39*rB(3,19) + crLeftHandSideMatrix40*rB(4,19) + crLeftHandSideMatrix41*rB(5,19));
// rLeftHandSideMatrix(6,20)+=(crLeftHandSideMatrix36*rB(0,20) + crLeftHandSideMatrix37*rB(1,20) + crLeftHandSideMatrix38*rB(2,20) + crLeftHandSideMatrix39*rB(3,20) + crLeftHandSideMatrix40*rB(4,20) + crLeftHandSideMatrix41*rB(5,20));
// rLeftHandSideMatrix(6,21)+=(crLeftHandSideMatrix36*rB(0,21) + crLeftHandSideMatrix37*rB(1,21) + crLeftHandSideMatrix38*rB(2,21) + crLeftHandSideMatrix39*rB(3,21) + crLeftHandSideMatrix40*rB(4,21) + crLeftHandSideMatrix41*rB(5,21));
// rLeftHandSideMatrix(6,22)+=(crLeftHandSideMatrix36*rB(0,22) + crLeftHandSideMatrix37*rB(1,22) + crLeftHandSideMatrix38*rB(2,22) + crLeftHandSideMatrix39*rB(3,22) + crLeftHandSideMatrix40*rB(4,22) + crLeftHandSideMatrix41*rB(5,22));
// rLeftHandSideMatrix(6,23)+=(crLeftHandSideMatrix36*rB(0,23) + crLeftHandSideMatrix37*rB(1,23) + crLeftHandSideMatrix38*rB(2,23) + crLeftHandSideMatrix39*rB(3,23) + crLeftHandSideMatrix40*rB(4,23) + crLeftHandSideMatrix41*rB(5,23));

// // rLeftHandSideMatrix(7,0)+=(crLeftHandSideMatrix42*rB(0,0) + crLeftHandSideMatrix43*rB(1,0) + crLeftHandSideMatrix44*rB(2,0) + crLeftHandSideMatrix45*rB(3,0) + crLeftHandSideMatrix46*rB(4,0) + crLeftHandSideMatrix47*rB(5,0));
// // rLeftHandSideMatrix(7,1)+=(crLeftHandSideMatrix42*rB(0,1) + crLeftHandSideMatrix43*rB(1,1) + crLeftHandSideMatrix44*rB(2,1) + crLeftHandSideMatrix45*rB(3,1) + crLeftHandSideMatrix46*rB(4,1) + crLeftHandSideMatrix47*rB(5,1));
// // rLeftHandSideMatrix(7,2)+=(crLeftHandSideMatrix42*rB(0,2) + crLeftHandSideMatrix43*rB(1,2) + crLeftHandSideMatrix44*rB(2,2) + crLeftHandSideMatrix45*rB(3,2) + crLeftHandSideMatrix46*rB(4,2) + crLeftHandSideMatrix47*rB(5,2));
// // rLeftHandSideMatrix(7,3)+=(crLeftHandSideMatrix42*rB(0,3) + crLeftHandSideMatrix43*rB(1,3) + crLeftHandSideMatrix44*rB(2,3) + crLeftHandSideMatrix45*rB(3,3) + crLeftHandSideMatrix46*rB(4,3) + crLeftHandSideMatrix47*rB(5,3));
// // rLeftHandSideMatrix(7,4)+=(crLeftHandSideMatrix42*rB(0,4) + crLeftHandSideMatrix43*rB(1,4) + crLeftHandSideMatrix44*rB(2,4) + crLeftHandSideMatrix45*rB(3,4) + crLeftHandSideMatrix46*rB(4,4) + crLeftHandSideMatrix47*rB(5,4));
// // rLeftHandSideMatrix(7,5)+=(crLeftHandSideMatrix42*rB(0,5) + crLeftHandSideMatrix43*rB(1,5) + crLeftHandSideMatrix44*rB(2,5) + crLeftHandSideMatrix45*rB(3,5) + crLeftHandSideMatrix46*rB(4,5) + crLeftHandSideMatrix47*rB(5,5));
// // rLeftHandSideMatrix(7,6)+=(crLeftHandSideMatrix42*rB(0,6) + crLeftHandSideMatrix43*rB(1,6) + crLeftHandSideMatrix44*rB(2,6) + crLeftHandSideMatrix45*rB(3,6) + crLeftHandSideMatrix46*rB(4,6) + crLeftHandSideMatrix47*rB(5,6));
// rLeftHandSideMatrix(7,7)+=(crLeftHandSideMatrix42*rB(0,7) + crLeftHandSideMatrix43*rB(1,7) + crLeftHandSideMatrix44*rB(2,7) + crLeftHandSideMatrix45*rB(3,7) + crLeftHandSideMatrix46*rB(4,7) + crLeftHandSideMatrix47*rB(5,7));
// rLeftHandSideMatrix(7,8)+=(crLeftHandSideMatrix42*rB(0,8) + crLeftHandSideMatrix43*rB(1,8) + crLeftHandSideMatrix44*rB(2,8) + crLeftHandSideMatrix45*rB(3,8) + crLeftHandSideMatrix46*rB(4,8) + crLeftHandSideMatrix47*rB(5,8));
// rLeftHandSideMatrix(7,9)+=(crLeftHandSideMatrix42*rB(0,9) + crLeftHandSideMatrix43*rB(1,9) + crLeftHandSideMatrix44*rB(2,9) + crLeftHandSideMatrix45*rB(3,9) + crLeftHandSideMatrix46*rB(4,9) + crLeftHandSideMatrix47*rB(5,9));
// rLeftHandSideMatrix(7,10)+=(crLeftHandSideMatrix42*rB(0,10) + crLeftHandSideMatrix43*rB(1,10) + crLeftHandSideMatrix44*rB(2,10) + crLeftHandSideMatrix45*rB(3,10) + crLeftHandSideMatrix46*rB(4,10) + crLeftHandSideMatrix47*rB(5,10));
// rLeftHandSideMatrix(7,11)+=(crLeftHandSideMatrix42*rB(0,11) + crLeftHandSideMatrix43*rB(1,11) + crLeftHandSideMatrix44*rB(2,11) + crLeftHandSideMatrix45*rB(3,11) + crLeftHandSideMatrix46*rB(4,11) + crLeftHandSideMatrix47*rB(5,11));
// rLeftHandSideMatrix(7,12)+=(crLeftHandSideMatrix42*rB(0,12) + crLeftHandSideMatrix43*rB(1,12) + crLeftHandSideMatrix44*rB(2,12) + crLeftHandSideMatrix45*rB(3,12) + crLeftHandSideMatrix46*rB(4,12) + crLeftHandSideMatrix47*rB(5,12));
// rLeftHandSideMatrix(7,13)+=(crLeftHandSideMatrix42*rB(0,13) + crLeftHandSideMatrix43*rB(1,13) + crLeftHandSideMatrix44*rB(2,13) + crLeftHandSideMatrix45*rB(3,13) + crLeftHandSideMatrix46*rB(4,13) + crLeftHandSideMatrix47*rB(5,13));
// rLeftHandSideMatrix(7,14)+=(crLeftHandSideMatrix42*rB(0,14) + crLeftHandSideMatrix43*rB(1,14) + crLeftHandSideMatrix44*rB(2,14) + crLeftHandSideMatrix45*rB(3,14) + crLeftHandSideMatrix46*rB(4,14) + crLeftHandSideMatrix47*rB(5,14));
// rLeftHandSideMatrix(7,15)+=(crLeftHandSideMatrix42*rB(0,15) + crLeftHandSideMatrix43*rB(1,15) + crLeftHandSideMatrix44*rB(2,15) + crLeftHandSideMatrix45*rB(3,15) + crLeftHandSideMatrix46*rB(4,15) + crLeftHandSideMatrix47*rB(5,15));
// rLeftHandSideMatrix(7,16)+=(crLeftHandSideMatrix42*rB(0,16) + crLeftHandSideMatrix43*rB(1,16) + crLeftHandSideMatrix44*rB(2,16) + crLeftHandSideMatrix45*rB(3,16) + crLeftHandSideMatrix46*rB(4,16) + crLeftHandSideMatrix47*rB(5,16));
// rLeftHandSideMatrix(7,17)+=(crLeftHandSideMatrix42*rB(0,17) + crLeftHandSideMatrix43*rB(1,17) + crLeftHandSideMatrix44*rB(2,17) + crLeftHandSideMatrix45*rB(3,17) + crLeftHandSideMatrix46*rB(4,17) + crLeftHandSideMatrix47*rB(5,17));
// rLeftHandSideMatrix(7,18)+=(crLeftHandSideMatrix42*rB(0,18) + crLeftHandSideMatrix43*rB(1,18) + crLeftHandSideMatrix44*rB(2,18) + crLeftHandSideMatrix45*rB(3,18) + crLeftHandSideMatrix46*rB(4,18) + crLeftHandSideMatrix47*rB(5,18));
// rLeftHandSideMatrix(7,19)+=(crLeftHandSideMatrix42*rB(0,19) + crLeftHandSideMatrix43*rB(1,19) + crLeftHandSideMatrix44*rB(2,19) + crLeftHandSideMatrix45*rB(3,19) + crLeftHandSideMatrix46*rB(4,19) + crLeftHandSideMatrix47*rB(5,19));
// rLeftHandSideMatrix(7,20)+=(crLeftHandSideMatrix42*rB(0,20) + crLeftHandSideMatrix43*rB(1,20) + crLeftHandSideMatrix44*rB(2,20) + crLeftHandSideMatrix45*rB(3,20) + crLeftHandSideMatrix46*rB(4,20) + crLeftHandSideMatrix47*rB(5,20));
// rLeftHandSideMatrix(7,21)+=(crLeftHandSideMatrix42*rB(0,21) + crLeftHandSideMatrix43*rB(1,21) + crLeftHandSideMatrix44*rB(2,21) + crLeftHandSideMatrix45*rB(3,21) + crLeftHandSideMatrix46*rB(4,21) + crLeftHandSideMatrix47*rB(5,21));
// rLeftHandSideMatrix(7,22)+=(crLeftHandSideMatrix42*rB(0,22) + crLeftHandSideMatrix43*rB(1,22) + crLeftHandSideMatrix44*rB(2,22) + crLeftHandSideMatrix45*rB(3,22) + crLeftHandSideMatrix46*rB(4,22) + crLeftHandSideMatrix47*rB(5,22));
// rLeftHandSideMatrix(7,23)+=(crLeftHandSideMatrix42*rB(0,23) + crLeftHandSideMatrix43*rB(1,23) + crLeftHandSideMatrix44*rB(2,23) + crLeftHandSideMatrix45*rB(3,23) + crLeftHandSideMatrix46*rB(4,23) + crLeftHandSideMatrix47*rB(5,23));

// // rLeftHandSideMatrix(8,0)+=(crLeftHandSideMatrix48*rB(0,0) + crLeftHandSideMatrix49*rB(1,0) + crLeftHandSideMatrix50*rB(2,0) + crLeftHandSideMatrix51*rB(3,0) + crLeftHandSideMatrix52*rB(4,0) + crLeftHandSideMatrix53*rB(5,0));
// // rLeftHandSideMatrix(8,1)+=(crLeftHandSideMatrix48*rB(0,1) + crLeftHandSideMatrix49*rB(1,1) + crLeftHandSideMatrix50*rB(2,1) + crLeftHandSideMatrix51*rB(3,1) + crLeftHandSideMatrix52*rB(4,1) + crLeftHandSideMatrix53*rB(5,1));
// // rLeftHandSideMatrix(8,2)+=(crLeftHandSideMatrix48*rB(0,2) + crLeftHandSideMatrix49*rB(1,2) + crLeftHandSideMatrix50*rB(2,2) + crLeftHandSideMatrix51*rB(3,2) + crLeftHandSideMatrix52*rB(4,2) + crLeftHandSideMatrix53*rB(5,2));
// // rLeftHandSideMatrix(8,3)+=(crLeftHandSideMatrix48*rB(0,3) + crLeftHandSideMatrix49*rB(1,3) + crLeftHandSideMatrix50*rB(2,3) + crLeftHandSideMatrix51*rB(3,3) + crLeftHandSideMatrix52*rB(4,3) + crLeftHandSideMatrix53*rB(5,3));
// // rLeftHandSideMatrix(8,4)+=(crLeftHandSideMatrix48*rB(0,4) + crLeftHandSideMatrix49*rB(1,4) + crLeftHandSideMatrix50*rB(2,4) + crLeftHandSideMatrix51*rB(3,4) + crLeftHandSideMatrix52*rB(4,4) + crLeftHandSideMatrix53*rB(5,4));
// // rLeftHandSideMatrix(8,5)+=(crLeftHandSideMatrix48*rB(0,5) + crLeftHandSideMatrix49*rB(1,5) + crLeftHandSideMatrix50*rB(2,5) + crLeftHandSideMatrix51*rB(3,5) + crLeftHandSideMatrix52*rB(4,5) + crLeftHandSideMatrix53*rB(5,5));
// // rLeftHandSideMatrix(8,6)+=(crLeftHandSideMatrix48*rB(0,6) + crLeftHandSideMatrix49*rB(1,6) + crLeftHandSideMatrix50*rB(2,6) + crLeftHandSideMatrix51*rB(3,6) + crLeftHandSideMatrix52*rB(4,6) + crLeftHandSideMatrix53*rB(5,6));
// // rLeftHandSideMatrix(8,7)+=(crLeftHandSideMatrix48*rB(0,7) + crLeftHandSideMatrix49*rB(1,7) + crLeftHandSideMatrix50*rB(2,7) + crLeftHandSideMatrix51*rB(3,7) + crLeftHandSideMatrix52*rB(4,7) + crLeftHandSideMatrix53*rB(5,7));
// rLeftHandSideMatrix(8,8)+=(crLeftHandSideMatrix48*rB(0,8) + crLeftHandSideMatrix49*rB(1,8) + crLeftHandSideMatrix50*rB(2,8) + crLeftHandSideMatrix51*rB(3,8) + crLeftHandSideMatrix52*rB(4,8) + crLeftHandSideMatrix53*rB(5,8));
// rLeftHandSideMatrix(8,9)+=(crLeftHandSideMatrix48*rB(0,9) + crLeftHandSideMatrix49*rB(1,9) + crLeftHandSideMatrix50*rB(2,9) + crLeftHandSideMatrix51*rB(3,9) + crLeftHandSideMatrix52*rB(4,9) + crLeftHandSideMatrix53*rB(5,9));
// rLeftHandSideMatrix(8,10)+=(crLeftHandSideMatrix48*rB(0,10) + crLeftHandSideMatrix49*rB(1,10) + crLeftHandSideMatrix50*rB(2,10) + crLeftHandSideMatrix51*rB(3,10) + crLeftHandSideMatrix52*rB(4,10) + crLeftHandSideMatrix53*rB(5,10));
// rLeftHandSideMatrix(8,11)+=(crLeftHandSideMatrix48*rB(0,11) + crLeftHandSideMatrix49*rB(1,11) + crLeftHandSideMatrix50*rB(2,11) + crLeftHandSideMatrix51*rB(3,11) + crLeftHandSideMatrix52*rB(4,11) + crLeftHandSideMatrix53*rB(5,11));
// rLeftHandSideMatrix(8,12)+=(crLeftHandSideMatrix48*rB(0,12) + crLeftHandSideMatrix49*rB(1,12) + crLeftHandSideMatrix50*rB(2,12) + crLeftHandSideMatrix51*rB(3,12) + crLeftHandSideMatrix52*rB(4,12) + crLeftHandSideMatrix53*rB(5,12));
// rLeftHandSideMatrix(8,13)+=(crLeftHandSideMatrix48*rB(0,13) + crLeftHandSideMatrix49*rB(1,13) + crLeftHandSideMatrix50*rB(2,13) + crLeftHandSideMatrix51*rB(3,13) + crLeftHandSideMatrix52*rB(4,13) + crLeftHandSideMatrix53*rB(5,13));
// rLeftHandSideMatrix(8,14)+=(crLeftHandSideMatrix48*rB(0,14) + crLeftHandSideMatrix49*rB(1,14) + crLeftHandSideMatrix50*rB(2,14) + crLeftHandSideMatrix51*rB(3,14) + crLeftHandSideMatrix52*rB(4,14) + crLeftHandSideMatrix53*rB(5,14));
// rLeftHandSideMatrix(8,15)+=(crLeftHandSideMatrix48*rB(0,15) + crLeftHandSideMatrix49*rB(1,15) + crLeftHandSideMatrix50*rB(2,15) + crLeftHandSideMatrix51*rB(3,15) + crLeftHandSideMatrix52*rB(4,15) + crLeftHandSideMatrix53*rB(5,15));
// rLeftHandSideMatrix(8,16)+=(crLeftHandSideMatrix48*rB(0,16) + crLeftHandSideMatrix49*rB(1,16) + crLeftHandSideMatrix50*rB(2,16) + crLeftHandSideMatrix51*rB(3,16) + crLeftHandSideMatrix52*rB(4,16) + crLeftHandSideMatrix53*rB(5,16));
// rLeftHandSideMatrix(8,17)+=(crLeftHandSideMatrix48*rB(0,17) + crLeftHandSideMatrix49*rB(1,17) + crLeftHandSideMatrix50*rB(2,17) + crLeftHandSideMatrix51*rB(3,17) + crLeftHandSideMatrix52*rB(4,17) + crLeftHandSideMatrix53*rB(5,17));
// rLeftHandSideMatrix(8,18)+=(crLeftHandSideMatrix48*rB(0,18) + crLeftHandSideMatrix49*rB(1,18) + crLeftHandSideMatrix50*rB(2,18) + crLeftHandSideMatrix51*rB(3,18) + crLeftHandSideMatrix52*rB(4,18) + crLeftHandSideMatrix53*rB(5,18));
// rLeftHandSideMatrix(8,19)+=(crLeftHandSideMatrix48*rB(0,19) + crLeftHandSideMatrix49*rB(1,19) + crLeftHandSideMatrix50*rB(2,19) + crLeftHandSideMatrix51*rB(3,19) + crLeftHandSideMatrix52*rB(4,19) + crLeftHandSideMatrix53*rB(5,19));
// rLeftHandSideMatrix(8,20)+=(crLeftHandSideMatrix48*rB(0,20) + crLeftHandSideMatrix49*rB(1,20) + crLeftHandSideMatrix50*rB(2,20) + crLeftHandSideMatrix51*rB(3,20) + crLeftHandSideMatrix52*rB(4,20) + crLeftHandSideMatrix53*rB(5,20));
// rLeftHandSideMatrix(8,21)+=(crLeftHandSideMatrix48*rB(0,21) + crLeftHandSideMatrix49*rB(1,21) + crLeftHandSideMatrix50*rB(2,21) + crLeftHandSideMatrix51*rB(3,21) + crLeftHandSideMatrix52*rB(4,21) + crLeftHandSideMatrix53*rB(5,21));
// rLeftHandSideMatrix(8,22)+=(crLeftHandSideMatrix48*rB(0,22) + crLeftHandSideMatrix49*rB(1,22) + crLeftHandSideMatrix50*rB(2,22) + crLeftHandSideMatrix51*rB(3,22) + crLeftHandSideMatrix52*rB(4,22) + crLeftHandSideMatrix53*rB(5,22));
// rLeftHandSideMatrix(8,23)+=(crLeftHandSideMatrix48*rB(0,23) + crLeftHandSideMatrix49*rB(1,23) + crLeftHandSideMatrix50*rB(2,23) + crLeftHandSideMatrix51*rB(3,23) + crLeftHandSideMatrix52*rB(4,23) + crLeftHandSideMatrix53*rB(5,23));

// // rLeftHandSideMatrix(9,0)+=(crLeftHandSideMatrix54*rB(0,0) + crLeftHandSideMatrix55*rB(1,0) + crLeftHandSideMatrix56*rB(2,0) + crLeftHandSideMatrix57*rB(3,0) + crLeftHandSideMatrix58*rB(4,0) + crLeftHandSideMatrix59*rB(5,0));
// // rLeftHandSideMatrix(9,1)+=(crLeftHandSideMatrix54*rB(0,1) + crLeftHandSideMatrix55*rB(1,1) + crLeftHandSideMatrix56*rB(2,1) + crLeftHandSideMatrix57*rB(3,1) + crLeftHandSideMatrix58*rB(4,1) + crLeftHandSideMatrix59*rB(5,1));
// // rLeftHandSideMatrix(9,2)+=(crLeftHandSideMatrix54*rB(0,2) + crLeftHandSideMatrix55*rB(1,2) + crLeftHandSideMatrix56*rB(2,2) + crLeftHandSideMatrix57*rB(3,2) + crLeftHandSideMatrix58*rB(4,2) + crLeftHandSideMatrix59*rB(5,2));
// // rLeftHandSideMatrix(9,3)+=(crLeftHandSideMatrix54*rB(0,3) + crLeftHandSideMatrix55*rB(1,3) + crLeftHandSideMatrix56*rB(2,3) + crLeftHandSideMatrix57*rB(3,3) + crLeftHandSideMatrix58*rB(4,3) + crLeftHandSideMatrix59*rB(5,3));
// // rLeftHandSideMatrix(9,4)+=(crLeftHandSideMatrix54*rB(0,4) + crLeftHandSideMatrix55*rB(1,4) + crLeftHandSideMatrix56*rB(2,4) + crLeftHandSideMatrix57*rB(3,4) + crLeftHandSideMatrix58*rB(4,4) + crLeftHandSideMatrix59*rB(5,4));
// // rLeftHandSideMatrix(9,5)+=(crLeftHandSideMatrix54*rB(0,5) + crLeftHandSideMatrix55*rB(1,5) + crLeftHandSideMatrix56*rB(2,5) + crLeftHandSideMatrix57*rB(3,5) + crLeftHandSideMatrix58*rB(4,5) + crLeftHandSideMatrix59*rB(5,5));
// // rLeftHandSideMatrix(9,6)+=(crLeftHandSideMatrix54*rB(0,6) + crLeftHandSideMatrix55*rB(1,6) + crLeftHandSideMatrix56*rB(2,6) + crLeftHandSideMatrix57*rB(3,6) + crLeftHandSideMatrix58*rB(4,6) + crLeftHandSideMatrix59*rB(5,6));
// // rLeftHandSideMatrix(9,7)+=(crLeftHandSideMatrix54*rB(0,7) + crLeftHandSideMatrix55*rB(1,7) + crLeftHandSideMatrix56*rB(2,7) + crLeftHandSideMatrix57*rB(3,7) + crLeftHandSideMatrix58*rB(4,7) + crLeftHandSideMatrix59*rB(5,7));
// // rLeftHandSideMatrix(9,8)+=(crLeftHandSideMatrix54*rB(0,8) + crLeftHandSideMatrix55*rB(1,8) + crLeftHandSideMatrix56*rB(2,8) + crLeftHandSideMatrix57*rB(3,8) + crLeftHandSideMatrix58*rB(4,8) + crLeftHandSideMatrix59*rB(5,8));
// rLeftHandSideMatrix(9,9)+=(crLeftHandSideMatrix54*rB(0,9) + crLeftHandSideMatrix55*rB(1,9) + crLeftHandSideMatrix56*rB(2,9) + crLeftHandSideMatrix57*rB(3,9) + crLeftHandSideMatrix58*rB(4,9) + crLeftHandSideMatrix59*rB(5,9));
// rLeftHandSideMatrix(9,10)+=(crLeftHandSideMatrix54*rB(0,10) + crLeftHandSideMatrix55*rB(1,10) + crLeftHandSideMatrix56*rB(2,10) + crLeftHandSideMatrix57*rB(3,10) + crLeftHandSideMatrix58*rB(4,10) + crLeftHandSideMatrix59*rB(5,10));
// rLeftHandSideMatrix(9,11)+=(crLeftHandSideMatrix54*rB(0,11) + crLeftHandSideMatrix55*rB(1,11) + crLeftHandSideMatrix56*rB(2,11) + crLeftHandSideMatrix57*rB(3,11) + crLeftHandSideMatrix58*rB(4,11) + crLeftHandSideMatrix59*rB(5,11));
// rLeftHandSideMatrix(9,12)+=(crLeftHandSideMatrix54*rB(0,12) + crLeftHandSideMatrix55*rB(1,12) + crLeftHandSideMatrix56*rB(2,12) + crLeftHandSideMatrix57*rB(3,12) + crLeftHandSideMatrix58*rB(4,12) + crLeftHandSideMatrix59*rB(5,12));
// rLeftHandSideMatrix(9,13)+=(crLeftHandSideMatrix54*rB(0,13) + crLeftHandSideMatrix55*rB(1,13) + crLeftHandSideMatrix56*rB(2,13) + crLeftHandSideMatrix57*rB(3,13) + crLeftHandSideMatrix58*rB(4,13) + crLeftHandSideMatrix59*rB(5,13));
// rLeftHandSideMatrix(9,14)+=(crLeftHandSideMatrix54*rB(0,14) + crLeftHandSideMatrix55*rB(1,14) + crLeftHandSideMatrix56*rB(2,14) + crLeftHandSideMatrix57*rB(3,14) + crLeftHandSideMatrix58*rB(4,14) + crLeftHandSideMatrix59*rB(5,14));
// rLeftHandSideMatrix(9,15)+=(crLeftHandSideMatrix54*rB(0,15) + crLeftHandSideMatrix55*rB(1,15) + crLeftHandSideMatrix56*rB(2,15) + crLeftHandSideMatrix57*rB(3,15) + crLeftHandSideMatrix58*rB(4,15) + crLeftHandSideMatrix59*rB(5,15));
// rLeftHandSideMatrix(9,16)+=(crLeftHandSideMatrix54*rB(0,16) + crLeftHandSideMatrix55*rB(1,16) + crLeftHandSideMatrix56*rB(2,16) + crLeftHandSideMatrix57*rB(3,16) + crLeftHandSideMatrix58*rB(4,16) + crLeftHandSideMatrix59*rB(5,16));
// rLeftHandSideMatrix(9,17)+=(crLeftHandSideMatrix54*rB(0,17) + crLeftHandSideMatrix55*rB(1,17) + crLeftHandSideMatrix56*rB(2,17) + crLeftHandSideMatrix57*rB(3,17) + crLeftHandSideMatrix58*rB(4,17) + crLeftHandSideMatrix59*rB(5,17));
// rLeftHandSideMatrix(9,18)+=(crLeftHandSideMatrix54*rB(0,18) + crLeftHandSideMatrix55*rB(1,18) + crLeftHandSideMatrix56*rB(2,18) + crLeftHandSideMatrix57*rB(3,18) + crLeftHandSideMatrix58*rB(4,18) + crLeftHandSideMatrix59*rB(5,18));
// rLeftHandSideMatrix(9,19)+=(crLeftHandSideMatrix54*rB(0,19) + crLeftHandSideMatrix55*rB(1,19) + crLeftHandSideMatrix56*rB(2,19) + crLeftHandSideMatrix57*rB(3,19) + crLeftHandSideMatrix58*rB(4,19) + crLeftHandSideMatrix59*rB(5,19));
// rLeftHandSideMatrix(9,20)+=(crLeftHandSideMatrix54*rB(0,20) + crLeftHandSideMatrix55*rB(1,20) + crLeftHandSideMatrix56*rB(2,20) + crLeftHandSideMatrix57*rB(3,20) + crLeftHandSideMatrix58*rB(4,20) + crLeftHandSideMatrix59*rB(5,20));
// rLeftHandSideMatrix(9,21)+=(crLeftHandSideMatrix54*rB(0,21) + crLeftHandSideMatrix55*rB(1,21) + crLeftHandSideMatrix56*rB(2,21) + crLeftHandSideMatrix57*rB(3,21) + crLeftHandSideMatrix58*rB(4,21) + crLeftHandSideMatrix59*rB(5,21));
// rLeftHandSideMatrix(9,22)+=(crLeftHandSideMatrix54*rB(0,22) + crLeftHandSideMatrix55*rB(1,22) + crLeftHandSideMatrix56*rB(2,22) + crLeftHandSideMatrix57*rB(3,22) + crLeftHandSideMatrix58*rB(4,22) + crLeftHandSideMatrix59*rB(5,22));
// rLeftHandSideMatrix(9,23)+=(crLeftHandSideMatrix54*rB(0,23) + crLeftHandSideMatrix55*rB(1,23) + crLeftHandSideMatrix56*rB(2,23) + crLeftHandSideMatrix57*rB(3,23) + crLeftHandSideMatrix58*rB(4,23) + crLeftHandSideMatrix59*rB(5,23));

// // rLeftHandSideMatrix(10,0)+=(crLeftHandSideMatrix60*rB(0,0) + crLeftHandSideMatrix61*rB(1,0) + crLeftHandSideMatrix62*rB(2,0) + crLeftHandSideMatrix63*rB(3,0) + crLeftHandSideMatrix64*rB(4,0) + crLeftHandSideMatrix65*rB(5,0));
// // rLeftHandSideMatrix(10,1)+=(crLeftHandSideMatrix60*rB(0,1) + crLeftHandSideMatrix61*rB(1,1) + crLeftHandSideMatrix62*rB(2,1) + crLeftHandSideMatrix63*rB(3,1) + crLeftHandSideMatrix64*rB(4,1) + crLeftHandSideMatrix65*rB(5,1));
// // rLeftHandSideMatrix(10,2)+=(crLeftHandSideMatrix60*rB(0,2) + crLeftHandSideMatrix61*rB(1,2) + crLeftHandSideMatrix62*rB(2,2) + crLeftHandSideMatrix63*rB(3,2) + crLeftHandSideMatrix64*rB(4,2) + crLeftHandSideMatrix65*rB(5,2));
// // rLeftHandSideMatrix(10,3)+=(crLeftHandSideMatrix60*rB(0,3) + crLeftHandSideMatrix61*rB(1,3) + crLeftHandSideMatrix62*rB(2,3) + crLeftHandSideMatrix63*rB(3,3) + crLeftHandSideMatrix64*rB(4,3) + crLeftHandSideMatrix65*rB(5,3));
// // rLeftHandSideMatrix(10,4)+=(crLeftHandSideMatrix60*rB(0,4) + crLeftHandSideMatrix61*rB(1,4) + crLeftHandSideMatrix62*rB(2,4) + crLeftHandSideMatrix63*rB(3,4) + crLeftHandSideMatrix64*rB(4,4) + crLeftHandSideMatrix65*rB(5,4));
// // rLeftHandSideMatrix(10,5)+=(crLeftHandSideMatrix60*rB(0,5) + crLeftHandSideMatrix61*rB(1,5) + crLeftHandSideMatrix62*rB(2,5) + crLeftHandSideMatrix63*rB(3,5) + crLeftHandSideMatrix64*rB(4,5) + crLeftHandSideMatrix65*rB(5,5));
// // rLeftHandSideMatrix(10,6)+=(crLeftHandSideMatrix60*rB(0,6) + crLeftHandSideMatrix61*rB(1,6) + crLeftHandSideMatrix62*rB(2,6) + crLeftHandSideMatrix63*rB(3,6) + crLeftHandSideMatrix64*rB(4,6) + crLeftHandSideMatrix65*rB(5,6));
// // rLeftHandSideMatrix(10,7)+=(crLeftHandSideMatrix60*rB(0,7) + crLeftHandSideMatrix61*rB(1,7) + crLeftHandSideMatrix62*rB(2,7) + crLeftHandSideMatrix63*rB(3,7) + crLeftHandSideMatrix64*rB(4,7) + crLeftHandSideMatrix65*rB(5,7));
// // rLeftHandSideMatrix(10,8)+=(crLeftHandSideMatrix60*rB(0,8) + crLeftHandSideMatrix61*rB(1,8) + crLeftHandSideMatrix62*rB(2,8) + crLeftHandSideMatrix63*rB(3,8) + crLeftHandSideMatrix64*rB(4,8) + crLeftHandSideMatrix65*rB(5,8));
// // rLeftHandSideMatrix(10,9)+=(crLeftHandSideMatrix60*rB(0,9) + crLeftHandSideMatrix61*rB(1,9) + crLeftHandSideMatrix62*rB(2,9) + crLeftHandSideMatrix63*rB(3,9) + crLeftHandSideMatrix64*rB(4,9) + crLeftHandSideMatrix65*rB(5,9));
// rLeftHandSideMatrix(10,10)+=(crLeftHandSideMatrix60*rB(0,10) + crLeftHandSideMatrix61*rB(1,10) + crLeftHandSideMatrix62*rB(2,10) + crLeftHandSideMatrix63*rB(3,10) + crLeftHandSideMatrix64*rB(4,10) + crLeftHandSideMatrix65*rB(5,10));
// rLeftHandSideMatrix(10,11)+=(crLeftHandSideMatrix60*rB(0,11) + crLeftHandSideMatrix61*rB(1,11) + crLeftHandSideMatrix62*rB(2,11) + crLeftHandSideMatrix63*rB(3,11) + crLeftHandSideMatrix64*rB(4,11) + crLeftHandSideMatrix65*rB(5,11));
// rLeftHandSideMatrix(10,12)+=(crLeftHandSideMatrix60*rB(0,12) + crLeftHandSideMatrix61*rB(1,12) + crLeftHandSideMatrix62*rB(2,12) + crLeftHandSideMatrix63*rB(3,12) + crLeftHandSideMatrix64*rB(4,12) + crLeftHandSideMatrix65*rB(5,12));
// rLeftHandSideMatrix(10,13)+=(crLeftHandSideMatrix60*rB(0,13) + crLeftHandSideMatrix61*rB(1,13) + crLeftHandSideMatrix62*rB(2,13) + crLeftHandSideMatrix63*rB(3,13) + crLeftHandSideMatrix64*rB(4,13) + crLeftHandSideMatrix65*rB(5,13));
// rLeftHandSideMatrix(10,14)+=(crLeftHandSideMatrix60*rB(0,14) + crLeftHandSideMatrix61*rB(1,14) + crLeftHandSideMatrix62*rB(2,14) + crLeftHandSideMatrix63*rB(3,14) + crLeftHandSideMatrix64*rB(4,14) + crLeftHandSideMatrix65*rB(5,14));
// rLeftHandSideMatrix(10,15)+=(crLeftHandSideMatrix60*rB(0,15) + crLeftHandSideMatrix61*rB(1,15) + crLeftHandSideMatrix62*rB(2,15) + crLeftHandSideMatrix63*rB(3,15) + crLeftHandSideMatrix64*rB(4,15) + crLeftHandSideMatrix65*rB(5,15));
// rLeftHandSideMatrix(10,16)+=(crLeftHandSideMatrix60*rB(0,16) + crLeftHandSideMatrix61*rB(1,16) + crLeftHandSideMatrix62*rB(2,16) + crLeftHandSideMatrix63*rB(3,16) + crLeftHandSideMatrix64*rB(4,16) + crLeftHandSideMatrix65*rB(5,16));
// rLeftHandSideMatrix(10,17)+=(crLeftHandSideMatrix60*rB(0,17) + crLeftHandSideMatrix61*rB(1,17) + crLeftHandSideMatrix62*rB(2,17) + crLeftHandSideMatrix63*rB(3,17) + crLeftHandSideMatrix64*rB(4,17) + crLeftHandSideMatrix65*rB(5,17));
// rLeftHandSideMatrix(10,18)+=(crLeftHandSideMatrix60*rB(0,18) + crLeftHandSideMatrix61*rB(1,18) + crLeftHandSideMatrix62*rB(2,18) + crLeftHandSideMatrix63*rB(3,18) + crLeftHandSideMatrix64*rB(4,18) + crLeftHandSideMatrix65*rB(5,18));
// rLeftHandSideMatrix(10,19)+=(crLeftHandSideMatrix60*rB(0,19) + crLeftHandSideMatrix61*rB(1,19) + crLeftHandSideMatrix62*rB(2,19) + crLeftHandSideMatrix63*rB(3,19) + crLeftHandSideMatrix64*rB(4,19) + crLeftHandSideMatrix65*rB(5,19));
// rLeftHandSideMatrix(10,20)+=(crLeftHandSideMatrix60*rB(0,20) + crLeftHandSideMatrix61*rB(1,20) + crLeftHandSideMatrix62*rB(2,20) + crLeftHandSideMatrix63*rB(3,20) + crLeftHandSideMatrix64*rB(4,20) + crLeftHandSideMatrix65*rB(5,20));
// rLeftHandSideMatrix(10,21)+=(crLeftHandSideMatrix60*rB(0,21) + crLeftHandSideMatrix61*rB(1,21) + crLeftHandSideMatrix62*rB(2,21) + crLeftHandSideMatrix63*rB(3,21) + crLeftHandSideMatrix64*rB(4,21) + crLeftHandSideMatrix65*rB(5,21));
// rLeftHandSideMatrix(10,22)+=(crLeftHandSideMatrix60*rB(0,22) + crLeftHandSideMatrix61*rB(1,22) + crLeftHandSideMatrix62*rB(2,22) + crLeftHandSideMatrix63*rB(3,22) + crLeftHandSideMatrix64*rB(4,22) + crLeftHandSideMatrix65*rB(5,22));
// rLeftHandSideMatrix(10,23)+=(crLeftHandSideMatrix60*rB(0,23) + crLeftHandSideMatrix61*rB(1,23) + crLeftHandSideMatrix62*rB(2,23) + crLeftHandSideMatrix63*rB(3,23) + crLeftHandSideMatrix64*rB(4,23) + crLeftHandSideMatrix65*rB(5,23));

// // rLeftHandSideMatrix(11,0)+=(crLeftHandSideMatrix66*rB(0,0) + crLeftHandSideMatrix67*rB(1,0) + crLeftHandSideMatrix68*rB(2,0) + crLeftHandSideMatrix69*rB(3,0) + crLeftHandSideMatrix70*rB(4,0) + crLeftHandSideMatrix71*rB(5,0));
// // rLeftHandSideMatrix(11,1)+=(crLeftHandSideMatrix66*rB(0,1) + crLeftHandSideMatrix67*rB(1,1) + crLeftHandSideMatrix68*rB(2,1) + crLeftHandSideMatrix69*rB(3,1) + crLeftHandSideMatrix70*rB(4,1) + crLeftHandSideMatrix71*rB(5,1));
// // rLeftHandSideMatrix(11,2)+=(crLeftHandSideMatrix66*rB(0,2) + crLeftHandSideMatrix67*rB(1,2) + crLeftHandSideMatrix68*rB(2,2) + crLeftHandSideMatrix69*rB(3,2) + crLeftHandSideMatrix70*rB(4,2) + crLeftHandSideMatrix71*rB(5,2));
// // rLeftHandSideMatrix(11,3)+=(crLeftHandSideMatrix66*rB(0,3) + crLeftHandSideMatrix67*rB(1,3) + crLeftHandSideMatrix68*rB(2,3) + crLeftHandSideMatrix69*rB(3,3) + crLeftHandSideMatrix70*rB(4,3) + crLeftHandSideMatrix71*rB(5,3));
// // rLeftHandSideMatrix(11,4)+=(crLeftHandSideMatrix66*rB(0,4) + crLeftHandSideMatrix67*rB(1,4) + crLeftHandSideMatrix68*rB(2,4) + crLeftHandSideMatrix69*rB(3,4) + crLeftHandSideMatrix70*rB(4,4) + crLeftHandSideMatrix71*rB(5,4));
// // rLeftHandSideMatrix(11,5)+=(crLeftHandSideMatrix66*rB(0,5) + crLeftHandSideMatrix67*rB(1,5) + crLeftHandSideMatrix68*rB(2,5) + crLeftHandSideMatrix69*rB(3,5) + crLeftHandSideMatrix70*rB(4,5) + crLeftHandSideMatrix71*rB(5,5));
// // rLeftHandSideMatrix(11,6)+=(crLeftHandSideMatrix66*rB(0,6) + crLeftHandSideMatrix67*rB(1,6) + crLeftHandSideMatrix68*rB(2,6) + crLeftHandSideMatrix69*rB(3,6) + crLeftHandSideMatrix70*rB(4,6) + crLeftHandSideMatrix71*rB(5,6));
// // rLeftHandSideMatrix(11,7)+=(crLeftHandSideMatrix66*rB(0,7) + crLeftHandSideMatrix67*rB(1,7) + crLeftHandSideMatrix68*rB(2,7) + crLeftHandSideMatrix69*rB(3,7) + crLeftHandSideMatrix70*rB(4,7) + crLeftHandSideMatrix71*rB(5,7));
// // rLeftHandSideMatrix(11,8)+=(crLeftHandSideMatrix66*rB(0,8) + crLeftHandSideMatrix67*rB(1,8) + crLeftHandSideMatrix68*rB(2,8) + crLeftHandSideMatrix69*rB(3,8) + crLeftHandSideMatrix70*rB(4,8) + crLeftHandSideMatrix71*rB(5,8));
// // rLeftHandSideMatrix(11,9)+=(crLeftHandSideMatrix66*rB(0,9) + crLeftHandSideMatrix67*rB(1,9) + crLeftHandSideMatrix68*rB(2,9) + crLeftHandSideMatrix69*rB(3,9) + crLeftHandSideMatrix70*rB(4,9) + crLeftHandSideMatrix71*rB(5,9));
// // rLeftHandSideMatrix(11,10)+=(crLeftHandSideMatrix66*rB(0,10) + crLeftHandSideMatrix67*rB(1,10) + crLeftHandSideMatrix68*rB(2,10) + crLeftHandSideMatrix69*rB(3,10) + crLeftHandSideMatrix70*rB(4,10) + crLeftHandSideMatrix71*rB(5,10));
// rLeftHandSideMatrix(11,11)+=(crLeftHandSideMatrix66*rB(0,11) + crLeftHandSideMatrix67*rB(1,11) + crLeftHandSideMatrix68*rB(2,11) + crLeftHandSideMatrix69*rB(3,11) + crLeftHandSideMatrix70*rB(4,11) + crLeftHandSideMatrix71*rB(5,11));
// rLeftHandSideMatrix(11,12)+=(crLeftHandSideMatrix66*rB(0,12) + crLeftHandSideMatrix67*rB(1,12) + crLeftHandSideMatrix68*rB(2,12) + crLeftHandSideMatrix69*rB(3,12) + crLeftHandSideMatrix70*rB(4,12) + crLeftHandSideMatrix71*rB(5,12));
// rLeftHandSideMatrix(11,13)+=(crLeftHandSideMatrix66*rB(0,13) + crLeftHandSideMatrix67*rB(1,13) + crLeftHandSideMatrix68*rB(2,13) + crLeftHandSideMatrix69*rB(3,13) + crLeftHandSideMatrix70*rB(4,13) + crLeftHandSideMatrix71*rB(5,13));
// rLeftHandSideMatrix(11,14)+=(crLeftHandSideMatrix66*rB(0,14) + crLeftHandSideMatrix67*rB(1,14) + crLeftHandSideMatrix68*rB(2,14) + crLeftHandSideMatrix69*rB(3,14) + crLeftHandSideMatrix70*rB(4,14) + crLeftHandSideMatrix71*rB(5,14));
// rLeftHandSideMatrix(11,15)+=(crLeftHandSideMatrix66*rB(0,15) + crLeftHandSideMatrix67*rB(1,15) + crLeftHandSideMatrix68*rB(2,15) + crLeftHandSideMatrix69*rB(3,15) + crLeftHandSideMatrix70*rB(4,15) + crLeftHandSideMatrix71*rB(5,15));
// rLeftHandSideMatrix(11,16)+=(crLeftHandSideMatrix66*rB(0,16) + crLeftHandSideMatrix67*rB(1,16) + crLeftHandSideMatrix68*rB(2,16) + crLeftHandSideMatrix69*rB(3,16) + crLeftHandSideMatrix70*rB(4,16) + crLeftHandSideMatrix71*rB(5,16));
// rLeftHandSideMatrix(11,17)+=(crLeftHandSideMatrix66*rB(0,17) + crLeftHandSideMatrix67*rB(1,17) + crLeftHandSideMatrix68*rB(2,17) + crLeftHandSideMatrix69*rB(3,17) + crLeftHandSideMatrix70*rB(4,17) + crLeftHandSideMatrix71*rB(5,17));
// rLeftHandSideMatrix(11,18)+=(crLeftHandSideMatrix66*rB(0,18) + crLeftHandSideMatrix67*rB(1,18) + crLeftHandSideMatrix68*rB(2,18) + crLeftHandSideMatrix69*rB(3,18) + crLeftHandSideMatrix70*rB(4,18) + crLeftHandSideMatrix71*rB(5,18));
// rLeftHandSideMatrix(11,19)+=(crLeftHandSideMatrix66*rB(0,19) + crLeftHandSideMatrix67*rB(1,19) + crLeftHandSideMatrix68*rB(2,19) + crLeftHandSideMatrix69*rB(3,19) + crLeftHandSideMatrix70*rB(4,19) + crLeftHandSideMatrix71*rB(5,19));
// rLeftHandSideMatrix(11,20)+=(crLeftHandSideMatrix66*rB(0,20) + crLeftHandSideMatrix67*rB(1,20) + crLeftHandSideMatrix68*rB(2,20) + crLeftHandSideMatrix69*rB(3,20) + crLeftHandSideMatrix70*rB(4,20) + crLeftHandSideMatrix71*rB(5,20));
// rLeftHandSideMatrix(11,21)+=(crLeftHandSideMatrix66*rB(0,21) + crLeftHandSideMatrix67*rB(1,21) + crLeftHandSideMatrix68*rB(2,21) + crLeftHandSideMatrix69*rB(3,21) + crLeftHandSideMatrix70*rB(4,21) + crLeftHandSideMatrix71*rB(5,21));
// rLeftHandSideMatrix(11,22)+=(crLeftHandSideMatrix66*rB(0,22) + crLeftHandSideMatrix67*rB(1,22) + crLeftHandSideMatrix68*rB(2,22) + crLeftHandSideMatrix69*rB(3,22) + crLeftHandSideMatrix70*rB(4,22) + crLeftHandSideMatrix71*rB(5,22));
// rLeftHandSideMatrix(11,23)+=(crLeftHandSideMatrix66*rB(0,23) + crLeftHandSideMatrix67*rB(1,23) + crLeftHandSideMatrix68*rB(2,23) + crLeftHandSideMatrix69*rB(3,23) + crLeftHandSideMatrix70*rB(4,23) + crLeftHandSideMatrix71*rB(5,23));

// // rLeftHandSideMatrix(12,0)+=(crLeftHandSideMatrix72*rB(0,0) + crLeftHandSideMatrix73*rB(1,0) + crLeftHandSideMatrix74*rB(2,0) + crLeftHandSideMatrix75*rB(3,0) + crLeftHandSideMatrix76*rB(4,0) + crLeftHandSideMatrix77*rB(5,0));
// // rLeftHandSideMatrix(12,1)+=(crLeftHandSideMatrix72*rB(0,1) + crLeftHandSideMatrix73*rB(1,1) + crLeftHandSideMatrix74*rB(2,1) + crLeftHandSideMatrix75*rB(3,1) + crLeftHandSideMatrix76*rB(4,1) + crLeftHandSideMatrix77*rB(5,1));
// // rLeftHandSideMatrix(12,2)+=(crLeftHandSideMatrix72*rB(0,2) + crLeftHandSideMatrix73*rB(1,2) + crLeftHandSideMatrix74*rB(2,2) + crLeftHandSideMatrix75*rB(3,2) + crLeftHandSideMatrix76*rB(4,2) + crLeftHandSideMatrix77*rB(5,2));
// // rLeftHandSideMatrix(12,3)+=(crLeftHandSideMatrix72*rB(0,3) + crLeftHandSideMatrix73*rB(1,3) + crLeftHandSideMatrix74*rB(2,3) + crLeftHandSideMatrix75*rB(3,3) + crLeftHandSideMatrix76*rB(4,3) + crLeftHandSideMatrix77*rB(5,3));
// // rLeftHandSideMatrix(12,4)+=(crLeftHandSideMatrix72*rB(0,4) + crLeftHandSideMatrix73*rB(1,4) + crLeftHandSideMatrix74*rB(2,4) + crLeftHandSideMatrix75*rB(3,4) + crLeftHandSideMatrix76*rB(4,4) + crLeftHandSideMatrix77*rB(5,4));
// // rLeftHandSideMatrix(12,5)+=(crLeftHandSideMatrix72*rB(0,5) + crLeftHandSideMatrix73*rB(1,5) + crLeftHandSideMatrix74*rB(2,5) + crLeftHandSideMatrix75*rB(3,5) + crLeftHandSideMatrix76*rB(4,5) + crLeftHandSideMatrix77*rB(5,5));
// // rLeftHandSideMatrix(12,6)+=(crLeftHandSideMatrix72*rB(0,6) + crLeftHandSideMatrix73*rB(1,6) + crLeftHandSideMatrix74*rB(2,6) + crLeftHandSideMatrix75*rB(3,6) + crLeftHandSideMatrix76*rB(4,6) + crLeftHandSideMatrix77*rB(5,6));
// // rLeftHandSideMatrix(12,7)+=(crLeftHandSideMatrix72*rB(0,7) + crLeftHandSideMatrix73*rB(1,7) + crLeftHandSideMatrix74*rB(2,7) + crLeftHandSideMatrix75*rB(3,7) + crLeftHandSideMatrix76*rB(4,7) + crLeftHandSideMatrix77*rB(5,7));
// // rLeftHandSideMatrix(12,8)+=(crLeftHandSideMatrix72*rB(0,8) + crLeftHandSideMatrix73*rB(1,8) + crLeftHandSideMatrix74*rB(2,8) + crLeftHandSideMatrix75*rB(3,8) + crLeftHandSideMatrix76*rB(4,8) + crLeftHandSideMatrix77*rB(5,8));
// // rLeftHandSideMatrix(12,9)+=(crLeftHandSideMatrix72*rB(0,9) + crLeftHandSideMatrix73*rB(1,9) + crLeftHandSideMatrix74*rB(2,9) + crLeftHandSideMatrix75*rB(3,9) + crLeftHandSideMatrix76*rB(4,9) + crLeftHandSideMatrix77*rB(5,9));
// // rLeftHandSideMatrix(12,10)+=(crLeftHandSideMatrix72*rB(0,10) + crLeftHandSideMatrix73*rB(1,10) + crLeftHandSideMatrix74*rB(2,10) + crLeftHandSideMatrix75*rB(3,10) + crLeftHandSideMatrix76*rB(4,10) + crLeftHandSideMatrix77*rB(5,10));
// // rLeftHandSideMatrix(12,11)+=(crLeftHandSideMatrix72*rB(0,11) + crLeftHandSideMatrix73*rB(1,11) + crLeftHandSideMatrix74*rB(2,11) + crLeftHandSideMatrix75*rB(3,11) + crLeftHandSideMatrix76*rB(4,11) + crLeftHandSideMatrix77*rB(5,11));
// rLeftHandSideMatrix(12,12)+=(crLeftHandSideMatrix72*rB(0,12) + crLeftHandSideMatrix73*rB(1,12) + crLeftHandSideMatrix74*rB(2,12) + crLeftHandSideMatrix75*rB(3,12) + crLeftHandSideMatrix76*rB(4,12) + crLeftHandSideMatrix77*rB(5,12));
// rLeftHandSideMatrix(12,13)+=(crLeftHandSideMatrix72*rB(0,13) + crLeftHandSideMatrix73*rB(1,13) + crLeftHandSideMatrix74*rB(2,13) + crLeftHandSideMatrix75*rB(3,13) + crLeftHandSideMatrix76*rB(4,13) + crLeftHandSideMatrix77*rB(5,13));
// rLeftHandSideMatrix(12,14)+=(crLeftHandSideMatrix72*rB(0,14) + crLeftHandSideMatrix73*rB(1,14) + crLeftHandSideMatrix74*rB(2,14) + crLeftHandSideMatrix75*rB(3,14) + crLeftHandSideMatrix76*rB(4,14) + crLeftHandSideMatrix77*rB(5,14));
// rLeftHandSideMatrix(12,15)+=(crLeftHandSideMatrix72*rB(0,15) + crLeftHandSideMatrix73*rB(1,15) + crLeftHandSideMatrix74*rB(2,15) + crLeftHandSideMatrix75*rB(3,15) + crLeftHandSideMatrix76*rB(4,15) + crLeftHandSideMatrix77*rB(5,15));
// rLeftHandSideMatrix(12,16)+=(crLeftHandSideMatrix72*rB(0,16) + crLeftHandSideMatrix73*rB(1,16) + crLeftHandSideMatrix74*rB(2,16) + crLeftHandSideMatrix75*rB(3,16) + crLeftHandSideMatrix76*rB(4,16) + crLeftHandSideMatrix77*rB(5,16));
// rLeftHandSideMatrix(12,17)+=(crLeftHandSideMatrix72*rB(0,17) + crLeftHandSideMatrix73*rB(1,17) + crLeftHandSideMatrix74*rB(2,17) + crLeftHandSideMatrix75*rB(3,17) + crLeftHandSideMatrix76*rB(4,17) + crLeftHandSideMatrix77*rB(5,17));
// rLeftHandSideMatrix(12,18)+=(crLeftHandSideMatrix72*rB(0,18) + crLeftHandSideMatrix73*rB(1,18) + crLeftHandSideMatrix74*rB(2,18) + crLeftHandSideMatrix75*rB(3,18) + crLeftHandSideMatrix76*rB(4,18) + crLeftHandSideMatrix77*rB(5,18));
// rLeftHandSideMatrix(12,19)+=(crLeftHandSideMatrix72*rB(0,19) + crLeftHandSideMatrix73*rB(1,19) + crLeftHandSideMatrix74*rB(2,19) + crLeftHandSideMatrix75*rB(3,19) + crLeftHandSideMatrix76*rB(4,19) + crLeftHandSideMatrix77*rB(5,19));
// rLeftHandSideMatrix(12,20)+=(crLeftHandSideMatrix72*rB(0,20) + crLeftHandSideMatrix73*rB(1,20) + crLeftHandSideMatrix74*rB(2,20) + crLeftHandSideMatrix75*rB(3,20) + crLeftHandSideMatrix76*rB(4,20) + crLeftHandSideMatrix77*rB(5,20));
// rLeftHandSideMatrix(12,21)+=(crLeftHandSideMatrix72*rB(0,21) + crLeftHandSideMatrix73*rB(1,21) + crLeftHandSideMatrix74*rB(2,21) + crLeftHandSideMatrix75*rB(3,21) + crLeftHandSideMatrix76*rB(4,21) + crLeftHandSideMatrix77*rB(5,21));
// rLeftHandSideMatrix(12,22)+=(crLeftHandSideMatrix72*rB(0,22) + crLeftHandSideMatrix73*rB(1,22) + crLeftHandSideMatrix74*rB(2,22) + crLeftHandSideMatrix75*rB(3,22) + crLeftHandSideMatrix76*rB(4,22) + crLeftHandSideMatrix77*rB(5,22));
// rLeftHandSideMatrix(12,23)+=(crLeftHandSideMatrix72*rB(0,23) + crLeftHandSideMatrix73*rB(1,23) + crLeftHandSideMatrix74*rB(2,23) + crLeftHandSideMatrix75*rB(3,23) + crLeftHandSideMatrix76*rB(4,23) + crLeftHandSideMatrix77*rB(5,23));

// // rLeftHandSideMatrix(13,0)+=(crLeftHandSideMatrix78*rB(0,0) + crLeftHandSideMatrix79*rB(1,0) + crLeftHandSideMatrix80*rB(2,0) + crLeftHandSideMatrix81*rB(3,0) + crLeftHandSideMatrix82*rB(4,0) + crLeftHandSideMatrix83*rB(5,0));
// // rLeftHandSideMatrix(13,1)+=(crLeftHandSideMatrix78*rB(0,1) + crLeftHandSideMatrix79*rB(1,1) + crLeftHandSideMatrix80*rB(2,1) + crLeftHandSideMatrix81*rB(3,1) + crLeftHandSideMatrix82*rB(4,1) + crLeftHandSideMatrix83*rB(5,1));
// // rLeftHandSideMatrix(13,2)+=(crLeftHandSideMatrix78*rB(0,2) + crLeftHandSideMatrix79*rB(1,2) + crLeftHandSideMatrix80*rB(2,2) + crLeftHandSideMatrix81*rB(3,2) + crLeftHandSideMatrix82*rB(4,2) + crLeftHandSideMatrix83*rB(5,2));
// // rLeftHandSideMatrix(13,3)+=(crLeftHandSideMatrix78*rB(0,3) + crLeftHandSideMatrix79*rB(1,3) + crLeftHandSideMatrix80*rB(2,3) + crLeftHandSideMatrix81*rB(3,3) + crLeftHandSideMatrix82*rB(4,3) + crLeftHandSideMatrix83*rB(5,3));
// // rLeftHandSideMatrix(13,4)+=(crLeftHandSideMatrix78*rB(0,4) + crLeftHandSideMatrix79*rB(1,4) + crLeftHandSideMatrix80*rB(2,4) + crLeftHandSideMatrix81*rB(3,4) + crLeftHandSideMatrix82*rB(4,4) + crLeftHandSideMatrix83*rB(5,4));
// // rLeftHandSideMatrix(13,5)+=(crLeftHandSideMatrix78*rB(0,5) + crLeftHandSideMatrix79*rB(1,5) + crLeftHandSideMatrix80*rB(2,5) + crLeftHandSideMatrix81*rB(3,5) + crLeftHandSideMatrix82*rB(4,5) + crLeftHandSideMatrix83*rB(5,5));
// // rLeftHandSideMatrix(13,6)+=(crLeftHandSideMatrix78*rB(0,6) + crLeftHandSideMatrix79*rB(1,6) + crLeftHandSideMatrix80*rB(2,6) + crLeftHandSideMatrix81*rB(3,6) + crLeftHandSideMatrix82*rB(4,6) + crLeftHandSideMatrix83*rB(5,6));
// // rLeftHandSideMatrix(13,7)+=(crLeftHandSideMatrix78*rB(0,7) + crLeftHandSideMatrix79*rB(1,7) + crLeftHandSideMatrix80*rB(2,7) + crLeftHandSideMatrix81*rB(3,7) + crLeftHandSideMatrix82*rB(4,7) + crLeftHandSideMatrix83*rB(5,7));
// // rLeftHandSideMatrix(13,8)+=(crLeftHandSideMatrix78*rB(0,8) + crLeftHandSideMatrix79*rB(1,8) + crLeftHandSideMatrix80*rB(2,8) + crLeftHandSideMatrix81*rB(3,8) + crLeftHandSideMatrix82*rB(4,8) + crLeftHandSideMatrix83*rB(5,8));
// // rLeftHandSideMatrix(13,9)+=(crLeftHandSideMatrix78*rB(0,9) + crLeftHandSideMatrix79*rB(1,9) + crLeftHandSideMatrix80*rB(2,9) + crLeftHandSideMatrix81*rB(3,9) + crLeftHandSideMatrix82*rB(4,9) + crLeftHandSideMatrix83*rB(5,9));
// // rLeftHandSideMatrix(13,10)+=(crLeftHandSideMatrix78*rB(0,10) + crLeftHandSideMatrix79*rB(1,10) + crLeftHandSideMatrix80*rB(2,10) + crLeftHandSideMatrix81*rB(3,10) + crLeftHandSideMatrix82*rB(4,10) + crLeftHandSideMatrix83*rB(5,10));
// // rLeftHandSideMatrix(13,11)+=(crLeftHandSideMatrix78*rB(0,11) + crLeftHandSideMatrix79*rB(1,11) + crLeftHandSideMatrix80*rB(2,11) + crLeftHandSideMatrix81*rB(3,11) + crLeftHandSideMatrix82*rB(4,11) + crLeftHandSideMatrix83*rB(5,11));
// // rLeftHandSideMatrix(13,12)+=(crLeftHandSideMatrix78*rB(0,12) + crLeftHandSideMatrix79*rB(1,12) + crLeftHandSideMatrix80*rB(2,12) + crLeftHandSideMatrix81*rB(3,12) + crLeftHandSideMatrix82*rB(4,12) + crLeftHandSideMatrix83*rB(5,12));
// rLeftHandSideMatrix(13,13)+=(crLeftHandSideMatrix78*rB(0,13) + crLeftHandSideMatrix79*rB(1,13) + crLeftHandSideMatrix80*rB(2,13) + crLeftHandSideMatrix81*rB(3,13) + crLeftHandSideMatrix82*rB(4,13) + crLeftHandSideMatrix83*rB(5,13));
// rLeftHandSideMatrix(13,14)+=(crLeftHandSideMatrix78*rB(0,14) + crLeftHandSideMatrix79*rB(1,14) + crLeftHandSideMatrix80*rB(2,14) + crLeftHandSideMatrix81*rB(3,14) + crLeftHandSideMatrix82*rB(4,14) + crLeftHandSideMatrix83*rB(5,14));
// rLeftHandSideMatrix(13,15)+=(crLeftHandSideMatrix78*rB(0,15) + crLeftHandSideMatrix79*rB(1,15) + crLeftHandSideMatrix80*rB(2,15) + crLeftHandSideMatrix81*rB(3,15) + crLeftHandSideMatrix82*rB(4,15) + crLeftHandSideMatrix83*rB(5,15));
// rLeftHandSideMatrix(13,16)+=(crLeftHandSideMatrix78*rB(0,16) + crLeftHandSideMatrix79*rB(1,16) + crLeftHandSideMatrix80*rB(2,16) + crLeftHandSideMatrix81*rB(3,16) + crLeftHandSideMatrix82*rB(4,16) + crLeftHandSideMatrix83*rB(5,16));
// rLeftHandSideMatrix(13,17)+=(crLeftHandSideMatrix78*rB(0,17) + crLeftHandSideMatrix79*rB(1,17) + crLeftHandSideMatrix80*rB(2,17) + crLeftHandSideMatrix81*rB(3,17) + crLeftHandSideMatrix82*rB(4,17) + crLeftHandSideMatrix83*rB(5,17));
// rLeftHandSideMatrix(13,18)+=(crLeftHandSideMatrix78*rB(0,18) + crLeftHandSideMatrix79*rB(1,18) + crLeftHandSideMatrix80*rB(2,18) + crLeftHandSideMatrix81*rB(3,18) + crLeftHandSideMatrix82*rB(4,18) + crLeftHandSideMatrix83*rB(5,18));
// rLeftHandSideMatrix(13,19)+=(crLeftHandSideMatrix78*rB(0,19) + crLeftHandSideMatrix79*rB(1,19) + crLeftHandSideMatrix80*rB(2,19) + crLeftHandSideMatrix81*rB(3,19) + crLeftHandSideMatrix82*rB(4,19) + crLeftHandSideMatrix83*rB(5,19));
// rLeftHandSideMatrix(13,20)+=(crLeftHandSideMatrix78*rB(0,20) + crLeftHandSideMatrix79*rB(1,20) + crLeftHandSideMatrix80*rB(2,20) + crLeftHandSideMatrix81*rB(3,20) + crLeftHandSideMatrix82*rB(4,20) + crLeftHandSideMatrix83*rB(5,20));
// rLeftHandSideMatrix(13,21)+=(crLeftHandSideMatrix78*rB(0,21) + crLeftHandSideMatrix79*rB(1,21) + crLeftHandSideMatrix80*rB(2,21) + crLeftHandSideMatrix81*rB(3,21) + crLeftHandSideMatrix82*rB(4,21) + crLeftHandSideMatrix83*rB(5,21));
// rLeftHandSideMatrix(13,22)+=(crLeftHandSideMatrix78*rB(0,22) + crLeftHandSideMatrix79*rB(1,22) + crLeftHandSideMatrix80*rB(2,22) + crLeftHandSideMatrix81*rB(3,22) + crLeftHandSideMatrix82*rB(4,22) + crLeftHandSideMatrix83*rB(5,22));
// rLeftHandSideMatrix(13,23)+=(crLeftHandSideMatrix78*rB(0,23) + crLeftHandSideMatrix79*rB(1,23) + crLeftHandSideMatrix80*rB(2,23) + crLeftHandSideMatrix81*rB(3,23) + crLeftHandSideMatrix82*rB(4,23) + crLeftHandSideMatrix83*rB(5,23));

// // rLeftHandSideMatrix(14,0)+=(crLeftHandSideMatrix84*rB(0,0) + crLeftHandSideMatrix85*rB(1,0) + crLeftHandSideMatrix86*rB(2,0) + crLeftHandSideMatrix87*rB(3,0) + crLeftHandSideMatrix88*rB(4,0) + crLeftHandSideMatrix89*rB(5,0));
// // rLeftHandSideMatrix(14,1)+=(crLeftHandSideMatrix84*rB(0,1) + crLeftHandSideMatrix85*rB(1,1) + crLeftHandSideMatrix86*rB(2,1) + crLeftHandSideMatrix87*rB(3,1) + crLeftHandSideMatrix88*rB(4,1) + crLeftHandSideMatrix89*rB(5,1));
// // rLeftHandSideMatrix(14,2)+=(crLeftHandSideMatrix84*rB(0,2) + crLeftHandSideMatrix85*rB(1,2) + crLeftHandSideMatrix86*rB(2,2) + crLeftHandSideMatrix87*rB(3,2) + crLeftHandSideMatrix88*rB(4,2) + crLeftHandSideMatrix89*rB(5,2));
// // rLeftHandSideMatrix(14,3)+=(crLeftHandSideMatrix84*rB(0,3) + crLeftHandSideMatrix85*rB(1,3) + crLeftHandSideMatrix86*rB(2,3) + crLeftHandSideMatrix87*rB(3,3) + crLeftHandSideMatrix88*rB(4,3) + crLeftHandSideMatrix89*rB(5,3));
// // rLeftHandSideMatrix(14,4)+=(crLeftHandSideMatrix84*rB(0,4) + crLeftHandSideMatrix85*rB(1,4) + crLeftHandSideMatrix86*rB(2,4) + crLeftHandSideMatrix87*rB(3,4) + crLeftHandSideMatrix88*rB(4,4) + crLeftHandSideMatrix89*rB(5,4));
// // rLeftHandSideMatrix(14,5)+=(crLeftHandSideMatrix84*rB(0,5) + crLeftHandSideMatrix85*rB(1,5) + crLeftHandSideMatrix86*rB(2,5) + crLeftHandSideMatrix87*rB(3,5) + crLeftHandSideMatrix88*rB(4,5) + crLeftHandSideMatrix89*rB(5,5));
// // rLeftHandSideMatrix(14,6)+=(crLeftHandSideMatrix84*rB(0,6) + crLeftHandSideMatrix85*rB(1,6) + crLeftHandSideMatrix86*rB(2,6) + crLeftHandSideMatrix87*rB(3,6) + crLeftHandSideMatrix88*rB(4,6) + crLeftHandSideMatrix89*rB(5,6));
// // rLeftHandSideMatrix(14,7)+=(crLeftHandSideMatrix84*rB(0,7) + crLeftHandSideMatrix85*rB(1,7) + crLeftHandSideMatrix86*rB(2,7) + crLeftHandSideMatrix87*rB(3,7) + crLeftHandSideMatrix88*rB(4,7) + crLeftHandSideMatrix89*rB(5,7));
// // rLeftHandSideMatrix(14,8)+=(crLeftHandSideMatrix84*rB(0,8) + crLeftHandSideMatrix85*rB(1,8) + crLeftHandSideMatrix86*rB(2,8) + crLeftHandSideMatrix87*rB(3,8) + crLeftHandSideMatrix88*rB(4,8) + crLeftHandSideMatrix89*rB(5,8));
// // rLeftHandSideMatrix(14,9)+=(crLeftHandSideMatrix84*rB(0,9) + crLeftHandSideMatrix85*rB(1,9) + crLeftHandSideMatrix86*rB(2,9) + crLeftHandSideMatrix87*rB(3,9) + crLeftHandSideMatrix88*rB(4,9) + crLeftHandSideMatrix89*rB(5,9));
// // rLeftHandSideMatrix(14,10)+=(crLeftHandSideMatrix84*rB(0,10) + crLeftHandSideMatrix85*rB(1,10) + crLeftHandSideMatrix86*rB(2,10) + crLeftHandSideMatrix87*rB(3,10) + crLeftHandSideMatrix88*rB(4,10) + crLeftHandSideMatrix89*rB(5,10));
// // rLeftHandSideMatrix(14,11)+=(crLeftHandSideMatrix84*rB(0,11) + crLeftHandSideMatrix85*rB(1,11) + crLeftHandSideMatrix86*rB(2,11) + crLeftHandSideMatrix87*rB(3,11) + crLeftHandSideMatrix88*rB(4,11) + crLeftHandSideMatrix89*rB(5,11));
// // rLeftHandSideMatrix(14,12)+=(crLeftHandSideMatrix84*rB(0,12) + crLeftHandSideMatrix85*rB(1,12) + crLeftHandSideMatrix86*rB(2,12) + crLeftHandSideMatrix87*rB(3,12) + crLeftHandSideMatrix88*rB(4,12) + crLeftHandSideMatrix89*rB(5,12));
// // rLeftHandSideMatrix(14,13)+=(crLeftHandSideMatrix84*rB(0,13) + crLeftHandSideMatrix85*rB(1,13) + crLeftHandSideMatrix86*rB(2,13) + crLeftHandSideMatrix87*rB(3,13) + crLeftHandSideMatrix88*rB(4,13) + crLeftHandSideMatrix89*rB(5,13));
// rLeftHandSideMatrix(14,14)+=(crLeftHandSideMatrix84*rB(0,14) + crLeftHandSideMatrix85*rB(1,14) + crLeftHandSideMatrix86*rB(2,14) + crLeftHandSideMatrix87*rB(3,14) + crLeftHandSideMatrix88*rB(4,14) + crLeftHandSideMatrix89*rB(5,14));
// rLeftHandSideMatrix(14,15)+=(crLeftHandSideMatrix84*rB(0,15) + crLeftHandSideMatrix85*rB(1,15) + crLeftHandSideMatrix86*rB(2,15) + crLeftHandSideMatrix87*rB(3,15) + crLeftHandSideMatrix88*rB(4,15) + crLeftHandSideMatrix89*rB(5,15));
// rLeftHandSideMatrix(14,16)+=(crLeftHandSideMatrix84*rB(0,16) + crLeftHandSideMatrix85*rB(1,16) + crLeftHandSideMatrix86*rB(2,16) + crLeftHandSideMatrix87*rB(3,16) + crLeftHandSideMatrix88*rB(4,16) + crLeftHandSideMatrix89*rB(5,16));
// rLeftHandSideMatrix(14,17)+=(crLeftHandSideMatrix84*rB(0,17) + crLeftHandSideMatrix85*rB(1,17) + crLeftHandSideMatrix86*rB(2,17) + crLeftHandSideMatrix87*rB(3,17) + crLeftHandSideMatrix88*rB(4,17) + crLeftHandSideMatrix89*rB(5,17));
// rLeftHandSideMatrix(14,18)+=(crLeftHandSideMatrix84*rB(0,18) + crLeftHandSideMatrix85*rB(1,18) + crLeftHandSideMatrix86*rB(2,18) + crLeftHandSideMatrix87*rB(3,18) + crLeftHandSideMatrix88*rB(4,18) + crLeftHandSideMatrix89*rB(5,18));
// rLeftHandSideMatrix(14,19)+=(crLeftHandSideMatrix84*rB(0,19) + crLeftHandSideMatrix85*rB(1,19) + crLeftHandSideMatrix86*rB(2,19) + crLeftHandSideMatrix87*rB(3,19) + crLeftHandSideMatrix88*rB(4,19) + crLeftHandSideMatrix89*rB(5,19));
// rLeftHandSideMatrix(14,20)+=(crLeftHandSideMatrix84*rB(0,20) + crLeftHandSideMatrix85*rB(1,20) + crLeftHandSideMatrix86*rB(2,20) + crLeftHandSideMatrix87*rB(3,20) + crLeftHandSideMatrix88*rB(4,20) + crLeftHandSideMatrix89*rB(5,20));
// rLeftHandSideMatrix(14,21)+=(crLeftHandSideMatrix84*rB(0,21) + crLeftHandSideMatrix85*rB(1,21) + crLeftHandSideMatrix86*rB(2,21) + crLeftHandSideMatrix87*rB(3,21) + crLeftHandSideMatrix88*rB(4,21) + crLeftHandSideMatrix89*rB(5,21));
// rLeftHandSideMatrix(14,22)+=(crLeftHandSideMatrix84*rB(0,22) + crLeftHandSideMatrix85*rB(1,22) + crLeftHandSideMatrix86*rB(2,22) + crLeftHandSideMatrix87*rB(3,22) + crLeftHandSideMatrix88*rB(4,22) + crLeftHandSideMatrix89*rB(5,22));
// rLeftHandSideMatrix(14,23)+=(crLeftHandSideMatrix84*rB(0,23) + crLeftHandSideMatrix85*rB(1,23) + crLeftHandSideMatrix86*rB(2,23) + crLeftHandSideMatrix87*rB(3,23) + crLeftHandSideMatrix88*rB(4,23) + crLeftHandSideMatrix89*rB(5,23));

// // rLeftHandSideMatrix(15,0)+=(crLeftHandSideMatrix90*rB(0,0) + crLeftHandSideMatrix91*rB(1,0) + crLeftHandSideMatrix92*rB(2,0) + crLeftHandSideMatrix93*rB(3,0) + crLeftHandSideMatrix94*rB(4,0) + crLeftHandSideMatrix95*rB(5,0));
// // rLeftHandSideMatrix(15,1)+=(crLeftHandSideMatrix90*rB(0,1) + crLeftHandSideMatrix91*rB(1,1) + crLeftHandSideMatrix92*rB(2,1) + crLeftHandSideMatrix93*rB(3,1) + crLeftHandSideMatrix94*rB(4,1) + crLeftHandSideMatrix95*rB(5,1));
// // rLeftHandSideMatrix(15,2)+=(crLeftHandSideMatrix90*rB(0,2) + crLeftHandSideMatrix91*rB(1,2) + crLeftHandSideMatrix92*rB(2,2) + crLeftHandSideMatrix93*rB(3,2) + crLeftHandSideMatrix94*rB(4,2) + crLeftHandSideMatrix95*rB(5,2));
// // rLeftHandSideMatrix(15,3)+=(crLeftHandSideMatrix90*rB(0,3) + crLeftHandSideMatrix91*rB(1,3) + crLeftHandSideMatrix92*rB(2,3) + crLeftHandSideMatrix93*rB(3,3) + crLeftHandSideMatrix94*rB(4,3) + crLeftHandSideMatrix95*rB(5,3));
// // rLeftHandSideMatrix(15,4)+=(crLeftHandSideMatrix90*rB(0,4) + crLeftHandSideMatrix91*rB(1,4) + crLeftHandSideMatrix92*rB(2,4) + crLeftHandSideMatrix93*rB(3,4) + crLeftHandSideMatrix94*rB(4,4) + crLeftHandSideMatrix95*rB(5,4));
// // rLeftHandSideMatrix(15,5)+=(crLeftHandSideMatrix90*rB(0,5) + crLeftHandSideMatrix91*rB(1,5) + crLeftHandSideMatrix92*rB(2,5) + crLeftHandSideMatrix93*rB(3,5) + crLeftHandSideMatrix94*rB(4,5) + crLeftHandSideMatrix95*rB(5,5));
// // rLeftHandSideMatrix(15,6)+=(crLeftHandSideMatrix90*rB(0,6) + crLeftHandSideMatrix91*rB(1,6) + crLeftHandSideMatrix92*rB(2,6) + crLeftHandSideMatrix93*rB(3,6) + crLeftHandSideMatrix94*rB(4,6) + crLeftHandSideMatrix95*rB(5,6));
// // rLeftHandSideMatrix(15,7)+=(crLeftHandSideMatrix90*rB(0,7) + crLeftHandSideMatrix91*rB(1,7) + crLeftHandSideMatrix92*rB(2,7) + crLeftHandSideMatrix93*rB(3,7) + crLeftHandSideMatrix94*rB(4,7) + crLeftHandSideMatrix95*rB(5,7));
// // rLeftHandSideMatrix(15,8)+=(crLeftHandSideMatrix90*rB(0,8) + crLeftHandSideMatrix91*rB(1,8) + crLeftHandSideMatrix92*rB(2,8) + crLeftHandSideMatrix93*rB(3,8) + crLeftHandSideMatrix94*rB(4,8) + crLeftHandSideMatrix95*rB(5,8));
// // rLeftHandSideMatrix(15,9)+=(crLeftHandSideMatrix90*rB(0,9) + crLeftHandSideMatrix91*rB(1,9) + crLeftHandSideMatrix92*rB(2,9) + crLeftHandSideMatrix93*rB(3,9) + crLeftHandSideMatrix94*rB(4,9) + crLeftHandSideMatrix95*rB(5,9));
// // rLeftHandSideMatrix(15,10)+=(crLeftHandSideMatrix90*rB(0,10) + crLeftHandSideMatrix91*rB(1,10) + crLeftHandSideMatrix92*rB(2,10) + crLeftHandSideMatrix93*rB(3,10) + crLeftHandSideMatrix94*rB(4,10) + crLeftHandSideMatrix95*rB(5,10));
// // rLeftHandSideMatrix(15,11)+=(crLeftHandSideMatrix90*rB(0,11) + crLeftHandSideMatrix91*rB(1,11) + crLeftHandSideMatrix92*rB(2,11) + crLeftHandSideMatrix93*rB(3,11) + crLeftHandSideMatrix94*rB(4,11) + crLeftHandSideMatrix95*rB(5,11));
// // rLeftHandSideMatrix(15,12)+=(crLeftHandSideMatrix90*rB(0,12) + crLeftHandSideMatrix91*rB(1,12) + crLeftHandSideMatrix92*rB(2,12) + crLeftHandSideMatrix93*rB(3,12) + crLeftHandSideMatrix94*rB(4,12) + crLeftHandSideMatrix95*rB(5,12));
// // rLeftHandSideMatrix(15,13)+=(crLeftHandSideMatrix90*rB(0,13) + crLeftHandSideMatrix91*rB(1,13) + crLeftHandSideMatrix92*rB(2,13) + crLeftHandSideMatrix93*rB(3,13) + crLeftHandSideMatrix94*rB(4,13) + crLeftHandSideMatrix95*rB(5,13));
// // rLeftHandSideMatrix(15,14)+=(crLeftHandSideMatrix90*rB(0,14) + crLeftHandSideMatrix91*rB(1,14) + crLeftHandSideMatrix92*rB(2,14) + crLeftHandSideMatrix93*rB(3,14) + crLeftHandSideMatrix94*rB(4,14) + crLeftHandSideMatrix95*rB(5,14));
// rLeftHandSideMatrix(15,15)+=(crLeftHandSideMatrix90*rB(0,15) + crLeftHandSideMatrix91*rB(1,15) + crLeftHandSideMatrix92*rB(2,15) + crLeftHandSideMatrix93*rB(3,15) + crLeftHandSideMatrix94*rB(4,15) + crLeftHandSideMatrix95*rB(5,15));
// rLeftHandSideMatrix(15,16)+=(crLeftHandSideMatrix90*rB(0,16) + crLeftHandSideMatrix91*rB(1,16) + crLeftHandSideMatrix92*rB(2,16) + crLeftHandSideMatrix93*rB(3,16) + crLeftHandSideMatrix94*rB(4,16) + crLeftHandSideMatrix95*rB(5,16));
// rLeftHandSideMatrix(15,17)+=(crLeftHandSideMatrix90*rB(0,17) + crLeftHandSideMatrix91*rB(1,17) + crLeftHandSideMatrix92*rB(2,17) + crLeftHandSideMatrix93*rB(3,17) + crLeftHandSideMatrix94*rB(4,17) + crLeftHandSideMatrix95*rB(5,17));
// rLeftHandSideMatrix(15,18)+=(crLeftHandSideMatrix90*rB(0,18) + crLeftHandSideMatrix91*rB(1,18) + crLeftHandSideMatrix92*rB(2,18) + crLeftHandSideMatrix93*rB(3,18) + crLeftHandSideMatrix94*rB(4,18) + crLeftHandSideMatrix95*rB(5,18));
// rLeftHandSideMatrix(15,19)+=(crLeftHandSideMatrix90*rB(0,19) + crLeftHandSideMatrix91*rB(1,19) + crLeftHandSideMatrix92*rB(2,19) + crLeftHandSideMatrix93*rB(3,19) + crLeftHandSideMatrix94*rB(4,19) + crLeftHandSideMatrix95*rB(5,19));
// rLeftHandSideMatrix(15,20)+=(crLeftHandSideMatrix90*rB(0,20) + crLeftHandSideMatrix91*rB(1,20) + crLeftHandSideMatrix92*rB(2,20) + crLeftHandSideMatrix93*rB(3,20) + crLeftHandSideMatrix94*rB(4,20) + crLeftHandSideMatrix95*rB(5,20));
// rLeftHandSideMatrix(15,21)+=(crLeftHandSideMatrix90*rB(0,21) + crLeftHandSideMatrix91*rB(1,21) + crLeftHandSideMatrix92*rB(2,21) + crLeftHandSideMatrix93*rB(3,21) + crLeftHandSideMatrix94*rB(4,21) + crLeftHandSideMatrix95*rB(5,21));
// rLeftHandSideMatrix(15,22)+=(crLeftHandSideMatrix90*rB(0,22) + crLeftHandSideMatrix91*rB(1,22) + crLeftHandSideMatrix92*rB(2,22) + crLeftHandSideMatrix93*rB(3,22) + crLeftHandSideMatrix94*rB(4,22) + crLeftHandSideMatrix95*rB(5,22));
// rLeftHandSideMatrix(15,23)+=(crLeftHandSideMatrix90*rB(0,23) + crLeftHandSideMatrix91*rB(1,23) + crLeftHandSideMatrix92*rB(2,23) + crLeftHandSideMatrix93*rB(3,23) + crLeftHandSideMatrix94*rB(4,23) + crLeftHandSideMatrix95*rB(5,23));

// // rLeftHandSideMatrix(16,0)+=(crLeftHandSideMatrix100*rB(4,0) + crLeftHandSideMatrix101*rB(5,0) + crLeftHandSideMatrix96*rB(0,0) + crLeftHandSideMatrix97*rB(1,0) + crLeftHandSideMatrix98*rB(2,0) + crLeftHandSideMatrix99*rB(3,0));
// // rLeftHandSideMatrix(16,1)+=(crLeftHandSideMatrix100*rB(4,1) + crLeftHandSideMatrix101*rB(5,1) + crLeftHandSideMatrix96*rB(0,1) + crLeftHandSideMatrix97*rB(1,1) + crLeftHandSideMatrix98*rB(2,1) + crLeftHandSideMatrix99*rB(3,1));
// // rLeftHandSideMatrix(16,2)+=(crLeftHandSideMatrix100*rB(4,2) + crLeftHandSideMatrix101*rB(5,2) + crLeftHandSideMatrix96*rB(0,2) + crLeftHandSideMatrix97*rB(1,2) + crLeftHandSideMatrix98*rB(2,2) + crLeftHandSideMatrix99*rB(3,2));
// // rLeftHandSideMatrix(16,3)+=(crLeftHandSideMatrix100*rB(4,3) + crLeftHandSideMatrix101*rB(5,3) + crLeftHandSideMatrix96*rB(0,3) + crLeftHandSideMatrix97*rB(1,3) + crLeftHandSideMatrix98*rB(2,3) + crLeftHandSideMatrix99*rB(3,3));
// // rLeftHandSideMatrix(16,4)+=(crLeftHandSideMatrix100*rB(4,4) + crLeftHandSideMatrix101*rB(5,4) + crLeftHandSideMatrix96*rB(0,4) + crLeftHandSideMatrix97*rB(1,4) + crLeftHandSideMatrix98*rB(2,4) + crLeftHandSideMatrix99*rB(3,4));
// // rLeftHandSideMatrix(16,5)+=(crLeftHandSideMatrix100*rB(4,5) + crLeftHandSideMatrix101*rB(5,5) + crLeftHandSideMatrix96*rB(0,5) + crLeftHandSideMatrix97*rB(1,5) + crLeftHandSideMatrix98*rB(2,5) + crLeftHandSideMatrix99*rB(3,5));
// // rLeftHandSideMatrix(16,6)+=(crLeftHandSideMatrix100*rB(4,6) + crLeftHandSideMatrix101*rB(5,6) + crLeftHandSideMatrix96*rB(0,6) + crLeftHandSideMatrix97*rB(1,6) + crLeftHandSideMatrix98*rB(2,6) + crLeftHandSideMatrix99*rB(3,6));
// // rLeftHandSideMatrix(16,7)+=(crLeftHandSideMatrix100*rB(4,7) + crLeftHandSideMatrix101*rB(5,7) + crLeftHandSideMatrix96*rB(0,7) + crLeftHandSideMatrix97*rB(1,7) + crLeftHandSideMatrix98*rB(2,7) + crLeftHandSideMatrix99*rB(3,7));
// // rLeftHandSideMatrix(16,8)+=(crLeftHandSideMatrix100*rB(4,8) + crLeftHandSideMatrix101*rB(5,8) + crLeftHandSideMatrix96*rB(0,8) + crLeftHandSideMatrix97*rB(1,8) + crLeftHandSideMatrix98*rB(2,8) + crLeftHandSideMatrix99*rB(3,8));
// // rLeftHandSideMatrix(16,9)+=(crLeftHandSideMatrix100*rB(4,9) + crLeftHandSideMatrix101*rB(5,9) + crLeftHandSideMatrix96*rB(0,9) + crLeftHandSideMatrix97*rB(1,9) + crLeftHandSideMatrix98*rB(2,9) + crLeftHandSideMatrix99*rB(3,9));
// // rLeftHandSideMatrix(16,10)+=(crLeftHandSideMatrix100*rB(4,10) + crLeftHandSideMatrix101*rB(5,10) + crLeftHandSideMatrix96*rB(0,10) + crLeftHandSideMatrix97*rB(1,10) + crLeftHandSideMatrix98*rB(2,10) + crLeftHandSideMatrix99*rB(3,10));
// // rLeftHandSideMatrix(16,11)+=(crLeftHandSideMatrix100*rB(4,11) + crLeftHandSideMatrix101*rB(5,11) + crLeftHandSideMatrix96*rB(0,11) + crLeftHandSideMatrix97*rB(1,11) + crLeftHandSideMatrix98*rB(2,11) + crLeftHandSideMatrix99*rB(3,11));
// // rLeftHandSideMatrix(16,12)+=(crLeftHandSideMatrix100*rB(4,12) + crLeftHandSideMatrix101*rB(5,12) + crLeftHandSideMatrix96*rB(0,12) + crLeftHandSideMatrix97*rB(1,12) + crLeftHandSideMatrix98*rB(2,12) + crLeftHandSideMatrix99*rB(3,12));
// // rLeftHandSideMatrix(16,13)+=(crLeftHandSideMatrix100*rB(4,13) + crLeftHandSideMatrix101*rB(5,13) + crLeftHandSideMatrix96*rB(0,13) + crLeftHandSideMatrix97*rB(1,13) + crLeftHandSideMatrix98*rB(2,13) + crLeftHandSideMatrix99*rB(3,13));
// // rLeftHandSideMatrix(16,14)+=(crLeftHandSideMatrix100*rB(4,14) + crLeftHandSideMatrix101*rB(5,14) + crLeftHandSideMatrix96*rB(0,14) + crLeftHandSideMatrix97*rB(1,14) + crLeftHandSideMatrix98*rB(2,14) + crLeftHandSideMatrix99*rB(3,14));
// // rLeftHandSideMatrix(16,15)+=(crLeftHandSideMatrix100*rB(4,15) + crLeftHandSideMatrix101*rB(5,15) + crLeftHandSideMatrix96*rB(0,15) + crLeftHandSideMatrix97*rB(1,15) + crLeftHandSideMatrix98*rB(2,15) + crLeftHandSideMatrix99*rB(3,15));
// rLeftHandSideMatrix(16,16)+=(crLeftHandSideMatrix100*rB(4,16) + crLeftHandSideMatrix101*rB(5,16) + crLeftHandSideMatrix96*rB(0,16) + crLeftHandSideMatrix97*rB(1,16) + crLeftHandSideMatrix98*rB(2,16) + crLeftHandSideMatrix99*rB(3,16));
// rLeftHandSideMatrix(16,17)+=(crLeftHandSideMatrix100*rB(4,17) + crLeftHandSideMatrix101*rB(5,17) + crLeftHandSideMatrix96*rB(0,17) + crLeftHandSideMatrix97*rB(1,17) + crLeftHandSideMatrix98*rB(2,17) + crLeftHandSideMatrix99*rB(3,17));
// rLeftHandSideMatrix(16,18)+=(crLeftHandSideMatrix100*rB(4,18) + crLeftHandSideMatrix101*rB(5,18) + crLeftHandSideMatrix96*rB(0,18) + crLeftHandSideMatrix97*rB(1,18) + crLeftHandSideMatrix98*rB(2,18) + crLeftHandSideMatrix99*rB(3,18));
// rLeftHandSideMatrix(16,19)+=(crLeftHandSideMatrix100*rB(4,19) + crLeftHandSideMatrix101*rB(5,19) + crLeftHandSideMatrix96*rB(0,19) + crLeftHandSideMatrix97*rB(1,19) + crLeftHandSideMatrix98*rB(2,19) + crLeftHandSideMatrix99*rB(3,19));
// rLeftHandSideMatrix(16,20)+=(crLeftHandSideMatrix100*rB(4,20) + crLeftHandSideMatrix101*rB(5,20) + crLeftHandSideMatrix96*rB(0,20) + crLeftHandSideMatrix97*rB(1,20) + crLeftHandSideMatrix98*rB(2,20) + crLeftHandSideMatrix99*rB(3,20));
// rLeftHandSideMatrix(16,21)+=(crLeftHandSideMatrix100*rB(4,21) + crLeftHandSideMatrix101*rB(5,21) + crLeftHandSideMatrix96*rB(0,21) + crLeftHandSideMatrix97*rB(1,21) + crLeftHandSideMatrix98*rB(2,21) + crLeftHandSideMatrix99*rB(3,21));
// rLeftHandSideMatrix(16,22)+=(crLeftHandSideMatrix100*rB(4,22) + crLeftHandSideMatrix101*rB(5,22) + crLeftHandSideMatrix96*rB(0,22) + crLeftHandSideMatrix97*rB(1,22) + crLeftHandSideMatrix98*rB(2,22) + crLeftHandSideMatrix99*rB(3,22));
// rLeftHandSideMatrix(16,23)+=(crLeftHandSideMatrix100*rB(4,23) + crLeftHandSideMatrix101*rB(5,23) + crLeftHandSideMatrix96*rB(0,23) + crLeftHandSideMatrix97*rB(1,23) + crLeftHandSideMatrix98*rB(2,23) + crLeftHandSideMatrix99*rB(3,23));

// // rLeftHandSideMatrix(17,0)+=(crLeftHandSideMatrix102*rB(0,0) + crLeftHandSideMatrix103*rB(1,0) + crLeftHandSideMatrix104*rB(2,0) + crLeftHandSideMatrix105*rB(3,0) + crLeftHandSideMatrix106*rB(4,0) + crLeftHandSideMatrix107*rB(5,0));
// // rLeftHandSideMatrix(17,1)+=(crLeftHandSideMatrix102*rB(0,1) + crLeftHandSideMatrix103*rB(1,1) + crLeftHandSideMatrix104*rB(2,1) + crLeftHandSideMatrix105*rB(3,1) + crLeftHandSideMatrix106*rB(4,1) + crLeftHandSideMatrix107*rB(5,1));
// // rLeftHandSideMatrix(17,2)+=(crLeftHandSideMatrix102*rB(0,2) + crLeftHandSideMatrix103*rB(1,2) + crLeftHandSideMatrix104*rB(2,2) + crLeftHandSideMatrix105*rB(3,2) + crLeftHandSideMatrix106*rB(4,2) + crLeftHandSideMatrix107*rB(5,2));
// // rLeftHandSideMatrix(17,3)+=(crLeftHandSideMatrix102*rB(0,3) + crLeftHandSideMatrix103*rB(1,3) + crLeftHandSideMatrix104*rB(2,3) + crLeftHandSideMatrix105*rB(3,3) + crLeftHandSideMatrix106*rB(4,3) + crLeftHandSideMatrix107*rB(5,3));
// // rLeftHandSideMatrix(17,4)+=(crLeftHandSideMatrix102*rB(0,4) + crLeftHandSideMatrix103*rB(1,4) + crLeftHandSideMatrix104*rB(2,4) + crLeftHandSideMatrix105*rB(3,4) + crLeftHandSideMatrix106*rB(4,4) + crLeftHandSideMatrix107*rB(5,4));
// // rLeftHandSideMatrix(17,5)+=(crLeftHandSideMatrix102*rB(0,5) + crLeftHandSideMatrix103*rB(1,5) + crLeftHandSideMatrix104*rB(2,5) + crLeftHandSideMatrix105*rB(3,5) + crLeftHandSideMatrix106*rB(4,5) + crLeftHandSideMatrix107*rB(5,5));
// // rLeftHandSideMatrix(17,6)+=(crLeftHandSideMatrix102*rB(0,6) + crLeftHandSideMatrix103*rB(1,6) + crLeftHandSideMatrix104*rB(2,6) + crLeftHandSideMatrix105*rB(3,6) + crLeftHandSideMatrix106*rB(4,6) + crLeftHandSideMatrix107*rB(5,6));
// // rLeftHandSideMatrix(17,7)+=(crLeftHandSideMatrix102*rB(0,7) + crLeftHandSideMatrix103*rB(1,7) + crLeftHandSideMatrix104*rB(2,7) + crLeftHandSideMatrix105*rB(3,7) + crLeftHandSideMatrix106*rB(4,7) + crLeftHandSideMatrix107*rB(5,7));
// // rLeftHandSideMatrix(17,8)+=(crLeftHandSideMatrix102*rB(0,8) + crLeftHandSideMatrix103*rB(1,8) + crLeftHandSideMatrix104*rB(2,8) + crLeftHandSideMatrix105*rB(3,8) + crLeftHandSideMatrix106*rB(4,8) + crLeftHandSideMatrix107*rB(5,8));
// // rLeftHandSideMatrix(17,9)+=(crLeftHandSideMatrix102*rB(0,9) + crLeftHandSideMatrix103*rB(1,9) + crLeftHandSideMatrix104*rB(2,9) + crLeftHandSideMatrix105*rB(3,9) + crLeftHandSideMatrix106*rB(4,9) + crLeftHandSideMatrix107*rB(5,9));
// // rLeftHandSideMatrix(17,10)+=(crLeftHandSideMatrix102*rB(0,10) + crLeftHandSideMatrix103*rB(1,10) + crLeftHandSideMatrix104*rB(2,10) + crLeftHandSideMatrix105*rB(3,10) + crLeftHandSideMatrix106*rB(4,10) + crLeftHandSideMatrix107*rB(5,10));
// // rLeftHandSideMatrix(17,11)+=(crLeftHandSideMatrix102*rB(0,11) + crLeftHandSideMatrix103*rB(1,11) + crLeftHandSideMatrix104*rB(2,11) + crLeftHandSideMatrix105*rB(3,11) + crLeftHandSideMatrix106*rB(4,11) + crLeftHandSideMatrix107*rB(5,11));
// // rLeftHandSideMatrix(17,12)+=(crLeftHandSideMatrix102*rB(0,12) + crLeftHandSideMatrix103*rB(1,12) + crLeftHandSideMatrix104*rB(2,12) + crLeftHandSideMatrix105*rB(3,12) + crLeftHandSideMatrix106*rB(4,12) + crLeftHandSideMatrix107*rB(5,12));
// // rLeftHandSideMatrix(17,13)+=(crLeftHandSideMatrix102*rB(0,13) + crLeftHandSideMatrix103*rB(1,13) + crLeftHandSideMatrix104*rB(2,13) + crLeftHandSideMatrix105*rB(3,13) + crLeftHandSideMatrix106*rB(4,13) + crLeftHandSideMatrix107*rB(5,13));
// // rLeftHandSideMatrix(17,14)+=(crLeftHandSideMatrix102*rB(0,14) + crLeftHandSideMatrix103*rB(1,14) + crLeftHandSideMatrix104*rB(2,14) + crLeftHandSideMatrix105*rB(3,14) + crLeftHandSideMatrix106*rB(4,14) + crLeftHandSideMatrix107*rB(5,14));
// // rLeftHandSideMatrix(17,15)+=(crLeftHandSideMatrix102*rB(0,15) + crLeftHandSideMatrix103*rB(1,15) + crLeftHandSideMatrix104*rB(2,15) + crLeftHandSideMatrix105*rB(3,15) + crLeftHandSideMatrix106*rB(4,15) + crLeftHandSideMatrix107*rB(5,15));
// // rLeftHandSideMatrix(17,16)+=(crLeftHandSideMatrix102*rB(0,16) + crLeftHandSideMatrix103*rB(1,16) + crLeftHandSideMatrix104*rB(2,16) + crLeftHandSideMatrix105*rB(3,16) + crLeftHandSideMatrix106*rB(4,16) + crLeftHandSideMatrix107*rB(5,16));
// rLeftHandSideMatrix(17,17)+=(crLeftHandSideMatrix102*rB(0,17) + crLeftHandSideMatrix103*rB(1,17) + crLeftHandSideMatrix104*rB(2,17) + crLeftHandSideMatrix105*rB(3,17) + crLeftHandSideMatrix106*rB(4,17) + crLeftHandSideMatrix107*rB(5,17));
// rLeftHandSideMatrix(17,18)+=(crLeftHandSideMatrix102*rB(0,18) + crLeftHandSideMatrix103*rB(1,18) + crLeftHandSideMatrix104*rB(2,18) + crLeftHandSideMatrix105*rB(3,18) + crLeftHandSideMatrix106*rB(4,18) + crLeftHandSideMatrix107*rB(5,18));
// rLeftHandSideMatrix(17,19)+=(crLeftHandSideMatrix102*rB(0,19) + crLeftHandSideMatrix103*rB(1,19) + crLeftHandSideMatrix104*rB(2,19) + crLeftHandSideMatrix105*rB(3,19) + crLeftHandSideMatrix106*rB(4,19) + crLeftHandSideMatrix107*rB(5,19));
// rLeftHandSideMatrix(17,20)+=(crLeftHandSideMatrix102*rB(0,20) + crLeftHandSideMatrix103*rB(1,20) + crLeftHandSideMatrix104*rB(2,20) + crLeftHandSideMatrix105*rB(3,20) + crLeftHandSideMatrix106*rB(4,20) + crLeftHandSideMatrix107*rB(5,20));
// rLeftHandSideMatrix(17,21)+=(crLeftHandSideMatrix102*rB(0,21) + crLeftHandSideMatrix103*rB(1,21) + crLeftHandSideMatrix104*rB(2,21) + crLeftHandSideMatrix105*rB(3,21) + crLeftHandSideMatrix106*rB(4,21) + crLeftHandSideMatrix107*rB(5,21));
// rLeftHandSideMatrix(17,22)+=(crLeftHandSideMatrix102*rB(0,22) + crLeftHandSideMatrix103*rB(1,22) + crLeftHandSideMatrix104*rB(2,22) + crLeftHandSideMatrix105*rB(3,22) + crLeftHandSideMatrix106*rB(4,22) + crLeftHandSideMatrix107*rB(5,22));
// rLeftHandSideMatrix(17,23)+=(crLeftHandSideMatrix102*rB(0,23) + crLeftHandSideMatrix103*rB(1,23) + crLeftHandSideMatrix104*rB(2,23) + crLeftHandSideMatrix105*rB(3,23) + crLeftHandSideMatrix106*rB(4,23) + crLeftHandSideMatrix107*rB(5,23));

// // rLeftHandSideMatrix(18,0)+=(crLeftHandSideMatrix108*rB(0,0) + crLeftHandSideMatrix109*rB(1,0) + crLeftHandSideMatrix110*rB(2,0) + crLeftHandSideMatrix111*rB(3,0) + crLeftHandSideMatrix112*rB(4,0) + crLeftHandSideMatrix113*rB(5,0));
// // rLeftHandSideMatrix(18,1)+=(crLeftHandSideMatrix108*rB(0,1) + crLeftHandSideMatrix109*rB(1,1) + crLeftHandSideMatrix110*rB(2,1) + crLeftHandSideMatrix111*rB(3,1) + crLeftHandSideMatrix112*rB(4,1) + crLeftHandSideMatrix113*rB(5,1));
// // rLeftHandSideMatrix(18,2)+=(crLeftHandSideMatrix108*rB(0,2) + crLeftHandSideMatrix109*rB(1,2) + crLeftHandSideMatrix110*rB(2,2) + crLeftHandSideMatrix111*rB(3,2) + crLeftHandSideMatrix112*rB(4,2) + crLeftHandSideMatrix113*rB(5,2));
// // rLeftHandSideMatrix(18,3)+=(crLeftHandSideMatrix108*rB(0,3) + crLeftHandSideMatrix109*rB(1,3) + crLeftHandSideMatrix110*rB(2,3) + crLeftHandSideMatrix111*rB(3,3) + crLeftHandSideMatrix112*rB(4,3) + crLeftHandSideMatrix113*rB(5,3));
// // rLeftHandSideMatrix(18,4)+=(crLeftHandSideMatrix108*rB(0,4) + crLeftHandSideMatrix109*rB(1,4) + crLeftHandSideMatrix110*rB(2,4) + crLeftHandSideMatrix111*rB(3,4) + crLeftHandSideMatrix112*rB(4,4) + crLeftHandSideMatrix113*rB(5,4));
// // rLeftHandSideMatrix(18,5)+=(crLeftHandSideMatrix108*rB(0,5) + crLeftHandSideMatrix109*rB(1,5) + crLeftHandSideMatrix110*rB(2,5) + crLeftHandSideMatrix111*rB(3,5) + crLeftHandSideMatrix112*rB(4,5) + crLeftHandSideMatrix113*rB(5,5));
// // rLeftHandSideMatrix(18,6)+=(crLeftHandSideMatrix108*rB(0,6) + crLeftHandSideMatrix109*rB(1,6) + crLeftHandSideMatrix110*rB(2,6) + crLeftHandSideMatrix111*rB(3,6) + crLeftHandSideMatrix112*rB(4,6) + crLeftHandSideMatrix113*rB(5,6));
// // rLeftHandSideMatrix(18,7)+=(crLeftHandSideMatrix108*rB(0,7) + crLeftHandSideMatrix109*rB(1,7) + crLeftHandSideMatrix110*rB(2,7) + crLeftHandSideMatrix111*rB(3,7) + crLeftHandSideMatrix112*rB(4,7) + crLeftHandSideMatrix113*rB(5,7));
// // rLeftHandSideMatrix(18,8)+=(crLeftHandSideMatrix108*rB(0,8) + crLeftHandSideMatrix109*rB(1,8) + crLeftHandSideMatrix110*rB(2,8) + crLeftHandSideMatrix111*rB(3,8) + crLeftHandSideMatrix112*rB(4,8) + crLeftHandSideMatrix113*rB(5,8));
// // rLeftHandSideMatrix(18,9)+=(crLeftHandSideMatrix108*rB(0,9) + crLeftHandSideMatrix109*rB(1,9) + crLeftHandSideMatrix110*rB(2,9) + crLeftHandSideMatrix111*rB(3,9) + crLeftHandSideMatrix112*rB(4,9) + crLeftHandSideMatrix113*rB(5,9));
// // rLeftHandSideMatrix(18,10)+=(crLeftHandSideMatrix108*rB(0,10) + crLeftHandSideMatrix109*rB(1,10) + crLeftHandSideMatrix110*rB(2,10) + crLeftHandSideMatrix111*rB(3,10) + crLeftHandSideMatrix112*rB(4,10) + crLeftHandSideMatrix113*rB(5,10));
// // rLeftHandSideMatrix(18,11)+=(crLeftHandSideMatrix108*rB(0,11) + crLeftHandSideMatrix109*rB(1,11) + crLeftHandSideMatrix110*rB(2,11) + crLeftHandSideMatrix111*rB(3,11) + crLeftHandSideMatrix112*rB(4,11) + crLeftHandSideMatrix113*rB(5,11));
// // rLeftHandSideMatrix(18,12)+=(crLeftHandSideMatrix108*rB(0,12) + crLeftHandSideMatrix109*rB(1,12) + crLeftHandSideMatrix110*rB(2,12) + crLeftHandSideMatrix111*rB(3,12) + crLeftHandSideMatrix112*rB(4,12) + crLeftHandSideMatrix113*rB(5,12));
// // rLeftHandSideMatrix(18,13)+=(crLeftHandSideMatrix108*rB(0,13) + crLeftHandSideMatrix109*rB(1,13) + crLeftHandSideMatrix110*rB(2,13) + crLeftHandSideMatrix111*rB(3,13) + crLeftHandSideMatrix112*rB(4,13) + crLeftHandSideMatrix113*rB(5,13));
// // rLeftHandSideMatrix(18,14)+=(crLeftHandSideMatrix108*rB(0,14) + crLeftHandSideMatrix109*rB(1,14) + crLeftHandSideMatrix110*rB(2,14) + crLeftHandSideMatrix111*rB(3,14) + crLeftHandSideMatrix112*rB(4,14) + crLeftHandSideMatrix113*rB(5,14));
// // rLeftHandSideMatrix(18,15)+=(crLeftHandSideMatrix108*rB(0,15) + crLeftHandSideMatrix109*rB(1,15) + crLeftHandSideMatrix110*rB(2,15) + crLeftHandSideMatrix111*rB(3,15) + crLeftHandSideMatrix112*rB(4,15) + crLeftHandSideMatrix113*rB(5,15));
// // rLeftHandSideMatrix(18,16)+=(crLeftHandSideMatrix108*rB(0,16) + crLeftHandSideMatrix109*rB(1,16) + crLeftHandSideMatrix110*rB(2,16) + crLeftHandSideMatrix111*rB(3,16) + crLeftHandSideMatrix112*rB(4,16) + crLeftHandSideMatrix113*rB(5,16));
// // rLeftHandSideMatrix(18,17)+=(crLeftHandSideMatrix108*rB(0,17) + crLeftHandSideMatrix109*rB(1,17) + crLeftHandSideMatrix110*rB(2,17) + crLeftHandSideMatrix111*rB(3,17) + crLeftHandSideMatrix112*rB(4,17) + crLeftHandSideMatrix113*rB(5,17));
// rLeftHandSideMatrix(18,18)+=(crLeftHandSideMatrix108*rB(0,18) + crLeftHandSideMatrix109*rB(1,18) + crLeftHandSideMatrix110*rB(2,18) + crLeftHandSideMatrix111*rB(3,18) + crLeftHandSideMatrix112*rB(4,18) + crLeftHandSideMatrix113*rB(5,18));
// rLeftHandSideMatrix(18,19)+=(crLeftHandSideMatrix108*rB(0,19) + crLeftHandSideMatrix109*rB(1,19) + crLeftHandSideMatrix110*rB(2,19) + crLeftHandSideMatrix111*rB(3,19) + crLeftHandSideMatrix112*rB(4,19) + crLeftHandSideMatrix113*rB(5,19));
// rLeftHandSideMatrix(18,20)+=(crLeftHandSideMatrix108*rB(0,20) + crLeftHandSideMatrix109*rB(1,20) + crLeftHandSideMatrix110*rB(2,20) + crLeftHandSideMatrix111*rB(3,20) + crLeftHandSideMatrix112*rB(4,20) + crLeftHandSideMatrix113*rB(5,20));
// rLeftHandSideMatrix(18,21)+=(crLeftHandSideMatrix108*rB(0,21) + crLeftHandSideMatrix109*rB(1,21) + crLeftHandSideMatrix110*rB(2,21) + crLeftHandSideMatrix111*rB(3,21) + crLeftHandSideMatrix112*rB(4,21) + crLeftHandSideMatrix113*rB(5,21));
// rLeftHandSideMatrix(18,22)+=(crLeftHandSideMatrix108*rB(0,22) + crLeftHandSideMatrix109*rB(1,22) + crLeftHandSideMatrix110*rB(2,22) + crLeftHandSideMatrix111*rB(3,22) + crLeftHandSideMatrix112*rB(4,22) + crLeftHandSideMatrix113*rB(5,22));
// rLeftHandSideMatrix(18,23)+=(crLeftHandSideMatrix108*rB(0,23) + crLeftHandSideMatrix109*rB(1,23) + crLeftHandSideMatrix110*rB(2,23) + crLeftHandSideMatrix111*rB(3,23) + crLeftHandSideMatrix112*rB(4,23) + crLeftHandSideMatrix113*rB(5,23));

// // rLeftHandSideMatrix(19,0)+=(crLeftHandSideMatrix114*rB(0,0) + crLeftHandSideMatrix115*rB(1,0) + crLeftHandSideMatrix116*rB(2,0) + crLeftHandSideMatrix117*rB(3,0) + crLeftHandSideMatrix118*rB(4,0) + crLeftHandSideMatrix119*rB(5,0));
// // rLeftHandSideMatrix(19,1)+=(crLeftHandSideMatrix114*rB(0,1) + crLeftHandSideMatrix115*rB(1,1) + crLeftHandSideMatrix116*rB(2,1) + crLeftHandSideMatrix117*rB(3,1) + crLeftHandSideMatrix118*rB(4,1) + crLeftHandSideMatrix119*rB(5,1));
// // rLeftHandSideMatrix(19,2)+=(crLeftHandSideMatrix114*rB(0,2) + crLeftHandSideMatrix115*rB(1,2) + crLeftHandSideMatrix116*rB(2,2) + crLeftHandSideMatrix117*rB(3,2) + crLeftHandSideMatrix118*rB(4,2) + crLeftHandSideMatrix119*rB(5,2));
// // rLeftHandSideMatrix(19,3)+=(crLeftHandSideMatrix114*rB(0,3) + crLeftHandSideMatrix115*rB(1,3) + crLeftHandSideMatrix116*rB(2,3) + crLeftHandSideMatrix117*rB(3,3) + crLeftHandSideMatrix118*rB(4,3) + crLeftHandSideMatrix119*rB(5,3));
// // rLeftHandSideMatrix(19,4)+=(crLeftHandSideMatrix114*rB(0,4) + crLeftHandSideMatrix115*rB(1,4) + crLeftHandSideMatrix116*rB(2,4) + crLeftHandSideMatrix117*rB(3,4) + crLeftHandSideMatrix118*rB(4,4) + crLeftHandSideMatrix119*rB(5,4));
// // rLeftHandSideMatrix(19,5)+=(crLeftHandSideMatrix114*rB(0,5) + crLeftHandSideMatrix115*rB(1,5) + crLeftHandSideMatrix116*rB(2,5) + crLeftHandSideMatrix117*rB(3,5) + crLeftHandSideMatrix118*rB(4,5) + crLeftHandSideMatrix119*rB(5,5));
// // rLeftHandSideMatrix(19,6)+=(crLeftHandSideMatrix114*rB(0,6) + crLeftHandSideMatrix115*rB(1,6) + crLeftHandSideMatrix116*rB(2,6) + crLeftHandSideMatrix117*rB(3,6) + crLeftHandSideMatrix118*rB(4,6) + crLeftHandSideMatrix119*rB(5,6));
// // rLeftHandSideMatrix(19,7)+=(crLeftHandSideMatrix114*rB(0,7) + crLeftHandSideMatrix115*rB(1,7) + crLeftHandSideMatrix116*rB(2,7) + crLeftHandSideMatrix117*rB(3,7) + crLeftHandSideMatrix118*rB(4,7) + crLeftHandSideMatrix119*rB(5,7));
// // rLeftHandSideMatrix(19,8)+=(crLeftHandSideMatrix114*rB(0,8) + crLeftHandSideMatrix115*rB(1,8) + crLeftHandSideMatrix116*rB(2,8) + crLeftHandSideMatrix117*rB(3,8) + crLeftHandSideMatrix118*rB(4,8) + crLeftHandSideMatrix119*rB(5,8));
// // rLeftHandSideMatrix(19,9)+=(crLeftHandSideMatrix114*rB(0,9) + crLeftHandSideMatrix115*rB(1,9) + crLeftHandSideMatrix116*rB(2,9) + crLeftHandSideMatrix117*rB(3,9) + crLeftHandSideMatrix118*rB(4,9) + crLeftHandSideMatrix119*rB(5,9));
// // rLeftHandSideMatrix(19,10)+=(crLeftHandSideMatrix114*rB(0,10) + crLeftHandSideMatrix115*rB(1,10) + crLeftHandSideMatrix116*rB(2,10) + crLeftHandSideMatrix117*rB(3,10) + crLeftHandSideMatrix118*rB(4,10) + crLeftHandSideMatrix119*rB(5,10));
// // rLeftHandSideMatrix(19,11)+=(crLeftHandSideMatrix114*rB(0,11) + crLeftHandSideMatrix115*rB(1,11) + crLeftHandSideMatrix116*rB(2,11) + crLeftHandSideMatrix117*rB(3,11) + crLeftHandSideMatrix118*rB(4,11) + crLeftHandSideMatrix119*rB(5,11));
// // rLeftHandSideMatrix(19,12)+=(crLeftHandSideMatrix114*rB(0,12) + crLeftHandSideMatrix115*rB(1,12) + crLeftHandSideMatrix116*rB(2,12) + crLeftHandSideMatrix117*rB(3,12) + crLeftHandSideMatrix118*rB(4,12) + crLeftHandSideMatrix119*rB(5,12));
// // rLeftHandSideMatrix(19,13)+=(crLeftHandSideMatrix114*rB(0,13) + crLeftHandSideMatrix115*rB(1,13) + crLeftHandSideMatrix116*rB(2,13) + crLeftHandSideMatrix117*rB(3,13) + crLeftHandSideMatrix118*rB(4,13) + crLeftHandSideMatrix119*rB(5,13));
// // rLeftHandSideMatrix(19,14)+=(crLeftHandSideMatrix114*rB(0,14) + crLeftHandSideMatrix115*rB(1,14) + crLeftHandSideMatrix116*rB(2,14) + crLeftHandSideMatrix117*rB(3,14) + crLeftHandSideMatrix118*rB(4,14) + crLeftHandSideMatrix119*rB(5,14));
// // rLeftHandSideMatrix(19,15)+=(crLeftHandSideMatrix114*rB(0,15) + crLeftHandSideMatrix115*rB(1,15) + crLeftHandSideMatrix116*rB(2,15) + crLeftHandSideMatrix117*rB(3,15) + crLeftHandSideMatrix118*rB(4,15) + crLeftHandSideMatrix119*rB(5,15));
// // rLeftHandSideMatrix(19,16)+=(crLeftHandSideMatrix114*rB(0,16) + crLeftHandSideMatrix115*rB(1,16) + crLeftHandSideMatrix116*rB(2,16) + crLeftHandSideMatrix117*rB(3,16) + crLeftHandSideMatrix118*rB(4,16) + crLeftHandSideMatrix119*rB(5,16));
// // rLeftHandSideMatrix(19,17)+=(crLeftHandSideMatrix114*rB(0,17) + crLeftHandSideMatrix115*rB(1,17) + crLeftHandSideMatrix116*rB(2,17) + crLeftHandSideMatrix117*rB(3,17) + crLeftHandSideMatrix118*rB(4,17) + crLeftHandSideMatrix119*rB(5,17));
// // rLeftHandSideMatrix(19,18)+=(crLeftHandSideMatrix114*rB(0,18) + crLeftHandSideMatrix115*rB(1,18) + crLeftHandSideMatrix116*rB(2,18) + crLeftHandSideMatrix117*rB(3,18) + crLeftHandSideMatrix118*rB(4,18) + crLeftHandSideMatrix119*rB(5,18));
// rLeftHandSideMatrix(19,19)+=(crLeftHandSideMatrix114*rB(0,19) + crLeftHandSideMatrix115*rB(1,19) + crLeftHandSideMatrix116*rB(2,19) + crLeftHandSideMatrix117*rB(3,19) + crLeftHandSideMatrix118*rB(4,19) + crLeftHandSideMatrix119*rB(5,19));
// rLeftHandSideMatrix(19,20)+=(crLeftHandSideMatrix114*rB(0,20) + crLeftHandSideMatrix115*rB(1,20) + crLeftHandSideMatrix116*rB(2,20) + crLeftHandSideMatrix117*rB(3,20) + crLeftHandSideMatrix118*rB(4,20) + crLeftHandSideMatrix119*rB(5,20));
// rLeftHandSideMatrix(19,21)+=(crLeftHandSideMatrix114*rB(0,21) + crLeftHandSideMatrix115*rB(1,21) + crLeftHandSideMatrix116*rB(2,21) + crLeftHandSideMatrix117*rB(3,21) + crLeftHandSideMatrix118*rB(4,21) + crLeftHandSideMatrix119*rB(5,21));
// rLeftHandSideMatrix(19,22)+=(crLeftHandSideMatrix114*rB(0,22) + crLeftHandSideMatrix115*rB(1,22) + crLeftHandSideMatrix116*rB(2,22) + crLeftHandSideMatrix117*rB(3,22) + crLeftHandSideMatrix118*rB(4,22) + crLeftHandSideMatrix119*rB(5,22));
// rLeftHandSideMatrix(19,23)+=(crLeftHandSideMatrix114*rB(0,23) + crLeftHandSideMatrix115*rB(1,23) + crLeftHandSideMatrix116*rB(2,23) + crLeftHandSideMatrix117*rB(3,23) + crLeftHandSideMatrix118*rB(4,23) + crLeftHandSideMatrix119*rB(5,23));

// // rLeftHandSideMatrix(20,0)+=(crLeftHandSideMatrix120*rB(0,0) + crLeftHandSideMatrix121*rB(1,0) + crLeftHandSideMatrix122*rB(2,0) + crLeftHandSideMatrix123*rB(3,0) + crLeftHandSideMatrix124*rB(4,0) + crLeftHandSideMatrix125*rB(5,0));
// // rLeftHandSideMatrix(20,1)+=(crLeftHandSideMatrix120*rB(0,1) + crLeftHandSideMatrix121*rB(1,1) + crLeftHandSideMatrix122*rB(2,1) + crLeftHandSideMatrix123*rB(3,1) + crLeftHandSideMatrix124*rB(4,1) + crLeftHandSideMatrix125*rB(5,1));
// // rLeftHandSideMatrix(20,2)+=(crLeftHandSideMatrix120*rB(0,2) + crLeftHandSideMatrix121*rB(1,2) + crLeftHandSideMatrix122*rB(2,2) + crLeftHandSideMatrix123*rB(3,2) + crLeftHandSideMatrix124*rB(4,2) + crLeftHandSideMatrix125*rB(5,2));
// // rLeftHandSideMatrix(20,3)+=(crLeftHandSideMatrix120*rB(0,3) + crLeftHandSideMatrix121*rB(1,3) + crLeftHandSideMatrix122*rB(2,3) + crLeftHandSideMatrix123*rB(3,3) + crLeftHandSideMatrix124*rB(4,3) + crLeftHandSideMatrix125*rB(5,3));
// // rLeftHandSideMatrix(20,4)+=(crLeftHandSideMatrix120*rB(0,4) + crLeftHandSideMatrix121*rB(1,4) + crLeftHandSideMatrix122*rB(2,4) + crLeftHandSideMatrix123*rB(3,4) + crLeftHandSideMatrix124*rB(4,4) + crLeftHandSideMatrix125*rB(5,4));
// // rLeftHandSideMatrix(20,5)+=(crLeftHandSideMatrix120*rB(0,5) + crLeftHandSideMatrix121*rB(1,5) + crLeftHandSideMatrix122*rB(2,5) + crLeftHandSideMatrix123*rB(3,5) + crLeftHandSideMatrix124*rB(4,5) + crLeftHandSideMatrix125*rB(5,5));
// // rLeftHandSideMatrix(20,6)+=(crLeftHandSideMatrix120*rB(0,6) + crLeftHandSideMatrix121*rB(1,6) + crLeftHandSideMatrix122*rB(2,6) + crLeftHandSideMatrix123*rB(3,6) + crLeftHandSideMatrix124*rB(4,6) + crLeftHandSideMatrix125*rB(5,6));
// // rLeftHandSideMatrix(20,7)+=(crLeftHandSideMatrix120*rB(0,7) + crLeftHandSideMatrix121*rB(1,7) + crLeftHandSideMatrix122*rB(2,7) + crLeftHandSideMatrix123*rB(3,7) + crLeftHandSideMatrix124*rB(4,7) + crLeftHandSideMatrix125*rB(5,7));
// // rLeftHandSideMatrix(20,8)+=(crLeftHandSideMatrix120*rB(0,8) + crLeftHandSideMatrix121*rB(1,8) + crLeftHandSideMatrix122*rB(2,8) + crLeftHandSideMatrix123*rB(3,8) + crLeftHandSideMatrix124*rB(4,8) + crLeftHandSideMatrix125*rB(5,8));
// // rLeftHandSideMatrix(20,9)+=(crLeftHandSideMatrix120*rB(0,9) + crLeftHandSideMatrix121*rB(1,9) + crLeftHandSideMatrix122*rB(2,9) + crLeftHandSideMatrix123*rB(3,9) + crLeftHandSideMatrix124*rB(4,9) + crLeftHandSideMatrix125*rB(5,9));
// // rLeftHandSideMatrix(20,10)+=(crLeftHandSideMatrix120*rB(0,10) + crLeftHandSideMatrix121*rB(1,10) + crLeftHandSideMatrix122*rB(2,10) + crLeftHandSideMatrix123*rB(3,10) + crLeftHandSideMatrix124*rB(4,10) + crLeftHandSideMatrix125*rB(5,10));
// // rLeftHandSideMatrix(20,11)+=(crLeftHandSideMatrix120*rB(0,11) + crLeftHandSideMatrix121*rB(1,11) + crLeftHandSideMatrix122*rB(2,11) + crLeftHandSideMatrix123*rB(3,11) + crLeftHandSideMatrix124*rB(4,11) + crLeftHandSideMatrix125*rB(5,11));
// // rLeftHandSideMatrix(20,12)+=(crLeftHandSideMatrix120*rB(0,12) + crLeftHandSideMatrix121*rB(1,12) + crLeftHandSideMatrix122*rB(2,12) + crLeftHandSideMatrix123*rB(3,12) + crLeftHandSideMatrix124*rB(4,12) + crLeftHandSideMatrix125*rB(5,12));
// // rLeftHandSideMatrix(20,13)+=(crLeftHandSideMatrix120*rB(0,13) + crLeftHandSideMatrix121*rB(1,13) + crLeftHandSideMatrix122*rB(2,13) + crLeftHandSideMatrix123*rB(3,13) + crLeftHandSideMatrix124*rB(4,13) + crLeftHandSideMatrix125*rB(5,13));
// // rLeftHandSideMatrix(20,14)+=(crLeftHandSideMatrix120*rB(0,14) + crLeftHandSideMatrix121*rB(1,14) + crLeftHandSideMatrix122*rB(2,14) + crLeftHandSideMatrix123*rB(3,14) + crLeftHandSideMatrix124*rB(4,14) + crLeftHandSideMatrix125*rB(5,14));
// // rLeftHandSideMatrix(20,15)+=(crLeftHandSideMatrix120*rB(0,15) + crLeftHandSideMatrix121*rB(1,15) + crLeftHandSideMatrix122*rB(2,15) + crLeftHandSideMatrix123*rB(3,15) + crLeftHandSideMatrix124*rB(4,15) + crLeftHandSideMatrix125*rB(5,15));
// // rLeftHandSideMatrix(20,16)+=(crLeftHandSideMatrix120*rB(0,16) + crLeftHandSideMatrix121*rB(1,16) + crLeftHandSideMatrix122*rB(2,16) + crLeftHandSideMatrix123*rB(3,16) + crLeftHandSideMatrix124*rB(4,16) + crLeftHandSideMatrix125*rB(5,16));
// // rLeftHandSideMatrix(20,17)+=(crLeftHandSideMatrix120*rB(0,17) + crLeftHandSideMatrix121*rB(1,17) + crLeftHandSideMatrix122*rB(2,17) + crLeftHandSideMatrix123*rB(3,17) + crLeftHandSideMatrix124*rB(4,17) + crLeftHandSideMatrix125*rB(5,17));
// // rLeftHandSideMatrix(20,18)+=(crLeftHandSideMatrix120*rB(0,18) + crLeftHandSideMatrix121*rB(1,18) + crLeftHandSideMatrix122*rB(2,18) + crLeftHandSideMatrix123*rB(3,18) + crLeftHandSideMatrix124*rB(4,18) + crLeftHandSideMatrix125*rB(5,18));
// // rLeftHandSideMatrix(20,19)+=(crLeftHandSideMatrix120*rB(0,19) + crLeftHandSideMatrix121*rB(1,19) + crLeftHandSideMatrix122*rB(2,19) + crLeftHandSideMatrix123*rB(3,19) + crLeftHandSideMatrix124*rB(4,19) + crLeftHandSideMatrix125*rB(5,19));
// rLeftHandSideMatrix(20,20)+=(crLeftHandSideMatrix120*rB(0,20) + crLeftHandSideMatrix121*rB(1,20) + crLeftHandSideMatrix122*rB(2,20) + crLeftHandSideMatrix123*rB(3,20) + crLeftHandSideMatrix124*rB(4,20) + crLeftHandSideMatrix125*rB(5,20));
// rLeftHandSideMatrix(20,21)+=(crLeftHandSideMatrix120*rB(0,21) + crLeftHandSideMatrix121*rB(1,21) + crLeftHandSideMatrix122*rB(2,21) + crLeftHandSideMatrix123*rB(3,21) + crLeftHandSideMatrix124*rB(4,21) + crLeftHandSideMatrix125*rB(5,21));
// rLeftHandSideMatrix(20,22)+=(crLeftHandSideMatrix120*rB(0,22) + crLeftHandSideMatrix121*rB(1,22) + crLeftHandSideMatrix122*rB(2,22) + crLeftHandSideMatrix123*rB(3,22) + crLeftHandSideMatrix124*rB(4,22) + crLeftHandSideMatrix125*rB(5,22));
// rLeftHandSideMatrix(20,23)+=(crLeftHandSideMatrix120*rB(0,23) + crLeftHandSideMatrix121*rB(1,23) + crLeftHandSideMatrix122*rB(2,23) + crLeftHandSideMatrix123*rB(3,23) + crLeftHandSideMatrix124*rB(4,23) + crLeftHandSideMatrix125*rB(5,23));

// // rLeftHandSideMatrix(21,0)+=(crLeftHandSideMatrix126*rB(0,0) + crLeftHandSideMatrix127*rB(1,0) + crLeftHandSideMatrix128*rB(2,0) + crLeftHandSideMatrix129*rB(3,0) + crLeftHandSideMatrix130*rB(4,0) + crLeftHandSideMatrix131*rB(5,0));
// // rLeftHandSideMatrix(21,1)+=(crLeftHandSideMatrix126*rB(0,1) + crLeftHandSideMatrix127*rB(1,1) + crLeftHandSideMatrix128*rB(2,1) + crLeftHandSideMatrix129*rB(3,1) + crLeftHandSideMatrix130*rB(4,1) + crLeftHandSideMatrix131*rB(5,1));
// // rLeftHandSideMatrix(21,2)+=(crLeftHandSideMatrix126*rB(0,2) + crLeftHandSideMatrix127*rB(1,2) + crLeftHandSideMatrix128*rB(2,2) + crLeftHandSideMatrix129*rB(3,2) + crLeftHandSideMatrix130*rB(4,2) + crLeftHandSideMatrix131*rB(5,2));
// // rLeftHandSideMatrix(21,3)+=(crLeftHandSideMatrix126*rB(0,3) + crLeftHandSideMatrix127*rB(1,3) + crLeftHandSideMatrix128*rB(2,3) + crLeftHandSideMatrix129*rB(3,3) + crLeftHandSideMatrix130*rB(4,3) + crLeftHandSideMatrix131*rB(5,3));
// // rLeftHandSideMatrix(21,4)+=(crLeftHandSideMatrix126*rB(0,4) + crLeftHandSideMatrix127*rB(1,4) + crLeftHandSideMatrix128*rB(2,4) + crLeftHandSideMatrix129*rB(3,4) + crLeftHandSideMatrix130*rB(4,4) + crLeftHandSideMatrix131*rB(5,4));
// // rLeftHandSideMatrix(21,5)+=(crLeftHandSideMatrix126*rB(0,5) + crLeftHandSideMatrix127*rB(1,5) + crLeftHandSideMatrix128*rB(2,5) + crLeftHandSideMatrix129*rB(3,5) + crLeftHandSideMatrix130*rB(4,5) + crLeftHandSideMatrix131*rB(5,5));
// // rLeftHandSideMatrix(21,6)+=(crLeftHandSideMatrix126*rB(0,6) + crLeftHandSideMatrix127*rB(1,6) + crLeftHandSideMatrix128*rB(2,6) + crLeftHandSideMatrix129*rB(3,6) + crLeftHandSideMatrix130*rB(4,6) + crLeftHandSideMatrix131*rB(5,6));
// // rLeftHandSideMatrix(21,7)+=(crLeftHandSideMatrix126*rB(0,7) + crLeftHandSideMatrix127*rB(1,7) + crLeftHandSideMatrix128*rB(2,7) + crLeftHandSideMatrix129*rB(3,7) + crLeftHandSideMatrix130*rB(4,7) + crLeftHandSideMatrix131*rB(5,7));
// // rLeftHandSideMatrix(21,8)+=(crLeftHandSideMatrix126*rB(0,8) + crLeftHandSideMatrix127*rB(1,8) + crLeftHandSideMatrix128*rB(2,8) + crLeftHandSideMatrix129*rB(3,8) + crLeftHandSideMatrix130*rB(4,8) + crLeftHandSideMatrix131*rB(5,8));
// // rLeftHandSideMatrix(21,9)+=(crLeftHandSideMatrix126*rB(0,9) + crLeftHandSideMatrix127*rB(1,9) + crLeftHandSideMatrix128*rB(2,9) + crLeftHandSideMatrix129*rB(3,9) + crLeftHandSideMatrix130*rB(4,9) + crLeftHandSideMatrix131*rB(5,9));
// // rLeftHandSideMatrix(21,10)+=(crLeftHandSideMatrix126*rB(0,10) + crLeftHandSideMatrix127*rB(1,10) + crLeftHandSideMatrix128*rB(2,10) + crLeftHandSideMatrix129*rB(3,10) + crLeftHandSideMatrix130*rB(4,10) + crLeftHandSideMatrix131*rB(5,10));
// // rLeftHandSideMatrix(21,11)+=(crLeftHandSideMatrix126*rB(0,11) + crLeftHandSideMatrix127*rB(1,11) + crLeftHandSideMatrix128*rB(2,11) + crLeftHandSideMatrix129*rB(3,11) + crLeftHandSideMatrix130*rB(4,11) + crLeftHandSideMatrix131*rB(5,11));
// // rLeftHandSideMatrix(21,12)+=(crLeftHandSideMatrix126*rB(0,12) + crLeftHandSideMatrix127*rB(1,12) + crLeftHandSideMatrix128*rB(2,12) + crLeftHandSideMatrix129*rB(3,12) + crLeftHandSideMatrix130*rB(4,12) + crLeftHandSideMatrix131*rB(5,12));
// // rLeftHandSideMatrix(21,13)+=(crLeftHandSideMatrix126*rB(0,13) + crLeftHandSideMatrix127*rB(1,13) + crLeftHandSideMatrix128*rB(2,13) + crLeftHandSideMatrix129*rB(3,13) + crLeftHandSideMatrix130*rB(4,13) + crLeftHandSideMatrix131*rB(5,13));
// // rLeftHandSideMatrix(21,14)+=(crLeftHandSideMatrix126*rB(0,14) + crLeftHandSideMatrix127*rB(1,14) + crLeftHandSideMatrix128*rB(2,14) + crLeftHandSideMatrix129*rB(3,14) + crLeftHandSideMatrix130*rB(4,14) + crLeftHandSideMatrix131*rB(5,14));
// // rLeftHandSideMatrix(21,15)+=(crLeftHandSideMatrix126*rB(0,15) + crLeftHandSideMatrix127*rB(1,15) + crLeftHandSideMatrix128*rB(2,15) + crLeftHandSideMatrix129*rB(3,15) + crLeftHandSideMatrix130*rB(4,15) + crLeftHandSideMatrix131*rB(5,15));
// // rLeftHandSideMatrix(21,16)+=(crLeftHandSideMatrix126*rB(0,16) + crLeftHandSideMatrix127*rB(1,16) + crLeftHandSideMatrix128*rB(2,16) + crLeftHandSideMatrix129*rB(3,16) + crLeftHandSideMatrix130*rB(4,16) + crLeftHandSideMatrix131*rB(5,16));
// // rLeftHandSideMatrix(21,17)+=(crLeftHandSideMatrix126*rB(0,17) + crLeftHandSideMatrix127*rB(1,17) + crLeftHandSideMatrix128*rB(2,17) + crLeftHandSideMatrix129*rB(3,17) + crLeftHandSideMatrix130*rB(4,17) + crLeftHandSideMatrix131*rB(5,17));
// // rLeftHandSideMatrix(21,18)+=(crLeftHandSideMatrix126*rB(0,18) + crLeftHandSideMatrix127*rB(1,18) + crLeftHandSideMatrix128*rB(2,18) + crLeftHandSideMatrix129*rB(3,18) + crLeftHandSideMatrix130*rB(4,18) + crLeftHandSideMatrix131*rB(5,18));
// // rLeftHandSideMatrix(21,19)+=(crLeftHandSideMatrix126*rB(0,19) + crLeftHandSideMatrix127*rB(1,19) + crLeftHandSideMatrix128*rB(2,19) + crLeftHandSideMatrix129*rB(3,19) + crLeftHandSideMatrix130*rB(4,19) + crLeftHandSideMatrix131*rB(5,19));
// // rLeftHandSideMatrix(21,20)+=(crLeftHandSideMatrix126*rB(0,20) + crLeftHandSideMatrix127*rB(1,20) + crLeftHandSideMatrix128*rB(2,20) + crLeftHandSideMatrix129*rB(3,20) + crLeftHandSideMatrix130*rB(4,20) + crLeftHandSideMatrix131*rB(5,20));
// rLeftHandSideMatrix(21,21)+=(crLeftHandSideMatrix126*rB(0,21) + crLeftHandSideMatrix127*rB(1,21) + crLeftHandSideMatrix128*rB(2,21) + crLeftHandSideMatrix129*rB(3,21) + crLeftHandSideMatrix130*rB(4,21) + crLeftHandSideMatrix131*rB(5,21));
// rLeftHandSideMatrix(21,22)+=(crLeftHandSideMatrix126*rB(0,22) + crLeftHandSideMatrix127*rB(1,22) + crLeftHandSideMatrix128*rB(2,22) + crLeftHandSideMatrix129*rB(3,22) + crLeftHandSideMatrix130*rB(4,22) + crLeftHandSideMatrix131*rB(5,22));
// rLeftHandSideMatrix(21,23)+=(crLeftHandSideMatrix126*rB(0,23) + crLeftHandSideMatrix127*rB(1,23) + crLeftHandSideMatrix128*rB(2,23) + crLeftHandSideMatrix129*rB(3,23) + crLeftHandSideMatrix130*rB(4,23) + crLeftHandSideMatrix131*rB(5,23));

// // rLeftHandSideMatrix(22,0)+=(crLeftHandSideMatrix132*rB(0,0) + crLeftHandSideMatrix133*rB(1,0) + crLeftHandSideMatrix134*rB(2,0) + crLeftHandSideMatrix135*rB(3,0) + crLeftHandSideMatrix136*rB(4,0) + crLeftHandSideMatrix137*rB(5,0));
// // rLeftHandSideMatrix(22,1)+=(crLeftHandSideMatrix132*rB(0,1) + crLeftHandSideMatrix133*rB(1,1) + crLeftHandSideMatrix134*rB(2,1) + crLeftHandSideMatrix135*rB(3,1) + crLeftHandSideMatrix136*rB(4,1) + crLeftHandSideMatrix137*rB(5,1));
// // rLeftHandSideMatrix(22,2)+=(crLeftHandSideMatrix132*rB(0,2) + crLeftHandSideMatrix133*rB(1,2) + crLeftHandSideMatrix134*rB(2,2) + crLeftHandSideMatrix135*rB(3,2) + crLeftHandSideMatrix136*rB(4,2) + crLeftHandSideMatrix137*rB(5,2));
// // rLeftHandSideMatrix(22,3)+=(crLeftHandSideMatrix132*rB(0,3) + crLeftHandSideMatrix133*rB(1,3) + crLeftHandSideMatrix134*rB(2,3) + crLeftHandSideMatrix135*rB(3,3) + crLeftHandSideMatrix136*rB(4,3) + crLeftHandSideMatrix137*rB(5,3));
// // rLeftHandSideMatrix(22,4)+=(crLeftHandSideMatrix132*rB(0,4) + crLeftHandSideMatrix133*rB(1,4) + crLeftHandSideMatrix134*rB(2,4) + crLeftHandSideMatrix135*rB(3,4) + crLeftHandSideMatrix136*rB(4,4) + crLeftHandSideMatrix137*rB(5,4));
// // rLeftHandSideMatrix(22,5)+=(crLeftHandSideMatrix132*rB(0,5) + crLeftHandSideMatrix133*rB(1,5) + crLeftHandSideMatrix134*rB(2,5) + crLeftHandSideMatrix135*rB(3,5) + crLeftHandSideMatrix136*rB(4,5) + crLeftHandSideMatrix137*rB(5,5));
// // rLeftHandSideMatrix(22,6)+=(crLeftHandSideMatrix132*rB(0,6) + crLeftHandSideMatrix133*rB(1,6) + crLeftHandSideMatrix134*rB(2,6) + crLeftHandSideMatrix135*rB(3,6) + crLeftHandSideMatrix136*rB(4,6) + crLeftHandSideMatrix137*rB(5,6));
// // rLeftHandSideMatrix(22,7)+=(crLeftHandSideMatrix132*rB(0,7) + crLeftHandSideMatrix133*rB(1,7) + crLeftHandSideMatrix134*rB(2,7) + crLeftHandSideMatrix135*rB(3,7) + crLeftHandSideMatrix136*rB(4,7) + crLeftHandSideMatrix137*rB(5,7));
// // rLeftHandSideMatrix(22,8)+=(crLeftHandSideMatrix132*rB(0,8) + crLeftHandSideMatrix133*rB(1,8) + crLeftHandSideMatrix134*rB(2,8) + crLeftHandSideMatrix135*rB(3,8) + crLeftHandSideMatrix136*rB(4,8) + crLeftHandSideMatrix137*rB(5,8));
// // rLeftHandSideMatrix(22,9)+=(crLeftHandSideMatrix132*rB(0,9) + crLeftHandSideMatrix133*rB(1,9) + crLeftHandSideMatrix134*rB(2,9) + crLeftHandSideMatrix135*rB(3,9) + crLeftHandSideMatrix136*rB(4,9) + crLeftHandSideMatrix137*rB(5,9));
// // rLeftHandSideMatrix(22,10)+=(crLeftHandSideMatrix132*rB(0,10) + crLeftHandSideMatrix133*rB(1,10) + crLeftHandSideMatrix134*rB(2,10) + crLeftHandSideMatrix135*rB(3,10) + crLeftHandSideMatrix136*rB(4,10) + crLeftHandSideMatrix137*rB(5,10));
// // rLeftHandSideMatrix(22,11)+=(crLeftHandSideMatrix132*rB(0,11) + crLeftHandSideMatrix133*rB(1,11) + crLeftHandSideMatrix134*rB(2,11) + crLeftHandSideMatrix135*rB(3,11) + crLeftHandSideMatrix136*rB(4,11) + crLeftHandSideMatrix137*rB(5,11));
// // rLeftHandSideMatrix(22,12)+=(crLeftHandSideMatrix132*rB(0,12) + crLeftHandSideMatrix133*rB(1,12) + crLeftHandSideMatrix134*rB(2,12) + crLeftHandSideMatrix135*rB(3,12) + crLeftHandSideMatrix136*rB(4,12) + crLeftHandSideMatrix137*rB(5,12));
// // rLeftHandSideMatrix(22,13)+=(crLeftHandSideMatrix132*rB(0,13) + crLeftHandSideMatrix133*rB(1,13) + crLeftHandSideMatrix134*rB(2,13) + crLeftHandSideMatrix135*rB(3,13) + crLeftHandSideMatrix136*rB(4,13) + crLeftHandSideMatrix137*rB(5,13));
// // rLeftHandSideMatrix(22,14)+=(crLeftHandSideMatrix132*rB(0,14) + crLeftHandSideMatrix133*rB(1,14) + crLeftHandSideMatrix134*rB(2,14) + crLeftHandSideMatrix135*rB(3,14) + crLeftHandSideMatrix136*rB(4,14) + crLeftHandSideMatrix137*rB(5,14));
// // rLeftHandSideMatrix(22,15)+=(crLeftHandSideMatrix132*rB(0,15) + crLeftHandSideMatrix133*rB(1,15) + crLeftHandSideMatrix134*rB(2,15) + crLeftHandSideMatrix135*rB(3,15) + crLeftHandSideMatrix136*rB(4,15) + crLeftHandSideMatrix137*rB(5,15));
// // rLeftHandSideMatrix(22,16)+=(crLeftHandSideMatrix132*rB(0,16) + crLeftHandSideMatrix133*rB(1,16) + crLeftHandSideMatrix134*rB(2,16) + crLeftHandSideMatrix135*rB(3,16) + crLeftHandSideMatrix136*rB(4,16) + crLeftHandSideMatrix137*rB(5,16));
// // rLeftHandSideMatrix(22,17)+=(crLeftHandSideMatrix132*rB(0,17) + crLeftHandSideMatrix133*rB(1,17) + crLeftHandSideMatrix134*rB(2,17) + crLeftHandSideMatrix135*rB(3,17) + crLeftHandSideMatrix136*rB(4,17) + crLeftHandSideMatrix137*rB(5,17));
// // rLeftHandSideMatrix(22,18)+=(crLeftHandSideMatrix132*rB(0,18) + crLeftHandSideMatrix133*rB(1,18) + crLeftHandSideMatrix134*rB(2,18) + crLeftHandSideMatrix135*rB(3,18) + crLeftHandSideMatrix136*rB(4,18) + crLeftHandSideMatrix137*rB(5,18));
// // rLeftHandSideMatrix(22,19)+=(crLeftHandSideMatrix132*rB(0,19) + crLeftHandSideMatrix133*rB(1,19) + crLeftHandSideMatrix134*rB(2,19) + crLeftHandSideMatrix135*rB(3,19) + crLeftHandSideMatrix136*rB(4,19) + crLeftHandSideMatrix137*rB(5,19));
// // rLeftHandSideMatrix(22,20)+=(crLeftHandSideMatrix132*rB(0,20) + crLeftHandSideMatrix133*rB(1,20) + crLeftHandSideMatrix134*rB(2,20) + crLeftHandSideMatrix135*rB(3,20) + crLeftHandSideMatrix136*rB(4,20) + crLeftHandSideMatrix137*rB(5,20));
// // rLeftHandSideMatrix(22,21)+=(crLeftHandSideMatrix132*rB(0,21) + crLeftHandSideMatrix133*rB(1,21) + crLeftHandSideMatrix134*rB(2,21) + crLeftHandSideMatrix135*rB(3,21) + crLeftHandSideMatrix136*rB(4,21) + crLeftHandSideMatrix137*rB(5,21));
// rLeftHandSideMatrix(22,22)+=(crLeftHandSideMatrix132*rB(0,22) + crLeftHandSideMatrix133*rB(1,22) + crLeftHandSideMatrix134*rB(2,22) + crLeftHandSideMatrix135*rB(3,22) + crLeftHandSideMatrix136*rB(4,22) + crLeftHandSideMatrix137*rB(5,22));
// rLeftHandSideMatrix(22,23)+=(crLeftHandSideMatrix132*rB(0,23) + crLeftHandSideMatrix133*rB(1,23) + crLeftHandSideMatrix134*rB(2,23) + crLeftHandSideMatrix135*rB(3,23) + crLeftHandSideMatrix136*rB(4,23) + crLeftHandSideMatrix137*rB(5,23));

// // rLeftHandSideMatrix(23,0)+=(crLeftHandSideMatrix138*rB(0,0) + crLeftHandSideMatrix139*rB(1,0) + crLeftHandSideMatrix140*rB(2,0) + crLeftHandSideMatrix141*rB(3,0) + crLeftHandSideMatrix142*rB(4,0) + crLeftHandSideMatrix143*rB(5,0));
// // rLeftHandSideMatrix(23,1)+=(crLeftHandSideMatrix138*rB(0,1) + crLeftHandSideMatrix139*rB(1,1) + crLeftHandSideMatrix140*rB(2,1) + crLeftHandSideMatrix141*rB(3,1) + crLeftHandSideMatrix142*rB(4,1) + crLeftHandSideMatrix143*rB(5,1));
// // rLeftHandSideMatrix(23,2)+=(crLeftHandSideMatrix138*rB(0,2) + crLeftHandSideMatrix139*rB(1,2) + crLeftHandSideMatrix140*rB(2,2) + crLeftHandSideMatrix141*rB(3,2) + crLeftHandSideMatrix142*rB(4,2) + crLeftHandSideMatrix143*rB(5,2));
// // rLeftHandSideMatrix(23,3)+=(crLeftHandSideMatrix138*rB(0,3) + crLeftHandSideMatrix139*rB(1,3) + crLeftHandSideMatrix140*rB(2,3) + crLeftHandSideMatrix141*rB(3,3) + crLeftHandSideMatrix142*rB(4,3) + crLeftHandSideMatrix143*rB(5,3));
// // rLeftHandSideMatrix(23,4)+=(crLeftHandSideMatrix138*rB(0,4) + crLeftHandSideMatrix139*rB(1,4) + crLeftHandSideMatrix140*rB(2,4) + crLeftHandSideMatrix141*rB(3,4) + crLeftHandSideMatrix142*rB(4,4) + crLeftHandSideMatrix143*rB(5,4));
// // rLeftHandSideMatrix(23,5)+=(crLeftHandSideMatrix138*rB(0,5) + crLeftHandSideMatrix139*rB(1,5) + crLeftHandSideMatrix140*rB(2,5) + crLeftHandSideMatrix141*rB(3,5) + crLeftHandSideMatrix142*rB(4,5) + crLeftHandSideMatrix143*rB(5,5));
// // rLeftHandSideMatrix(23,6)+=(crLeftHandSideMatrix138*rB(0,6) + crLeftHandSideMatrix139*rB(1,6) + crLeftHandSideMatrix140*rB(2,6) + crLeftHandSideMatrix141*rB(3,6) + crLeftHandSideMatrix142*rB(4,6) + crLeftHandSideMatrix143*rB(5,6));
// // rLeftHandSideMatrix(23,7)+=(crLeftHandSideMatrix138*rB(0,7) + crLeftHandSideMatrix139*rB(1,7) + crLeftHandSideMatrix140*rB(2,7) + crLeftHandSideMatrix141*rB(3,7) + crLeftHandSideMatrix142*rB(4,7) + crLeftHandSideMatrix143*rB(5,7));
// // rLeftHandSideMatrix(23,8)+=(crLeftHandSideMatrix138*rB(0,8) + crLeftHandSideMatrix139*rB(1,8) + crLeftHandSideMatrix140*rB(2,8) + crLeftHandSideMatrix141*rB(3,8) + crLeftHandSideMatrix142*rB(4,8) + crLeftHandSideMatrix143*rB(5,8));
// // rLeftHandSideMatrix(23,9)+=(crLeftHandSideMatrix138*rB(0,9) + crLeftHandSideMatrix139*rB(1,9) + crLeftHandSideMatrix140*rB(2,9) + crLeftHandSideMatrix141*rB(3,9) + crLeftHandSideMatrix142*rB(4,9) + crLeftHandSideMatrix143*rB(5,9));
// // rLeftHandSideMatrix(23,10)+=(crLeftHandSideMatrix138*rB(0,10) + crLeftHandSideMatrix139*rB(1,10) + crLeftHandSideMatrix140*rB(2,10) + crLeftHandSideMatrix141*rB(3,10) + crLeftHandSideMatrix142*rB(4,10) + crLeftHandSideMatrix143*rB(5,10));
// // rLeftHandSideMatrix(23,11)+=(crLeftHandSideMatrix138*rB(0,11) + crLeftHandSideMatrix139*rB(1,11) + crLeftHandSideMatrix140*rB(2,11) + crLeftHandSideMatrix141*rB(3,11) + crLeftHandSideMatrix142*rB(4,11) + crLeftHandSideMatrix143*rB(5,11));
// // rLeftHandSideMatrix(23,12)+=(crLeftHandSideMatrix138*rB(0,12) + crLeftHandSideMatrix139*rB(1,12) + crLeftHandSideMatrix140*rB(2,12) + crLeftHandSideMatrix141*rB(3,12) + crLeftHandSideMatrix142*rB(4,12) + crLeftHandSideMatrix143*rB(5,12));
// // rLeftHandSideMatrix(23,13)+=(crLeftHandSideMatrix138*rB(0,13) + crLeftHandSideMatrix139*rB(1,13) + crLeftHandSideMatrix140*rB(2,13) + crLeftHandSideMatrix141*rB(3,13) + crLeftHandSideMatrix142*rB(4,13) + crLeftHandSideMatrix143*rB(5,13));
// // rLeftHandSideMatrix(23,14)+=(crLeftHandSideMatrix138*rB(0,14) + crLeftHandSideMatrix139*rB(1,14) + crLeftHandSideMatrix140*rB(2,14) + crLeftHandSideMatrix141*rB(3,14) + crLeftHandSideMatrix142*rB(4,14) + crLeftHandSideMatrix143*rB(5,14));
// // rLeftHandSideMatrix(23,15)+=(crLeftHandSideMatrix138*rB(0,15) + crLeftHandSideMatrix139*rB(1,15) + crLeftHandSideMatrix140*rB(2,15) + crLeftHandSideMatrix141*rB(3,15) + crLeftHandSideMatrix142*rB(4,15) + crLeftHandSideMatrix143*rB(5,15));
// // rLeftHandSideMatrix(23,16)+=(crLeftHandSideMatrix138*rB(0,16) + crLeftHandSideMatrix139*rB(1,16) + crLeftHandSideMatrix140*rB(2,16) + crLeftHandSideMatrix141*rB(3,16) + crLeftHandSideMatrix142*rB(4,16) + crLeftHandSideMatrix143*rB(5,16));
// // rLeftHandSideMatrix(23,17)+=(crLeftHandSideMatrix138*rB(0,17) + crLeftHandSideMatrix139*rB(1,17) + crLeftHandSideMatrix140*rB(2,17) + crLeftHandSideMatrix141*rB(3,17) + crLeftHandSideMatrix142*rB(4,17) + crLeftHandSideMatrix143*rB(5,17));
// // rLeftHandSideMatrix(23,18)+=(crLeftHandSideMatrix138*rB(0,18) + crLeftHandSideMatrix139*rB(1,18) + crLeftHandSideMatrix140*rB(2,18) + crLeftHandSideMatrix141*rB(3,18) + crLeftHandSideMatrix142*rB(4,18) + crLeftHandSideMatrix143*rB(5,18));
// // rLeftHandSideMatrix(23,19)+=(crLeftHandSideMatrix138*rB(0,19) + crLeftHandSideMatrix139*rB(1,19) + crLeftHandSideMatrix140*rB(2,19) + crLeftHandSideMatrix141*rB(3,19) + crLeftHandSideMatrix142*rB(4,19) + crLeftHandSideMatrix143*rB(5,19));
// // rLeftHandSideMatrix(23,20)+=(crLeftHandSideMatrix138*rB(0,20) + crLeftHandSideMatrix139*rB(1,20) + crLeftHandSideMatrix140*rB(2,20) + crLeftHandSideMatrix141*rB(3,20) + crLeftHandSideMatrix142*rB(4,20) + crLeftHandSideMatrix143*rB(5,20));
// // rLeftHandSideMatrix(23,21)+=(crLeftHandSideMatrix138*rB(0,21) + crLeftHandSideMatrix139*rB(1,21) + crLeftHandSideMatrix140*rB(2,21) + crLeftHandSideMatrix141*rB(3,21) + crLeftHandSideMatrix142*rB(4,21) + crLeftHandSideMatrix143*rB(5,21));
// // rLeftHandSideMatrix(23,22)+=(crLeftHandSideMatrix138*rB(0,22) + crLeftHandSideMatrix139*rB(1,22) + crLeftHandSideMatrix140*rB(2,22) + crLeftHandSideMatrix141*rB(3,22) + crLeftHandSideMatrix142*rB(4,22) + crLeftHandSideMatrix143*rB(5,22));
// rLeftHandSideMatrix(23,23)+=(crLeftHandSideMatrix138*rB(0,23) + crLeftHandSideMatrix139*rB(1,23) + crLeftHandSideMatrix140*rB(2,23) + crLeftHandSideMatrix141*rB(3,23) + crLeftHandSideMatrix142*rB(4,23) + crLeftHandSideMatrix143*rB(5,23));

// for (int i = 0; i < 24; ++i)
//     for (int j = 0; j < i; ++j)
//             rLeftHandSideMatrix(i, j) += rLeftHandSideMatrix(j, i);

        KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddKg(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& DN_DX,
    const ConstitutiveLaw::StressVectorType& StressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix stress_tensor_x_weigth = IntegrationWeight * MathUtils<double>::StressVectorToTensor( StressVector );
    Matrix reduced_Kg(DN_DX.size1(), DN_DX.size1());
    MathUtils<double>::BDBtProductOperation(reduced_Kg, stress_tensor_x_weigth, DN_DX);
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, reduced_Kg, dimension );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const ConstitutiveLaw::StressVectorType& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    this->CalculateAndAddExtForceContribution( rThisKinematicVariables.N, rCurrentProcessInfo, rBodyForce, rRightHandSideVector, IntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    BoundedMatrix<double, 24, 6> aux2;
    noalias(aux2) = trans(rThisKinematicVariables.B);
    noalias( rRightHandSideVector ) -= prod( aux2, IntegrationWeight * rStressVector );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateAndAddExtForceContribution(
    const Vector& rN,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    VectorType& rRightHandSideVector,
    const double Weight
    ) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const SizeType index = dimension * i;

        for ( IndexType j = 0; j < dimension; ++j )
            rRightHandSideVector[index + j] += Weight * rN[i] * rBodyForce[j];
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    if(!UseGeometryIntegrationMethod())
        KRATOS_ERROR << "CalculateLumpedMassVector not implemented for element-based integration in base class" << std::endl;


    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rLumpedMassVector.size() != mat_size)
        rLumpedMassVector.resize( mat_size, false );

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = GetGeometry().DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < dimension; ++j ) {
            IndexType index = i * dimension + j;
            rLumpedMassVector[index] = temp;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::CalculateDampingMatrixWithLumpedMass(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // 1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0.0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) )
        alpha = GetProperties()[RAYLEIGH_ALPHA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) )
        alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];

    double beta  = 0.0;
    if( GetProperties().Has(RAYLEIGH_BETA) )
        beta = GetProperties()[RAYLEIGH_BETA];
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) )
        beta = rCurrentProcessInfo[RAYLEIGH_BETA];

    // Compose the Damping Matrix:
    // Rayleigh Damping Matrix: alpha*M + beta*K

    // 2.-Calculate mass matrix:
    if (alpha > std::numeric_limits<double>::epsilon()) {
        VectorType temp_vector(mat_size);
        this->CalculateLumpedMassVector(temp_vector, rCurrentProcessInfo);
        for (IndexType i = 0; i < mat_size; ++i)
            rDampingMatrix(i, i) += alpha * temp_vector[i];
    }

    // 3.-Calculate StiffnessMatrix:
    if (beta > std::numeric_limits<double>::epsilon()) {
        MatrixType stiffness_matrix( mat_size, mat_size );
        VectorType residual_vector( mat_size );

        this->CalculateAll(stiffness_matrix, residual_vector, rCurrentProcessInfo, true, false);

        noalias( rDampingMatrix ) += beta  * stiffness_matrix;
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters BaseSolidElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static","implicit","explicit"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["INTEGRATION_WEIGHT","STRAIN_ENERGY","ERROR_INTEGRATION_POINT","VON_MISES_STRESS","INSITU_STRESS","CAUCHY_STRESS_VECTOR","PK2_STRESS_VECTOR","GREEN_LAGRANGE_STRAIN_VECTOR","ALMANSI_STRAIN_VECTOR","CAUCHY_STRESS_TENSOR","PK2_STRESS_TENSOR","GREEN_LAGRANGE_STRAIN_TENSOR","ALMANSI_STRAIN_TENSOR","CONSTITUTIVE_MATRIX","DEFORMATION_GRADIENT","CONSTITUTIVE_LAW"],
            "nodal_historical"       : ["DISPLACEMENT","VELOCITY","ACCELERATION"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Triangle2D6", "Quadrilateral2D4", "Quadrilateral2D8", "Quadrilateral2D9","Tetrahedra3D4", "Prism3D6", "Prism3D15", "Hexahedra3D8", "Hexahedra3D20", "Hexahedra3D27", "Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : -1,
        "documentation"   : "This is a pure displacement element"
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

bool BaseSolidElement::IsElementRotated() const
{
    if (mConstitutiveLawVector[0]->GetStrainSize() == 6) {
        return (this->Has(LOCAL_AXIS_1) && this->Has(LOCAL_AXIS_2));
    } else if (mConstitutiveLawVector[0]->GetStrainSize() == 3) {
        return (this->Has(LOCAL_AXIS_1));
    }
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}
} // Namespace Kratos
