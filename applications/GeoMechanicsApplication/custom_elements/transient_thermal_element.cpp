// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

// Application includes
#include "custom_elements/transient_thermal_element.hpp"

namespace Kratos
{
    // ============================================================================================
	// ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    Element::Pointer TransientThermalElement<TDim, TNumNodes>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TransientThermalElement(NewId, pGeom, pProperties));
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY
    	//
    	const GeometryType& rGeom = this->GetGeometry();
        const unsigned int N_DOF = this->GetNumberOfDOF();
        unsigned int index = 0;
        //
        if (rElementalDofList.size() != N_DOF)
            rElementalDofList.resize(N_DOF);
        //
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            rElementalDofList[index++] = rGeom[i].pGetDof(TEMPERATURE);
        }
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY
    	//
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int N_DOF = this->GetNumberOfDOF();
        unsigned int index = 0;
        //
        if (rResult.size() != N_DOF)
            rResult.resize(N_DOF, false);
        //
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            rResult[index++] = rGeom[i].GetDof(TEMPERATURE).EquationId();
        }
        //
        KRATOS_CATCH("")
    }

	// ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        //Resizing mass matrix
        if (rMassMatrix.size1() != N_DOF)
            rMassMatrix.resize(N_DOF, N_DOF, false);
        noalias(rMassMatrix) = ZeroMatrix(N_DOF, N_DOF);
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        //
    	// Rayleigh Method (Damping Matrix = alpha*M + beta*K)
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        // Compute Damping Matrix
        if (rDampingMatrix.size1() != N_DOF)
            rDampingMatrix.resize(N_DOF, N_DOF, false);
        noalias(rDampingMatrix) = ZeroMatrix(N_DOF, N_DOF);
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::GetValuesVector(
        Vector& rValues,
        int Step) const
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        if (rValues.size() != N_DOF)
            rValues.resize(N_DOF, false);
        //
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[i] = 0.0;
        }
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::GetFirstDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        if (rValues.size() != N_DOF)
            rValues.resize(N_DOF, false);
        //
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[i] = 0.0;
        }
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::GetSecondDerivativesVector(
        Vector& rValues,
        int Step) const
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        if (rValues.size() != N_DOF)
            rValues.resize(N_DOF, false);
        //
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rValues[i] = 0.0;
        }
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::Initialize(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::Initialize()") << this->Id() << std::endl;
        //
        const GeometryType& rGeom = this->GetGeometry();
        const unsigned int NumGPoints = rGeom.IntegrationPointsNumber(this->GetIntegrationMethod());
        //
        // pointer to constitutive laws
        if (mConstitutiveLawVector.size() != NumGPoints)
            mConstitutiveLawVector.resize(NumGPoints);
        //
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i) {
            mConstitutiveLawVector[i] = nullptr;
        }
        //
        mIsInitialised = true;
        //
        // KRATOS_INFO("1-TransientThermalElement::Initialize()") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    int TransientThermalElement<TDim, TNumNodes>::Check(
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::Check()") << this->Id() << std::endl;
        //
    	const PropertiesType& Prop = this->GetProperties();
        const GeometryType& Geom = this->GetGeometry();
        //
        if (Geom.DomainSize() < 1.0e-15)
            KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;
        //
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            if (Geom[i].SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_ERROR << "missing variable TEMPERATURE on node " << Geom[i].Id() << std::endl;
            //
            if (Geom[i].SolutionStepsDataHas(DT_TEMPERATURE) == false)
                KRATOS_ERROR << "missing variable DT_TEMPERATURE on node " << Geom[i].Id() << std::endl;
            //
            if (Geom[i].HasDofFor(TEMPERATURE) == false)
                KRATOS_ERROR << "missing variable TEMPERATURE on node " << Geom[i].Id() << std::endl;
        }
        //
        // Verify properties
        if (Prop.Has(DENSITY_WATER) == false || Prop[DENSITY_WATER] < 0.0)
            KRATOS_ERROR << "DENSITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(POROSITY) == false || Prop[POROSITY] < 0.0 || Prop[POROSITY] > 1.0)
            KRATOS_ERROR << "POROSITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(SATURATION) == false || Prop[SATURATION] < 0.0 || Prop[SATURATION] > 1.0)
            KRATOS_ERROR << "SATURATION does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(DENSITY_SOLID) == false || Prop[DENSITY_SOLID] < 0.0)
            KRATOS_ERROR << "DENSITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(HEAT_CAPACITY_WATER) == false || Prop[HEAT_CAPACITY_WATER] < 0.0)
            KRATOS_ERROR << "HEAT_CAPACITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(HEAT_CAPACITY_SOLID) == false || Prop[HEAT_CAPACITY_SOLID] < 0.0)
            KRATOS_ERROR << "HEAT_CAPACITY_SOLID does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_WATER) == false || Prop[THERMAL_CONDUCTIVITY_WATER] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_WATER does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_XX) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_XX] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_YY) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_YY] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_YY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_XY) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_XY] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_XY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(THERMAL_CONDUCTIVITY_SOLID_YX) == false || Prop[THERMAL_CONDUCTIVITY_SOLID_YX] < 0.0)
            KRATOS_ERROR << "THERMAL_CONDUCTIVITY_SOLID_YX does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(LONGITUDINAL_DISPERSIVITY) == false || Prop[LONGITUDINAL_DISPERSIVITY] < 0.0)
            KRATOS_ERROR << "LONGITUDINAL_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(TRANSVERSE_DISPERSIVITY) == false || Prop[TRANSVERSE_DISPERSIVITY] < 0.0)
            KRATOS_ERROR << "TRANSVERSE_DISPERSIVITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (Prop.Has(SOLID_COMPRESSIBILITY) == false || Prop[SOLID_COMPRESSIBILITY] < 0.0)
            KRATOS_ERROR << "SOLID_COMPRESSIBILITY does not exist in the material properties or has an invalid value at element" << this->Id() << std::endl;
        //
        if (TDim == 2) {
            // If this is a 2D problem, nodes must be in XY plane
            for (unsigned int i = 0; i < TNumNodes; ++i) {
                if (Geom[i].Z() != 0.0)
                    KRATOS_ERROR << " Node with non-zero Z coordinate found. Id: " << Geom[i].Id() << std::endl;
            }
        }
        //
        // KRATOS_INFO("1-TransientThermalElement::Check()") << this->Id() << std::endl;
        KRATOS_CATCH("");
        //
        return 0;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::InitializeSolutionStep()") << this->Id() << std::endl;
        //
        if (!mIsInitialised) this->Initialize(rCurrentProcessInfo);
        //
        // KRATOS_INFO("1-TransientThermalElement::InitializeSolutionStep()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::InitializeNonLinearIteration(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        // nothing
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::FinalizeNonLinearIteration(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        // nothing
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::FinalizeSolutionStep(
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        // nothing
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::CalculateAll()") << std::endl;
        //
    	//Previous definitions
        const GeometryType& Geom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(this->GetIntegrationMethod());
        const unsigned int NumGPoints = IntegrationPoints.size();
        //
        //Element variables
        ElementVariables Variables;
        this->InitializeElementVariables(Variables, rCurrentProcessInfo);
        //
        //Loop over integration points
        for (unsigned int GPoint = 0; GPoint < NumGPoints; ++GPoint) {
            //Compute GradNpT, B and StrainVector
            this->CalculateKinematics(Variables, GPoint);
            //
            //Compute weighting coefficient for integration
            Variables.IntegrationCoefficient =
                this->CalculateIntegrationCoefficient(IntegrationPoints, GPoint, Variables.detJ);
            //
            //Contributions to the left hand side
            if (CalculateStiffnessMatrixFlag) this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);
            //
            //Contributions to the right hand side
            if (CalculateResidualVectorFlag) this->CalculateAndAddRHS(rRightHandSideVector, Variables, GPoint);
        }
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateAll()") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeElementVariables(
        ElementVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::InitializeElementVariables()") << std::endl;
        //
    	//Properties variables
    	this->InitializeProperties(rVariables);
        //
    	//ProcessInfo variables
        rVariables.DtTemperatureCoefficient = rCurrentProcessInfo[DT_TEMPERATURE_COEFFICIENT];
        //
        //Nodal Variables
        this->InitializeNodalTemperatureVariables(rVariables);
        //
        //Variables computed at each GP
        rVariables.Np.resize(TNumNodes, false);
        rVariables.GradNpT.resize(TNumNodes, TDim, false);
        //
        //General Variables
        const GeometryType& Geom = this->GetGeometry();
        const unsigned int NumGPoints = Geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        //
        // shape functions
        (rVariables.NContainer).resize(NumGPoints, TNumNodes, false);
        rVariables.NContainer = Geom.ShapeFunctionsValues(this->GetIntegrationMethod());
        //
        // gradient of shape functions and determinant of Jacobian
        rVariables.detJContainer.resize(NumGPoints, false);
        //
        Geom.ShapeFunctionsIntegrationPointsGradients(rVariables.DN_DXContainer,
            rVariables.detJContainer,
            this->GetIntegrationMethod());
        //
        // KRATOS_INFO("1-TransientThermalElement::InitializeElementVariables()") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddLHS(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateAndAddLHS()") << std::endl;
        this->CalculateAndAddConductivityMatrix(rLeftHandSideMatrix, rVariables);
        //
        this->CalculateAndAddCapacityMatrix(rLeftHandSideMatrix, rVariables);
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateAndAddLHS()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateAndAddConductivityMatrix()") << std::endl;
        //
        this->CalculateConductivityMatrix(rVariables.TMatrix, rVariables);
        //
        //Distribute compressibility block matrix into the elemental matrix
        GeoElementUtilities::
            AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rVariables.TMatrix);
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateAndAddConductivityMatrix()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityMatrix(
        MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateAndAddCapacityMatrix()") << std::endl;
        //
        this->CalculateCapacityMatrix(rVariables.TMatrix, rVariables);
        //
        //Distribute permeability block matrix into the elemental matrix
        GeoElementUtilities::AssemblePBlockMatrix<0, TNumNodes>(rLeftHandSideMatrix, rVariables.TMatrix);
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateAndAddCapacityMatrix()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddRHS(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables,
        unsigned int GPoint)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateAndAddRHS()") << std::endl;
        //
        this->CalculateAndAddCapacityVector(rRightHandSideVector, rVariables);
        //
        this->CalculateAndAddConductivityVector(rRightHandSideVector, rVariables);
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateAndAddRHS()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateKinematics(
        ElementVariables& rVariables,
        unsigned int PointNumber)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::CalculateKinematics") << std::endl;
        //
    	//Setting the vector of shape functions and the matrix of the shape functions global gradients
    	rVariables.Np = row(rVariables.NContainer, PointNumber);
        rVariables.GradNpT = rVariables.DN_DXContainer[PointNumber];
        //
        rVariables.detJ = rVariables.detJContainer[PointNumber];
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateKinematics") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    unsigned int TransientThermalElement<TDim, TNumNodes>::GetNumberOfDOF() const
    {
        return TNumNodes;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeNodalTemperatureVariables(
        ElementVariables& rVariables)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::InitializeNodalTemperatureVariables") << std::endl;
        //
    	const GeometryType& rGeom = this->GetGeometry();
        //
        //Nodal Variables
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rVariables.TemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(TEMPERATURE);
            rVariables.DtTemperatureVector[i] = rGeom[i].FastGetSolutionStepValue(DT_TEMPERATURE);
        }
        //
        // KRATOS_INFO("1-TransientThermalElement::InitializeNodalTemperatureVariables") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    double TransientThermalElement<TDim, TNumNodes>::CalculateIntegrationCoefficient(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        unsigned int PointNumber,
        double detJ)
    {
        return IntegrationPoints[PointNumber].Weight() * detJ;
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityMatrix(
        BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
        ElementVariables& rVariables) const
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateCapacityMatrix()") << std::endl;
        //
        const double c1 = rVariables.Porosity * rVariables.Saturation
    	                * rVariables.WaterDensity * rVariables.WaterHeatCapacity;
        const double c2 = (1.0 - rVariables.Porosity) * rVariables.SolidDensity
    	                * rVariables.SolidHeatCapacity;
        TMatrix = (c1 + c2) * outer_prod(rVariables.Np, rVariables.Np)
    	                * rVariables.IntegrationCoefficient
    	                * rVariables.DtTemperatureCoefficient;
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateCapacityMatrix()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityMatrix(
        BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateConductivityMatrix()") << std::endl;
        //
    	this->CalculateThermalDispersionMatrix(rVariables.ConstitutiveMatrix, rVariables);
        //
        BoundedMatrix<double, TDim, TNumNodes> Temp = prod(rVariables.ConstitutiveMatrix, trans(rVariables.GradNpT));
        TMatrix = prod(rVariables.GradNpT, Temp) * rVariables.IntegrationCoefficient;
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateConductivityMatrix()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddCapacityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateAndAddCapacityVector()") << std::endl;
        //
        this->CalculateCapacityVector(rVariables.TMatrix, rVariables.TVector, rVariables);
        //
        //Distribute permeability block vector into elemental vector
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.TVector);

        // KRATOS_INFO("1-TransientThermalElement::CalculateAndAddCapacityVector()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateCapacityVector(
        BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
        array_1d<double, TNumNodes>& TVector,
        ElementVariables& rVariables) const
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientThermalElement::CalculateCapacityVector()") << std::endl;
        //
        const double c1 = rVariables.Porosity * rVariables.Saturation
    	                * rVariables.WaterDensity * rVariables.WaterHeatCapacity;
        const double c2 = (1.0 - rVariables.Porosity) * rVariables.SolidDensity
    	                * rVariables.SolidHeatCapacity;
        TMatrix = (c1 + c2) * outer_prod(rVariables.Np, rVariables.Np)
    	                * rVariables.IntegrationCoefficient;
        //
        TVector = - prod(TMatrix, rVariables.DtTemperatureVector);
        //
        // KRATOS_INFO("1-TransientThermalElement::CalculateCapacityVector()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateAndAddConductivityVector(
        VectorType& rRightHandSideVector,
        ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-TransientPwElement::CalculateAndAddPermeabilityFlow()") << std::endl;
        //
        this->CalculateConductivityVector(rVariables.TDimMatrix, rVariables.TMatrix, rVariables.TVector, rVariables);
        //
        //Distribute permeability block vector into elemental vector
        GeoElementUtilities::AssemblePBlockVector<0, TNumNodes>(rRightHandSideVector, rVariables.TVector);
        //
        // KRATOS_INFO("1-TransientPwElement::CalculateAndAddPermeabilityFlow()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateConductivityVector(
        BoundedMatrix<double, TNumNodes, TDim>& TDimMatrix,
        BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
        array_1d<double, TNumNodes>& TVector,
        const ElementVariables& rVariables)
    {
        KRATOS_TRY;
        // KRATOS_INFO("0-UPwSmallStrainElement::CalculatePermeabilityFlow()") << std::endl;
        //
        //this->CalculateThermalDispersionMatrix(rVariables.ConstitutiveMatrix, rVariables);
        //
        TDimMatrix = prod(rVariables.GradNpT, rVariables.ConstitutiveMatrix);
        //
        TMatrix = prod(TDimMatrix, trans(rVariables.GradNpT))
            * rVariables.IntegrationCoefficient;
        //
        noalias(TVector) = - prod(TMatrix, rVariables.TemperatureVector);
        //
        // KRATOS_INFO("1-UPwSmallStrainElement::CalculatePermeabilityFlow()") << std::endl;
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::InitializeProperties(
        ElementVariables& rVariables)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::InitializeProperties") << std::endl;
        //
    	const PropertiesType& rProp = this->GetProperties();
        //
        rVariables.WaterDensity = rProp[DENSITY_WATER];
        rVariables.SolidDensity = rProp[DENSITY_SOLID];
        rVariables.Porosity = rProp[POROSITY];
        rVariables.WaterHeatCapacity = rProp[HEAT_CAPACITY_WATER];
        rVariables.SolidHeatCapacity = rProp[HEAT_CAPACITY_SOLID];
        rVariables.WaterThermalConductivity = rProp[THERMAL_CONDUCTIVITY_WATER];
        rVariables.SolidThermalConductivityXX = rProp[THERMAL_CONDUCTIVITY_SOLID_XX];
        rVariables.SolidThermalConductivityXY = rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
        rVariables.SolidThermalConductivityYX = rProp[THERMAL_CONDUCTIVITY_SOLID_YX];
        rVariables.SolidThermalConductivityYY = rProp[THERMAL_CONDUCTIVITY_SOLID_YY];
        rVariables.Saturation = rProp[SATURATION];
		rVariables.DtTemperatureCoefficient = rProp[DT_TEMPERATURE_COEFFICIENT];
        //
        // KRATOS_INFO("1-TransientThermalElement::InitializeProperties") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateThermalDispersionMatrix(
        BoundedMatrix<double, TDim, TDim>& C,
        ElementVariables& rVariables)
    {
        KRATOS_TRY
    	// KRATOS_INFO("0-TransientThermalElement::CalculateThermalDispersionMatrix") << std::endl;
        //
        const double c0 = rVariables.Porosity * rVariables.Saturation * rVariables.WaterThermalConductivity;
        const double c1 = 1.0 - rVariables.Porosity;
        //
        C(0, 0) = c1 * rVariables.SolidThermalConductivityXX + c0;
        C(0, 1) = c1 * rVariables.SolidThermalConductivityXY;
        C(1, 0) = c1 * rVariables.SolidThermalConductivityYX;
        C(1, 1) = c1 * rVariables.SolidThermalConductivityYY + c0;
        //
        // KRATOS_INFO("0-TransientThermalElement::CalculateThermalDispersionMatrix") << std::endl;
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template<unsigned int TDim, unsigned int TNumNodes>
    void TransientThermalElement<TDim, TNumNodes>::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        //Resetting the LHS
        if (rLeftHandSideMatrix.size1() != N_DOF)
            rLeftHandSideMatrix.resize(N_DOF, N_DOF, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(N_DOF, N_DOF);
        //
        //Resetting the RHS
        if (rRightHandSideVector.size() != N_DOF)
            rRightHandSideVector.resize(N_DOF, false);
        noalias(rRightHandSideVector) = ZeroVector(N_DOF);
        //
        //calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = true;
        //
        CalculateAll(rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        //
    	const unsigned int N_DOF = this->GetNumberOfDOF();
        //
        //Resetting the RHS
        if (rRightHandSideVector.size() != N_DOF)
            rRightHandSideVector.resize(N_DOF, false);
        noalias(rRightHandSideVector) = ZeroVector(N_DOF);
        //
        const bool CalculateStiffnessMatrixFlag = false;
        const bool CalculateResidualVectorFlag = true;
        MatrixType TempMatrix = Matrix();
        //
        CalculateAll(TempMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
        //
        KRATOS_CATCH("")
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    void TransientThermalElement<TDim, TNumNodes>::CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        //
        // Calculation flags
        const bool CalculateStiffnessMatrixFlag = true;
        const bool CalculateResidualVectorFlag = false;
        VectorType TempVector;
        //
        CalculateAll(rLeftHandSideMatrix,
            TempVector,
            rCurrentProcessInfo,
            CalculateStiffnessMatrixFlag,
            CalculateResidualVectorFlag);
        //
        KRATOS_CATCH("");
    }

    // ============================================================================================
    // ============================================================================================
    template< unsigned int TDim, unsigned int TNumNodes >
    GeometryData::IntegrationMethod TransientThermalElement<TDim, TNumNodes>::GetIntegrationMethod() const
    {
        return GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    // ============================================================================================
    // ============================================================================================
    template class TransientThermalElement<2,3>;
    template class TransientThermalElement<2,4>;
    template class TransientThermalElement<3,4>;
    template class TransientThermalElement<3,8>;
    //
    template class TransientThermalElement<2,6>;
    template class TransientThermalElement<2,8>;
    template class TransientThermalElement<2,9>;
    template class TransientThermalElement<3,10>;
    template class TransientThermalElement<3,20>;
    template class TransientThermalElement<3,27>;

} // Namespace Kratos