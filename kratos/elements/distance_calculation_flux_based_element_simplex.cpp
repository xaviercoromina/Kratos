//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez, Pablo Becker
//
#include "distance_calculation_flux_based_element_simplex.h"
#include "includes/checks.h"

namespace Kratos
{

//////////////////////////////////////////////////////////////////////////////
// Public Operations




template <unsigned int TDim >
void DistanceCalculationFluxBasedElement<TDim>::CalculateLocalSystem(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    
    if(step == 1){ //solve a transcient diffusion problem
        CalculatePotentialFlowSystem(rLeftHandSideMatrix, 
                                     rRightHandSideVector,
                                     rCurrentProcessInfo);
    }
    else if (step == 2){ // solve convection + source
        CalculateDistanceSystem(rLeftHandSideMatrix, 
                                     rRightHandSideVector,
                                     rCurrentProcessInfo);
    }
    else{
        KRATOS_ERROR << "FRACTIONAL_STEP must be 1 or 2" << std::endl;
    }

    KRATOS_CATCH("");
}


// template <unsigned int TDim >
// void DistanceCalculationFluxBasedElement<TDim>::CalculateLaplacianSystem(
//         MatrixType &rLeftHandSideMatrix,
//         VectorType &rRightHandSideVector,
//         const ProcessInfo &rCurrentProcessInfo)
// {
//     //getting data for the given geometry
//     double Area;
//     BoundedMatrix<double, NumNodes, TDim > DN_DX;
//     GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


//     //get the previous solution
//     for(unsigned int i=0; i<NumNodes; i++){
//         nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
//     }

//     //compute LHS
//     //add conductivity
//     noalias(rLeftHandSideMatrix) = Area*prod(DN_DX,trans(DN_DX));


//     //substracting previous solution
//     noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,nodal_values);
// }


template <unsigned int TDim >
void DistanceCalculationFluxBasedElement<TDim>::CalculatePotentialFlowSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
{
    //getting data for the given geometry
    double Volume;
    BoundedMatrix<double, NumNodes, TDim > DN_DX;
    array_1d<double, NumNodes > N;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

    array_1d<double, NumNodes > nodal_values;
    for(unsigned int i=0; i<NumNodes; i++){
        nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
    }

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(DENSITY)) << "DENSITY not defined" << std::endl;
    const double density  = rCurrentProcessInfo[DENSITY];

    //compute LHS
    //add conductivity: Like poiseuille flow, the velocity will depend on the square of thickness
    noalias(rLeftHandSideMatrix) = Volume*prod(DN_DX,trans(DN_DX));
    
    //add mass matrix: (Dt = 1)
    const double mass_factor = density*Volume/static_cast<double>(NumNodes);
    for(unsigned int i=0; i<NumNodes; i++){
        double d_gauss = inner_prod(N,nodal_values);
        double target_value = d_gauss > 0.0 ? 1.0 : -1.0;
        rLeftHandSideMatrix(i,i) += mass_factor;
        rRightHandSideVector[i] = mass_factor*target_value;
    }

    //substracting previous solution
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,nodal_values);
}

    
template <unsigned int TDim >
void DistanceCalculationFluxBasedElement<TDim>::CalculateDistanceSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
{
    double Volume;
    BoundedMatrix<double, NumNodes, TDim > DN_DX;
    array_1d<double, NumNodes > N;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

    //get the thickness, previous solution and velocity
    bool has_fixed_node = false; //near the "inlet" we increase conductivity

    array_1d<double, NumNodes > nodal_values;
    array_1d< array_1d<double, 3 >, NumNodes> v; //convection velocity
    array_1d< array_1d<double, 3 >, NumNodes> v_unit; //to decide if convection must be turned off
    array_1d<double, 3 > avg_v_unit = ZeroVector(3); //to decide if convection must be turned off

    for(unsigned int i=0; i<NumNodes; i++){
        nodal_values[i] = GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        v[i] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY); 
        v_unit[i] =  v[i]/MathUtils<double>::Norm3(v[i]); // |v_unit| = 1
        avg_v_unit += v_unit[i];

        if(GetGeometry()[i].IsFixed(DISTANCE)){
            has_fixed_node=true;
        }
    }
    avg_v_unit/=static_cast<double>(NumNodes);

    //computing element size(for Tau)
    const double h = ComputeH(DN_DX);

    //DECIDE IF CONVECTION IS ACTIVE/INACTIVE
   //we decide by checking if flow converges to /diverges from this element
   //if all velocity vectors are oriented in the same direction, then convection is valid
   //but if the element acts as a "sink" or "source", convection is turned off to avoid instabiliities.
    bool add_convection = true;
    double average_unit_vel_misaligment = 0.0;
    for(unsigned int i=0; i<NumNodes; i++) {
        average_unit_vel_misaligment += MathUtils<double>::Norm3( v_unit[i]  - avg_v_unit );
    }
    average_unit_vel_misaligment/=static_cast<double>(NumNodes);
    if(average_unit_vel_misaligment>0.33) add_convection=false;

    //TODO: REMOVE to debug. shows if convection is active or not.
    // this->SetValue(CONVECTION_COEFFICIENT,add_convection );

    //TERMS TO BE ASSEMBLED
    //terms which multiply the gradient of phi
    BoundedMatrix<double,NumNodes, NumNodes> aux2 = ZeroMatrix(NumNodes, NumNodes); //terms multiplying phi
    
    //source term
    array_1d<double, NumNodes> rhs_volumetric_heat = ZeroVector(NumNodes);

    //Gauss point container (using nodal integration)
    const bounded_matrix<double, NumNodes, NumNodes> Ncontainer=IdentityMatrix(NumNodes); 

    //Looping Gauss Points
    for (unsigned int igauss = 0; igauss<NumNodes; igauss++) {
        //Getting the correct GP
        array_1d<double, 4 > N = row(Ncontainer, igauss);

        //Velocity in GP
        array_1d<double, TDim > vel_gauss = ZeroVector(TDim);
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            for (unsigned int k = 0; k<TDim; k++)
                vel_gauss[k] += N[i] * v[i][k];
        }

        const double norm_vel = MathUtils<double>::Norm3(vel_gauss);

        array_1d<double, NumNodes > a_dot_grad = prod(DN_DX, vel_gauss);

        const double d_gauss = inner_prod(N,nodal_values);

        if (d_gauss <= 0.0) {
            vel_gauss = -vel_gauss;
        }

        const double tau_denom = std::max(2.0 * norm_vel / h ,  1e-3); 
        const double tau = 1.0 / (tau_denom);

        double source_term = d_gauss > 0 ? norm_vel : -norm_vel;
        // double source_term = norm_vel;  //to compute flowlength, the source term is exactly norm_vel to get dFlowLength/dx = 1 
        // if(rCurrentProcessInfo[FRACTIONAL_STEP]==3) //to compute pseudofilltime
        //     source_term=1.0; 

        if(add_convection){
            //convection + convection stabilization
            noalias(aux2) += outer_prod(N, a_dot_grad);
            noalias(aux2) += tau*outer_prod(a_dot_grad, a_dot_grad);

            //source
            rhs_volumetric_heat += source_term*N; 
            rhs_volumetric_heat += (tau*source_term)*a_dot_grad;
        }

        //spherical diffusion to enhance stabilization
        double spherical_diffusion = 0.1*tau*(norm_vel*norm_vel);
        noalias(aux2) += spherical_diffusion*prod(DN_DX,trans(DN_DX));
    }

    //adding to system LHS
    noalias(rLeftHandSideMatrix) += aux2;

    //adding to system RHS
    noalias(rRightHandSideVector) += rhs_volumetric_heat; //external forces


    //take out the dirichlet part to finish computing the residual
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, nodal_values);

    //GP weights
    rRightHandSideVector *= Volume/static_cast<double>(NumNodes); 
    rLeftHandSideMatrix *= Volume/static_cast<double>(NumNodes);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <unsigned int TDim >
int DistanceCalculationFluxBasedElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }

    const GeometryType& r_geometry = this->GetGeometry();

    for(unsigned int i=0; i<NumNodes; ++i){
        const Node<3>& rNode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rNode);

        // Check that required dofs exist
        KRATOS_CHECK_DOF_IN_NODE(DISTANCE,rNode);
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( TDim == 2){
        for (unsigned int i=0; i<NumNodes; ++i) {
            if (std::abs(r_geometry[i].Z())>1e-9)
                KRATOS_ERROR << "Node " << r_geometry[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <unsigned int TDim >
std::string DistanceCalculationFluxBasedElement<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "DistanceCalculationFluxBasedElement" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <unsigned int TDim >
void DistanceCalculationFluxBasedElement<TDim>::PrintInfo(
    std::ostream &rOStream) const
{
    rOStream << this->Info() << std::endl;
}

template< unsigned int TDim  >
void DistanceCalculationFluxBasedElement< TDim >::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int i = 0; i < NumNodes; ++i){
        rResult[LocalIndex++] = r_geometry[i].GetDof(DISTANCE).EquationId();
    }
}


template< unsigned int TDim  >
void DistanceCalculationFluxBasedElement< TDim >::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

     if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

     unsigned int LocalIndex = 0;

     for (unsigned int i = 0; i < NumNodes; ++i){
         rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(DISTANCE);
     }
}

template< unsigned int TDim  >
double DistanceCalculationFluxBasedElement< TDim >::ComputeH(BoundedMatrix<double,NumNodes, TDim>& DN_DX)
{
        double h=0.0;
        for(unsigned int i=0; i<NumNodes; i++){
            double h_inv = 0.0;
            for(unsigned int k=0; k<TDim; k++)
            {
                h_inv += DN_DX(i,k)*DN_DX(i,k);
            }
            h = std::max(h, 1.0 / h_inv);
        }
        h = std::sqrt(h);
        return h;
}

template < unsigned int TDim >
void DistanceCalculationFluxBasedElement< TDim >::AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo)
{
   	GeometryType& rGeometry = this->GetGeometry();

	//getting data for the given geometry
	BoundedMatrix<double, NumNodes, TDim > DN_DX;
    array_1d<double, NumNodes > N;
	double Volume;
	GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, Volume);

	//get the nodal values
	array_1d<double, NumNodes > values;
	for(unsigned int i=0; i< NumNodes; i++){
		values[i] = rGeometry[i].FastGetSolutionStepValue(DISTANCE);
	}
    //computing the gradient to see flow direction
    const array_1d<double, TDim> grad = prod(trans(DN_DX), values);
	const array_1d<double, TDim> unit_grad = (1.0/MathUtils<double>::Norm(grad)) * grad;
	const array_1d<double, TDim> vel = unit_grad;

	//saving data
	const double vol_factor =  Volume/static_cast<double>(NumNodes);
	for (unsigned int j = 0; j < NumNodes; j++){ //looping 4 nodes of the elem:
		rGeometry[j].SetLock();
		rGeometry[j].FastGetSolutionStepValue(NODAL_VOLUME)+=vol_factor;
		rGeometry[j].FastGetSolutionStepValue(VELOCITY)+=vel*vol_factor;
		rGeometry[j].UnSetLock();
	}
}


// Class template instantiation
template class  KRATOS_API(KRATOS_CORE) DistanceCalculationFluxBasedElement<2>;
template class  KRATOS_API(KRATOS_CORE) DistanceCalculationFluxBasedElement<3>;


} // namespace Kratos