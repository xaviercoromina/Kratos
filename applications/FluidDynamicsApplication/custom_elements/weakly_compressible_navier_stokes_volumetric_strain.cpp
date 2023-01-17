//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes

// Application includes
#include "weakly_compressible_navier_stokes_volumetric_strain.h"
#include "data_containers/weakly_compressible/weakly_compressible_navier_stokes_volumetric_strain_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::WeaklyCompressibleNavierStokesVolumetricStrain(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::WeaklyCompressibleNavierStokesVolumetricStrain(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::WeaklyCompressibleNavierStokesVolumetricStrain(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::WeaklyCompressibleNavierStokesVolumetricStrain(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::~WeaklyCompressibleNavierStokesVolumetricStrain() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokesVolumetricStrain>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokesVolumetricStrain>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template<class TElementData>
const Parameters WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["VORTICITY","Q_VALUE","VORTICITY_MAGNITUDE"],
            "nodal_historical"       : ["VELOCITY","VOLUMETRIC_STRAIN"],
            "nodal_non_historical"   : ["DENSITY","DYNAMIC_VISCOSITY","BULK_MODULUS","SOUND_VELOCITY"],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","IS_STRUCTURE","DISPLACEMENT","BODY_FORCE","NODAL_AREA","NODAL_H","ADVPROJ","DIVPROJ","REACTION","REACTION_WATER_PRESSURE","EXTERNAL_PRESSURE","NORMAL","Y_WALL","Q_VALUE"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["Newtonian2DLaw","Newtonian3DLaw","NewtonianTemperatureDependent2DLaw","NewtonianTemperatureDependent3DLaw","Euler2DLaw","Euler3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This implements a weakly compressible Navier-Stokes element with quasi-static Variational MultiScales (VMS) stabilization. Note that this formulation allows a weak coupling with a custom Equation of State that updates the nodal DENSITY from the obtained PRESSURE values. Also note that no viscous behavior is hardcoded, meaning that any fluid constitutive model can be used through a constitutive law."
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TElementData>
std::string WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "WeaklyCompressibleNavierStokesVolumetricStrain" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::AddBoundaryTraction(
    TElementData& rData,
    const Vector& rUnitNormal,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    // Set the current Gauss pt. Voigt notation normal projection matrix
    BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
    FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal, voigt_normal_projection_matrix);

    // Set the current Gauss pt. strain matrix
    BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX, B_matrix);

    // Compute some Gauss pt. auxiliar matrices
    const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, rData.C);
    const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

    // Fill the pressure to Voigt notation operator matrix
    BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            pres_to_voigt_matrix_op(comp, i*BlockSize+Dim) = rData.N[i];
        }
    }

    // Set the shape functions auxiliar transpose matrix
    BoundedMatrix<double, LocalSize, Dim> N_aux_trans = ZeroMatrix(LocalSize, Dim);
    for (unsigned int i=0; i<NumNodes; ++i) {
        for (unsigned int comp=0; comp<Dim; ++comp) {
            N_aux_trans(i*BlockSize+comp, comp) = rData.N[i];
        }
    }

    // Contribution coming fron the shear stress operator
    noalias(rData.lhs) = prod(N_aux_trans, aux_matrix_ACB);

    // Contribution coming from the pressure terms
    const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
    noalias(rData.lhs) -= prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);

    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    rData.lhs *= rData.Weight;
    noalias(rLHS) -= rData.lhs;
    noalias(rRHS) += prod(rData.lhs,values);
}

template <>
void WeaklyCompressibleNavierStokesVolumetricStrain< WeaklyCompressibleNavierStokesVolumetricStrainData<2,3> >::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesVolumetricStrainData<2,3>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double k = rData.BulkModulus;
    const double mu = rData.EffectiveViscosity;

    const double bdf0 = rData.bdf0;
    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,3,3>& C = rData.C;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs4 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs5 = DN(0,0)*clhs3 + DN(0,1)*clhs4;
const double clhs6 = N[0]*rho;
const double clhs7 = pow(N[0], 2)*bdf0;
const double clhs8 = N[0]*bdf0;
const double clhs9 = clhs5 + clhs8;
const double clhs10 = 1.0/(rho*stab_c2*sqrt(pow(clhs3, 2) + pow(clhs4, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs11 = clhs10*pow(rho, 2);
const double clhs12 = clhs11*clhs5;
const double clhs13 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs14 = clhs11*clhs13;
const double clhs15 = N[0]*clhs14;
const double clhs16 = clhs12*clhs9 + clhs15*clhs9 + clhs5*clhs6 + clhs7*rho;
const double clhs17 = C(0,1)*DN(0,1) + clhs1;
const double clhs18 = C(1,2)*DN(0,1);
const double clhs19 = C(2,2)*DN(0,0) + clhs18;
const double clhs20 = clhs10*clhs13;
const double clhs21 = clhs10*rho;
const double clhs22 = clhs21*clhs5;
const double clhs23 = k*(N[0] - clhs20*clhs6 - clhs22);
const double clhs24 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs25 = C(0,2)*DN(1,0);
const double clhs26 = C(2,2)*DN(1,1) + clhs25;
const double clhs27 = DN(1,0)*clhs3 + DN(1,1)*clhs4;
const double clhs28 = N[1]*clhs8;
const double clhs29 = clhs28*rho;
const double clhs30 = N[1]*bdf0;
const double clhs31 = clhs27 + clhs30;
const double clhs32 = clhs12*clhs31 + clhs15*clhs31 + clhs27*clhs6 + clhs29;
const double clhs33 = C(0,1)*DN(1,1) + clhs25;
const double clhs34 = C(1,2)*DN(1,1);
const double clhs35 = C(2,2)*DN(1,0) + clhs34;
const double clhs36 = DN(1,0)*N[0];
const double clhs37 = clhs13*clhs21;
const double clhs38 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs39 = C(0,2)*DN(2,0);
const double clhs40 = C(2,2)*DN(2,1) + clhs39;
const double clhs41 = DN(2,0)*clhs3 + DN(2,1)*clhs4;
const double clhs42 = N[2]*clhs8;
const double clhs43 = clhs42*rho;
const double clhs44 = N[2]*bdf0 + clhs41;
const double clhs45 = clhs12*clhs44 + clhs15*clhs44 + clhs41*clhs6 + clhs43;
const double clhs46 = C(0,1)*DN(2,1) + clhs39;
const double clhs47 = C(1,2)*DN(2,1);
const double clhs48 = C(2,2)*DN(2,0) + clhs47;
const double clhs49 = DN(2,0)*N[0];
const double clhs50 = C(0,1)*DN(0,0) + clhs18;
const double clhs51 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs52 = C(0,1)*DN(1,0) + clhs34;
const double clhs53 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs54 = DN(1,1)*N[0];
const double clhs55 = C(0,1)*DN(2,0) + clhs47;
const double clhs56 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs57 = DN(2,1)*N[0];
const double clhs58 = clhs21*clhs9;
const double clhs59 = N[0] + clhs58;
const double clhs60 = clhs10*k;
const double clhs61 = clhs21*clhs31;
const double clhs62 = DN(0,0)*clhs60;
const double clhs63 = DN(0,1)*clhs60;
const double clhs64 = DN(1,0)*clhs62 + DN(1,1)*clhs63 + clhs28;
const double clhs65 = clhs21*clhs44;
const double clhs66 = DN(2,0)*clhs62 + DN(2,1)*clhs63 + clhs42;
const double clhs67 = N[1]*rho;
const double clhs68 = clhs11*clhs27;
const double clhs69 = N[1]*clhs14;
const double clhs70 = clhs29 + clhs5*clhs67 + clhs68*clhs9 + clhs69*clhs9;
const double clhs71 = DN(0,0)*N[1];
const double clhs72 = clhs21*clhs27;
const double clhs73 = pow(N[1], 2)*bdf0;
const double clhs74 = clhs27*clhs67 + clhs31*clhs68 + clhs31*clhs69 + clhs73*rho;
const double clhs75 = k*(N[1] - clhs20*clhs67 - clhs72);
const double clhs76 = N[2]*clhs30;
const double clhs77 = clhs76*rho;
const double clhs78 = clhs41*clhs67 + clhs44*clhs68 + clhs44*clhs69 + clhs77;
const double clhs79 = DN(2,0)*N[1];
const double clhs80 = DN(0,1)*N[1];
const double clhs81 = DN(2,1)*N[1];
const double clhs82 = N[1] + clhs61;
const double clhs83 = DN(1,0)*DN(2,0)*clhs60 + DN(1,1)*DN(2,1)*clhs60 + clhs76;
const double clhs84 = N[2]*rho;
const double clhs85 = clhs11*clhs41;
const double clhs86 = N[2]*clhs14;
const double clhs87 = clhs43 + clhs5*clhs84 + clhs85*clhs9 + clhs86*clhs9;
const double clhs88 = DN(0,0)*N[2];
const double clhs89 = clhs21*clhs41;
const double clhs90 = clhs27*clhs84 + clhs31*clhs85 + clhs31*clhs86 + clhs77;
const double clhs91 = DN(1,0)*N[2];
const double clhs92 = pow(N[2], 2)*bdf0;
const double clhs93 = clhs41*clhs84 + clhs44*clhs85 + clhs44*clhs86 + clhs92*rho;
const double clhs94 = k*(N[2] - clhs20*clhs84 - clhs89);
const double clhs95 = DN(0,1)*N[2];
const double clhs96 = DN(1,1)*N[2];
const double clhs97 = N[2] + clhs65;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs16;
lhs(0,1)=DN(0,0)*clhs17 + DN(0,1)*clhs19;
lhs(0,2)=DN(0,0)*clhs23;
lhs(0,3)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + clhs32;
lhs(0,4)=DN(0,0)*clhs33 + DN(0,1)*clhs35;
lhs(0,5)=k*(DN(0,0)*N[1] - DN(1,0)*clhs22 - clhs36*clhs37);
lhs(0,6)=DN(0,0)*clhs38 + DN(0,1)*clhs40 + clhs45;
lhs(0,7)=DN(0,0)*clhs46 + DN(0,1)*clhs48;
lhs(0,8)=k*(DN(0,0)*N[2] - DN(2,0)*clhs22 - clhs37*clhs49);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs50;
lhs(1,1)=DN(0,0)*clhs19 + DN(0,1)*clhs51 + clhs16;
lhs(1,2)=DN(0,1)*clhs23;
lhs(1,3)=DN(0,0)*clhs26 + DN(0,1)*clhs52;
lhs(1,4)=DN(0,0)*clhs35 + DN(0,1)*clhs53 + clhs32;
lhs(1,5)=k*(DN(0,1)*N[1] - DN(1,1)*clhs22 - clhs37*clhs54);
lhs(1,6)=DN(0,0)*clhs40 + DN(0,1)*clhs55;
lhs(1,7)=DN(0,0)*clhs48 + DN(0,1)*clhs56 + clhs45;
lhs(1,8)=k*(DN(0,1)*N[2] - DN(2,1)*clhs22 - clhs37*clhs57);
lhs(2,0)=-DN(0,0)*clhs59;
lhs(2,1)=-DN(0,1)*clhs59;
lhs(2,2)=pow(DN(0,0), 2)*clhs60 + pow(DN(0,1), 2)*clhs60 + clhs7;
lhs(2,3)=-DN(0,0)*clhs61 - clhs36;
lhs(2,4)=-DN(0,1)*clhs61 - clhs54;
lhs(2,5)=clhs64;
lhs(2,6)=-DN(0,0)*clhs65 - clhs49;
lhs(2,7)=-DN(0,1)*clhs65 - clhs57;
lhs(2,8)=clhs66;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs70;
lhs(3,1)=DN(1,0)*clhs17 + DN(1,1)*clhs19;
lhs(3,2)=k*(-DN(0,0)*clhs72 + DN(1,0)*N[0] - clhs37*clhs71);
lhs(3,3)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + clhs74;
lhs(3,4)=DN(1,0)*clhs33 + DN(1,1)*clhs35;
lhs(3,5)=DN(1,0)*clhs75;
lhs(3,6)=DN(1,0)*clhs38 + DN(1,1)*clhs40 + clhs78;
lhs(3,7)=DN(1,0)*clhs46 + DN(1,1)*clhs48;
lhs(3,8)=k*(DN(1,0)*N[2] - DN(2,0)*clhs72 - clhs37*clhs79);
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs50;
lhs(4,1)=DN(1,0)*clhs19 + DN(1,1)*clhs51 + clhs70;
lhs(4,2)=k*(-DN(0,1)*clhs72 + DN(1,1)*N[0] - clhs37*clhs80);
lhs(4,3)=DN(1,0)*clhs26 + DN(1,1)*clhs52;
lhs(4,4)=DN(1,0)*clhs35 + DN(1,1)*clhs53 + clhs74;
lhs(4,5)=DN(1,1)*clhs75;
lhs(4,6)=DN(1,0)*clhs40 + DN(1,1)*clhs55;
lhs(4,7)=DN(1,0)*clhs48 + DN(1,1)*clhs56 + clhs78;
lhs(4,8)=k*(DN(1,1)*N[2] - DN(2,1)*clhs72 - clhs37*clhs81);
lhs(5,0)=-DN(1,0)*clhs58 - clhs71;
lhs(5,1)=-DN(1,1)*clhs58 - clhs80;
lhs(5,2)=clhs64;
lhs(5,3)=-DN(1,0)*clhs82;
lhs(5,4)=-DN(1,1)*clhs82;
lhs(5,5)=pow(DN(1,0), 2)*clhs60 + pow(DN(1,1), 2)*clhs60 + clhs73;
lhs(5,6)=-DN(1,0)*clhs65 - clhs79;
lhs(5,7)=-DN(1,1)*clhs65 - clhs81;
lhs(5,8)=clhs83;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs87;
lhs(6,1)=DN(2,0)*clhs17 + DN(2,1)*clhs19;
lhs(6,2)=k*(-DN(0,0)*clhs89 + DN(2,0)*N[0] - clhs37*clhs88);
lhs(6,3)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + clhs90;
lhs(6,4)=DN(2,0)*clhs33 + DN(2,1)*clhs35;
lhs(6,5)=k*(-DN(1,0)*clhs89 + DN(2,0)*N[1] - clhs37*clhs91);
lhs(6,6)=DN(2,0)*clhs38 + DN(2,1)*clhs40 + clhs93;
lhs(6,7)=DN(2,0)*clhs46 + DN(2,1)*clhs48;
lhs(6,8)=DN(2,0)*clhs94;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs50;
lhs(7,1)=DN(2,0)*clhs19 + DN(2,1)*clhs51 + clhs87;
lhs(7,2)=k*(-DN(0,1)*clhs89 + DN(2,1)*N[0] - clhs37*clhs95);
lhs(7,3)=DN(2,0)*clhs26 + DN(2,1)*clhs52;
lhs(7,4)=DN(2,0)*clhs35 + DN(2,1)*clhs53 + clhs90;
lhs(7,5)=k*(-DN(1,1)*clhs89 + DN(2,1)*N[1] - clhs37*clhs96);
lhs(7,6)=DN(2,0)*clhs40 + DN(2,1)*clhs55;
lhs(7,7)=DN(2,0)*clhs48 + DN(2,1)*clhs56 + clhs93;
lhs(7,8)=DN(2,1)*clhs94;
lhs(8,0)=-DN(2,0)*clhs58 - clhs88;
lhs(8,1)=-DN(2,1)*clhs58 - clhs95;
lhs(8,2)=clhs66;
lhs(8,3)=-DN(2,0)*clhs61 - clhs91;
lhs(8,4)=-DN(2,1)*clhs61 - clhs96;
lhs(8,5)=clhs83;
lhs(8,6)=-DN(2,0)*clhs97;
lhs(8,7)=-DN(2,1)*clhs97;
lhs(8,8)=pow(DN(2,0), 2)*clhs60 + pow(DN(2,1), 2)*clhs60 + clhs92;


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void WeaklyCompressibleNavierStokesVolumetricStrain<WeaklyCompressibleNavierStokesVolumetricStrainData<3,4>>::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesVolumetricStrainData<3,4>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double k = rData.BulkModulus;
    const double mu = rData.EffectiveViscosity;

    const double bdf0 = rData.bdf0;
    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,3,4> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,6,6>& C = rData.C;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rLeftHandSideMatrix
    auto& lhs = rData.lhs;

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs6 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs7 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs8 = DN(0,0)*clhs5 + DN(0,1)*clhs6 + DN(0,2)*clhs7;
const double clhs9 = N[0]*rho;
const double clhs10 = pow(N[0], 2)*bdf0;
const double clhs11 = N[0]*bdf0;
const double clhs12 = clhs11 + clhs8;
const double clhs13 = 1.0/(rho*stab_c2*sqrt(pow(clhs5, 2) + pow(clhs6, 2) + pow(clhs7, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs14 = clhs13*pow(rho, 2);
const double clhs15 = clhs14*clhs8;
const double clhs16 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs17 = clhs14*clhs16;
const double clhs18 = N[0]*clhs17;
const double clhs19 = clhs10*rho + clhs12*clhs15 + clhs12*clhs18 + clhs8*clhs9;
const double clhs20 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs21 = C(1,3)*DN(0,1);
const double clhs22 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs21;
const double clhs23 = C(3,5)*DN(0,0);
const double clhs24 = C(4,5)*DN(0,2);
const double clhs25 = C(1,5)*DN(0,1) + clhs23 + clhs24;
const double clhs26 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs27 = C(3,4)*DN(0,1);
const double clhs28 = C(2,3)*DN(0,2) + clhs23 + clhs27;
const double clhs29 = C(2,5)*DN(0,2);
const double clhs30 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs29;
const double clhs31 = clhs13*clhs16;
const double clhs32 = clhs13*rho;
const double clhs33 = clhs32*clhs8;
const double clhs34 = k*(N[0] - clhs31*clhs9 - clhs33);
const double clhs35 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs36 = C(0,3)*DN(1,0);
const double clhs37 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs36;
const double clhs38 = C(0,5)*DN(1,0);
const double clhs39 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs38;
const double clhs40 = DN(1,0)*clhs5 + DN(1,1)*clhs6 + DN(1,2)*clhs7;
const double clhs41 = N[1]*clhs11;
const double clhs42 = clhs41*rho;
const double clhs43 = N[1]*bdf0;
const double clhs44 = clhs40 + clhs43;
const double clhs45 = clhs15*clhs44 + clhs18*clhs44 + clhs40*clhs9 + clhs42;
const double clhs46 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs36;
const double clhs47 = C(1,3)*DN(1,1);
const double clhs48 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs47;
const double clhs49 = C(3,5)*DN(1,0);
const double clhs50 = C(4,5)*DN(1,2);
const double clhs51 = C(1,5)*DN(1,1) + clhs49 + clhs50;
const double clhs52 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs38;
const double clhs53 = C(3,4)*DN(1,1);
const double clhs54 = C(2,3)*DN(1,2) + clhs49 + clhs53;
const double clhs55 = C(2,5)*DN(1,2);
const double clhs56 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs55;
const double clhs57 = DN(1,0)*N[0];
const double clhs58 = clhs16*clhs32;
const double clhs59 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs60 = C(0,3)*DN(2,0);
const double clhs61 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs60;
const double clhs62 = C(0,5)*DN(2,0);
const double clhs63 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs62;
const double clhs64 = DN(2,0)*clhs5 + DN(2,1)*clhs6 + DN(2,2)*clhs7;
const double clhs65 = N[2]*clhs11;
const double clhs66 = clhs65*rho;
const double clhs67 = N[2]*bdf0;
const double clhs68 = clhs64 + clhs67;
const double clhs69 = clhs15*clhs68 + clhs18*clhs68 + clhs64*clhs9 + clhs66;
const double clhs70 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs60;
const double clhs71 = C(1,3)*DN(2,1);
const double clhs72 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs71;
const double clhs73 = C(3,5)*DN(2,0);
const double clhs74 = C(4,5)*DN(2,2);
const double clhs75 = C(1,5)*DN(2,1) + clhs73 + clhs74;
const double clhs76 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs62;
const double clhs77 = C(3,4)*DN(2,1);
const double clhs78 = C(2,3)*DN(2,2) + clhs73 + clhs77;
const double clhs79 = C(2,5)*DN(2,2);
const double clhs80 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs79;
const double clhs81 = DN(2,0)*N[0];
const double clhs82 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs83 = C(0,3)*DN(3,0);
const double clhs84 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs83;
const double clhs85 = C(0,5)*DN(3,0);
const double clhs86 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs85;
const double clhs87 = DN(3,0)*clhs5 + DN(3,1)*clhs6 + DN(3,2)*clhs7;
const double clhs88 = N[3]*clhs11;
const double clhs89 = clhs88*rho;
const double clhs90 = N[3]*bdf0 + clhs87;
const double clhs91 = clhs15*clhs90 + clhs18*clhs90 + clhs87*clhs9 + clhs89;
const double clhs92 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs83;
const double clhs93 = C(1,3)*DN(3,1);
const double clhs94 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs93;
const double clhs95 = C(3,5)*DN(3,0);
const double clhs96 = C(4,5)*DN(3,2);
const double clhs97 = C(1,5)*DN(3,1) + clhs95 + clhs96;
const double clhs98 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs85;
const double clhs99 = C(3,4)*DN(3,1);
const double clhs100 = C(2,3)*DN(3,2) + clhs95 + clhs99;
const double clhs101 = C(2,5)*DN(3,2);
const double clhs102 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs101;
const double clhs103 = DN(3,0)*N[0];
const double clhs104 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs21;
const double clhs105 = C(0,4)*DN(0,0) + clhs24 + clhs27;
const double clhs106 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs107 = C(1,4)*DN(0,1);
const double clhs108 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs107;
const double clhs109 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs107;
const double clhs110 = C(2,4)*DN(0,2);
const double clhs111 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs110;
const double clhs112 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs47;
const double clhs113 = C(0,4)*DN(1,0) + clhs50 + clhs53;
const double clhs114 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs115 = C(1,4)*DN(1,1);
const double clhs116 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs115;
const double clhs117 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs115;
const double clhs118 = C(2,4)*DN(1,2);
const double clhs119 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs118;
const double clhs120 = DN(1,1)*N[0];
const double clhs121 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs71;
const double clhs122 = C(0,4)*DN(2,0) + clhs74 + clhs77;
const double clhs123 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs124 = C(1,4)*DN(2,1);
const double clhs125 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs124;
const double clhs126 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs124;
const double clhs127 = C(2,4)*DN(2,2);
const double clhs128 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs127;
const double clhs129 = DN(2,1)*N[0];
const double clhs130 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs93;
const double clhs131 = C(0,4)*DN(3,0) + clhs96 + clhs99;
const double clhs132 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs133 = C(1,4)*DN(3,1);
const double clhs134 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs133;
const double clhs135 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs133;
const double clhs136 = C(2,4)*DN(3,2);
const double clhs137 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs136;
const double clhs138 = DN(3,1)*N[0];
const double clhs139 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs29;
const double clhs140 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs110;
const double clhs141 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs142 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs55;
const double clhs143 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs118;
const double clhs144 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs145 = DN(1,2)*N[0];
const double clhs146 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs79;
const double clhs147 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs127;
const double clhs148 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs149 = DN(2,2)*N[0];
const double clhs150 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs101;
const double clhs151 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs136;
const double clhs152 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs153 = DN(3,2)*N[0];
const double clhs154 = clhs12*clhs32;
const double clhs155 = N[0] + clhs154;
const double clhs156 = clhs13*k;
const double clhs157 = clhs32*clhs44;
const double clhs158 = DN(0,0)*clhs156;
const double clhs159 = DN(0,1)*clhs156;
const double clhs160 = DN(0,2)*clhs156;
const double clhs161 = DN(1,0)*clhs158 + DN(1,1)*clhs159 + DN(1,2)*clhs160 + clhs41;
const double clhs162 = clhs32*clhs68;
const double clhs163 = DN(2,0)*clhs158 + DN(2,1)*clhs159 + DN(2,2)*clhs160 + clhs65;
const double clhs164 = clhs32*clhs90;
const double clhs165 = DN(3,0)*clhs158 + DN(3,1)*clhs159 + DN(3,2)*clhs160 + clhs88;
const double clhs166 = N[1]*rho;
const double clhs167 = clhs14*clhs40;
const double clhs168 = N[1]*clhs17;
const double clhs169 = clhs12*clhs167 + clhs12*clhs168 + clhs166*clhs8 + clhs42;
const double clhs170 = DN(0,0)*N[1];
const double clhs171 = clhs32*clhs40;
const double clhs172 = pow(N[1], 2)*bdf0;
const double clhs173 = clhs166*clhs40 + clhs167*clhs44 + clhs168*clhs44 + clhs172*rho;
const double clhs174 = k*(N[1] - clhs166*clhs31 - clhs171);
const double clhs175 = N[2]*clhs43;
const double clhs176 = clhs175*rho;
const double clhs177 = clhs166*clhs64 + clhs167*clhs68 + clhs168*clhs68 + clhs176;
const double clhs178 = DN(2,0)*N[1];
const double clhs179 = N[3]*clhs43;
const double clhs180 = clhs179*rho;
const double clhs181 = clhs166*clhs87 + clhs167*clhs90 + clhs168*clhs90 + clhs180;
const double clhs182 = DN(3,0)*N[1];
const double clhs183 = DN(0,1)*N[1];
const double clhs184 = DN(2,1)*N[1];
const double clhs185 = DN(3,1)*N[1];
const double clhs186 = DN(0,2)*N[1];
const double clhs187 = DN(2,2)*N[1];
const double clhs188 = DN(3,2)*N[1];
const double clhs189 = N[1] + clhs157;
const double clhs190 = DN(1,0)*clhs156;
const double clhs191 = DN(1,1)*clhs156;
const double clhs192 = DN(1,2)*clhs156;
const double clhs193 = DN(2,0)*clhs190 + DN(2,1)*clhs191 + DN(2,2)*clhs192 + clhs175;
const double clhs194 = DN(3,0)*clhs190 + DN(3,1)*clhs191 + DN(3,2)*clhs192 + clhs179;
const double clhs195 = N[2]*rho;
const double clhs196 = clhs14*clhs64;
const double clhs197 = N[2]*clhs17;
const double clhs198 = clhs12*clhs196 + clhs12*clhs197 + clhs195*clhs8 + clhs66;
const double clhs199 = DN(0,0)*N[2];
const double clhs200 = clhs32*clhs64;
const double clhs201 = clhs176 + clhs195*clhs40 + clhs196*clhs44 + clhs197*clhs44;
const double clhs202 = DN(1,0)*N[2];
const double clhs203 = pow(N[2], 2)*bdf0;
const double clhs204 = clhs195*clhs64 + clhs196*clhs68 + clhs197*clhs68 + clhs203*rho;
const double clhs205 = k*(N[2] - clhs195*clhs31 - clhs200);
const double clhs206 = N[3]*clhs67;
const double clhs207 = clhs206*rho;
const double clhs208 = clhs195*clhs87 + clhs196*clhs90 + clhs197*clhs90 + clhs207;
const double clhs209 = DN(3,0)*N[2];
const double clhs210 = DN(0,1)*N[2];
const double clhs211 = DN(1,1)*N[2];
const double clhs212 = DN(3,1)*N[2];
const double clhs213 = DN(0,2)*N[2];
const double clhs214 = DN(1,2)*N[2];
const double clhs215 = DN(3,2)*N[2];
const double clhs216 = N[2] + clhs162;
const double clhs217 = DN(2,0)*DN(3,0)*clhs156 + DN(2,1)*DN(3,1)*clhs156 + DN(2,2)*DN(3,2)*clhs156 + clhs206;
const double clhs218 = N[3]*rho;
const double clhs219 = clhs14*clhs87;
const double clhs220 = N[3]*clhs17;
const double clhs221 = clhs12*clhs219 + clhs12*clhs220 + clhs218*clhs8 + clhs89;
const double clhs222 = DN(0,0)*N[3];
const double clhs223 = clhs32*clhs87;
const double clhs224 = clhs180 + clhs218*clhs40 + clhs219*clhs44 + clhs220*clhs44;
const double clhs225 = DN(1,0)*N[3];
const double clhs226 = clhs207 + clhs218*clhs64 + clhs219*clhs68 + clhs220*clhs68;
const double clhs227 = DN(2,0)*N[3];
const double clhs228 = pow(N[3], 2)*bdf0;
const double clhs229 = clhs218*clhs87 + clhs219*clhs90 + clhs220*clhs90 + clhs228*rho;
const double clhs230 = k*(N[3] - clhs218*clhs31 - clhs223);
const double clhs231 = DN(0,1)*N[3];
const double clhs232 = DN(1,1)*N[3];
const double clhs233 = DN(2,1)*N[3];
const double clhs234 = DN(0,2)*N[3];
const double clhs235 = DN(1,2)*N[3];
const double clhs236 = DN(2,2)*N[3];
const double clhs237 = N[3] + clhs164;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs19;
lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + DN(0,2)*clhs25;
lhs(0,2)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + DN(0,2)*clhs30;
lhs(0,3)=DN(0,0)*clhs34;
lhs(0,4)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + DN(0,2)*clhs39 + clhs45;
lhs(0,5)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs51;
lhs(0,6)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + DN(0,2)*clhs56;
lhs(0,7)=k*(DN(0,0)*N[1] - DN(1,0)*clhs33 - clhs57*clhs58);
lhs(0,8)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + DN(0,2)*clhs63 + clhs69;
lhs(0,9)=DN(0,0)*clhs70 + DN(0,1)*clhs72 + DN(0,2)*clhs75;
lhs(0,10)=DN(0,0)*clhs76 + DN(0,1)*clhs78 + DN(0,2)*clhs80;
lhs(0,11)=k*(DN(0,0)*N[2] - DN(2,0)*clhs33 - clhs58*clhs81);
lhs(0,12)=DN(0,0)*clhs82 + DN(0,1)*clhs84 + DN(0,2)*clhs86 + clhs91;
lhs(0,13)=DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs97;
lhs(0,14)=DN(0,0)*clhs98 + DN(0,1)*clhs100 + DN(0,2)*clhs102;
lhs(0,15)=k*(DN(0,0)*N[3] - DN(3,0)*clhs33 - clhs103*clhs58);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs104 + DN(0,2)*clhs105;
lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs106 + DN(0,2)*clhs108 + clhs19;
lhs(1,2)=DN(0,0)*clhs28 + DN(0,1)*clhs109 + DN(0,2)*clhs111;
lhs(1,3)=DN(0,1)*clhs34;
lhs(1,4)=DN(0,0)*clhs37 + DN(0,1)*clhs112 + DN(0,2)*clhs113;
lhs(1,5)=DN(0,0)*clhs48 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs45;
lhs(1,6)=DN(0,0)*clhs54 + DN(0,1)*clhs117 + DN(0,2)*clhs119;
lhs(1,7)=k*(DN(0,1)*N[1] - DN(1,1)*clhs33 - clhs120*clhs58);
lhs(1,8)=DN(0,0)*clhs61 + DN(0,1)*clhs121 + DN(0,2)*clhs122;
lhs(1,9)=DN(0,0)*clhs72 + DN(0,1)*clhs123 + DN(0,2)*clhs125 + clhs69;
lhs(1,10)=DN(0,0)*clhs78 + DN(0,1)*clhs126 + DN(0,2)*clhs128;
lhs(1,11)=k*(DN(0,1)*N[2] - DN(2,1)*clhs33 - clhs129*clhs58);
lhs(1,12)=DN(0,0)*clhs84 + DN(0,1)*clhs130 + DN(0,2)*clhs131;
lhs(1,13)=DN(0,0)*clhs94 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs91;
lhs(1,14)=DN(0,0)*clhs100 + DN(0,1)*clhs135 + DN(0,2)*clhs137;
lhs(1,15)=k*(DN(0,1)*N[3] - DN(3,1)*clhs33 - clhs138*clhs58);
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs105 + DN(0,2)*clhs139;
lhs(2,1)=DN(0,0)*clhs25 + DN(0,1)*clhs108 + DN(0,2)*clhs140;
lhs(2,2)=DN(0,0)*clhs30 + DN(0,1)*clhs111 + DN(0,2)*clhs141 + clhs19;
lhs(2,3)=DN(0,2)*clhs34;
lhs(2,4)=DN(0,0)*clhs39 + DN(0,1)*clhs113 + DN(0,2)*clhs142;
lhs(2,5)=DN(0,0)*clhs51 + DN(0,1)*clhs116 + DN(0,2)*clhs143;
lhs(2,6)=DN(0,0)*clhs56 + DN(0,1)*clhs119 + DN(0,2)*clhs144 + clhs45;
lhs(2,7)=k*(DN(0,2)*N[1] - DN(1,2)*clhs33 - clhs145*clhs58);
lhs(2,8)=DN(0,0)*clhs63 + DN(0,1)*clhs122 + DN(0,2)*clhs146;
lhs(2,9)=DN(0,0)*clhs75 + DN(0,1)*clhs125 + DN(0,2)*clhs147;
lhs(2,10)=DN(0,0)*clhs80 + DN(0,1)*clhs128 + DN(0,2)*clhs148 + clhs69;
lhs(2,11)=k*(DN(0,2)*N[2] - DN(2,2)*clhs33 - clhs149*clhs58);
lhs(2,12)=DN(0,0)*clhs86 + DN(0,1)*clhs131 + DN(0,2)*clhs150;
lhs(2,13)=DN(0,0)*clhs97 + DN(0,1)*clhs134 + DN(0,2)*clhs151;
lhs(2,14)=DN(0,0)*clhs102 + DN(0,1)*clhs137 + DN(0,2)*clhs152 + clhs91;
lhs(2,15)=k*(DN(0,2)*N[3] - DN(3,2)*clhs33 - clhs153*clhs58);
lhs(3,0)=-DN(0,0)*clhs155;
lhs(3,1)=-DN(0,1)*clhs155;
lhs(3,2)=-DN(0,2)*clhs155;
lhs(3,3)=pow(DN(0,0), 2)*clhs156 + pow(DN(0,1), 2)*clhs156 + pow(DN(0,2), 2)*clhs156 + clhs10;
lhs(3,4)=-DN(0,0)*clhs157 - clhs57;
lhs(3,5)=-DN(0,1)*clhs157 - clhs120;
lhs(3,6)=-DN(0,2)*clhs157 - clhs145;
lhs(3,7)=clhs161;
lhs(3,8)=-DN(0,0)*clhs162 - clhs81;
lhs(3,9)=-DN(0,1)*clhs162 - clhs129;
lhs(3,10)=-DN(0,2)*clhs162 - clhs149;
lhs(3,11)=clhs163;
lhs(3,12)=-DN(0,0)*clhs164 - clhs103;
lhs(3,13)=-DN(0,1)*clhs164 - clhs138;
lhs(3,14)=-DN(0,2)*clhs164 - clhs153;
lhs(3,15)=clhs165;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs169;
lhs(4,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + DN(1,2)*clhs25;
lhs(4,2)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + DN(1,2)*clhs30;
lhs(4,3)=k*(-DN(0,0)*clhs171 + DN(1,0)*N[0] - clhs170*clhs58);
lhs(4,4)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + DN(1,2)*clhs39 + clhs173;
lhs(4,5)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs51;
lhs(4,6)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + DN(1,2)*clhs56;
lhs(4,7)=DN(1,0)*clhs174;
lhs(4,8)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + DN(1,2)*clhs63 + clhs177;
lhs(4,9)=DN(1,0)*clhs70 + DN(1,1)*clhs72 + DN(1,2)*clhs75;
lhs(4,10)=DN(1,0)*clhs76 + DN(1,1)*clhs78 + DN(1,2)*clhs80;
lhs(4,11)=k*(DN(1,0)*N[2] - DN(2,0)*clhs171 - clhs178*clhs58);
lhs(4,12)=DN(1,0)*clhs82 + DN(1,1)*clhs84 + DN(1,2)*clhs86 + clhs181;
lhs(4,13)=DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs97;
lhs(4,14)=DN(1,0)*clhs98 + DN(1,1)*clhs100 + DN(1,2)*clhs102;
lhs(4,15)=k*(DN(1,0)*N[3] - DN(3,0)*clhs171 - clhs182*clhs58);
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs104 + DN(1,2)*clhs105;
lhs(5,1)=DN(1,0)*clhs22 + DN(1,1)*clhs106 + DN(1,2)*clhs108 + clhs169;
lhs(5,2)=DN(1,0)*clhs28 + DN(1,1)*clhs109 + DN(1,2)*clhs111;
lhs(5,3)=k*(-DN(0,1)*clhs171 + DN(1,1)*N[0] - clhs183*clhs58);
lhs(5,4)=DN(1,0)*clhs37 + DN(1,1)*clhs112 + DN(1,2)*clhs113;
lhs(5,5)=DN(1,0)*clhs48 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs173;
lhs(5,6)=DN(1,0)*clhs54 + DN(1,1)*clhs117 + DN(1,2)*clhs119;
lhs(5,7)=DN(1,1)*clhs174;
lhs(5,8)=DN(1,0)*clhs61 + DN(1,1)*clhs121 + DN(1,2)*clhs122;
lhs(5,9)=DN(1,0)*clhs72 + DN(1,1)*clhs123 + DN(1,2)*clhs125 + clhs177;
lhs(5,10)=DN(1,0)*clhs78 + DN(1,1)*clhs126 + DN(1,2)*clhs128;
lhs(5,11)=k*(DN(1,1)*N[2] - DN(2,1)*clhs171 - clhs184*clhs58);
lhs(5,12)=DN(1,0)*clhs84 + DN(1,1)*clhs130 + DN(1,2)*clhs131;
lhs(5,13)=DN(1,0)*clhs94 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs181;
lhs(5,14)=DN(1,0)*clhs100 + DN(1,1)*clhs135 + DN(1,2)*clhs137;
lhs(5,15)=k*(DN(1,1)*N[3] - DN(3,1)*clhs171 - clhs185*clhs58);
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs105 + DN(1,2)*clhs139;
lhs(6,1)=DN(1,0)*clhs25 + DN(1,1)*clhs108 + DN(1,2)*clhs140;
lhs(6,2)=DN(1,0)*clhs30 + DN(1,1)*clhs111 + DN(1,2)*clhs141 + clhs169;
lhs(6,3)=k*(-DN(0,2)*clhs171 + DN(1,2)*N[0] - clhs186*clhs58);
lhs(6,4)=DN(1,0)*clhs39 + DN(1,1)*clhs113 + DN(1,2)*clhs142;
lhs(6,5)=DN(1,0)*clhs51 + DN(1,1)*clhs116 + DN(1,2)*clhs143;
lhs(6,6)=DN(1,0)*clhs56 + DN(1,1)*clhs119 + DN(1,2)*clhs144 + clhs173;
lhs(6,7)=DN(1,2)*clhs174;
lhs(6,8)=DN(1,0)*clhs63 + DN(1,1)*clhs122 + DN(1,2)*clhs146;
lhs(6,9)=DN(1,0)*clhs75 + DN(1,1)*clhs125 + DN(1,2)*clhs147;
lhs(6,10)=DN(1,0)*clhs80 + DN(1,1)*clhs128 + DN(1,2)*clhs148 + clhs177;
lhs(6,11)=k*(DN(1,2)*N[2] - DN(2,2)*clhs171 - clhs187*clhs58);
lhs(6,12)=DN(1,0)*clhs86 + DN(1,1)*clhs131 + DN(1,2)*clhs150;
lhs(6,13)=DN(1,0)*clhs97 + DN(1,1)*clhs134 + DN(1,2)*clhs151;
lhs(6,14)=DN(1,0)*clhs102 + DN(1,1)*clhs137 + DN(1,2)*clhs152 + clhs181;
lhs(6,15)=k*(DN(1,2)*N[3] - DN(3,2)*clhs171 - clhs188*clhs58);
lhs(7,0)=-DN(1,0)*clhs154 - clhs170;
lhs(7,1)=-DN(1,1)*clhs154 - clhs183;
lhs(7,2)=-DN(1,2)*clhs154 - clhs186;
lhs(7,3)=clhs161;
lhs(7,4)=-DN(1,0)*clhs189;
lhs(7,5)=-DN(1,1)*clhs189;
lhs(7,6)=-DN(1,2)*clhs189;
lhs(7,7)=pow(DN(1,0), 2)*clhs156 + pow(DN(1,1), 2)*clhs156 + pow(DN(1,2), 2)*clhs156 + clhs172;
lhs(7,8)=-DN(1,0)*clhs162 - clhs178;
lhs(7,9)=-DN(1,1)*clhs162 - clhs184;
lhs(7,10)=-DN(1,2)*clhs162 - clhs187;
lhs(7,11)=clhs193;
lhs(7,12)=-DN(1,0)*clhs164 - clhs182;
lhs(7,13)=-DN(1,1)*clhs164 - clhs185;
lhs(7,14)=-DN(1,2)*clhs164 - clhs188;
lhs(7,15)=clhs194;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs198;
lhs(8,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + DN(2,2)*clhs25;
lhs(8,2)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + DN(2,2)*clhs30;
lhs(8,3)=k*(-DN(0,0)*clhs200 + DN(2,0)*N[0] - clhs199*clhs58);
lhs(8,4)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + DN(2,2)*clhs39 + clhs201;
lhs(8,5)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs51;
lhs(8,6)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + DN(2,2)*clhs56;
lhs(8,7)=k*(-DN(1,0)*clhs200 + DN(2,0)*N[1] - clhs202*clhs58);
lhs(8,8)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + DN(2,2)*clhs63 + clhs204;
lhs(8,9)=DN(2,0)*clhs70 + DN(2,1)*clhs72 + DN(2,2)*clhs75;
lhs(8,10)=DN(2,0)*clhs76 + DN(2,1)*clhs78 + DN(2,2)*clhs80;
lhs(8,11)=DN(2,0)*clhs205;
lhs(8,12)=DN(2,0)*clhs82 + DN(2,1)*clhs84 + DN(2,2)*clhs86 + clhs208;
lhs(8,13)=DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs97;
lhs(8,14)=DN(2,0)*clhs98 + DN(2,1)*clhs100 + DN(2,2)*clhs102;
lhs(8,15)=k*(DN(2,0)*N[3] - DN(3,0)*clhs200 - clhs209*clhs58);
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs104 + DN(2,2)*clhs105;
lhs(9,1)=DN(2,0)*clhs22 + DN(2,1)*clhs106 + DN(2,2)*clhs108 + clhs198;
lhs(9,2)=DN(2,0)*clhs28 + DN(2,1)*clhs109 + DN(2,2)*clhs111;
lhs(9,3)=k*(-DN(0,1)*clhs200 + DN(2,1)*N[0] - clhs210*clhs58);
lhs(9,4)=DN(2,0)*clhs37 + DN(2,1)*clhs112 + DN(2,2)*clhs113;
lhs(9,5)=DN(2,0)*clhs48 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs201;
lhs(9,6)=DN(2,0)*clhs54 + DN(2,1)*clhs117 + DN(2,2)*clhs119;
lhs(9,7)=k*(-DN(1,1)*clhs200 + DN(2,1)*N[1] - clhs211*clhs58);
lhs(9,8)=DN(2,0)*clhs61 + DN(2,1)*clhs121 + DN(2,2)*clhs122;
lhs(9,9)=DN(2,0)*clhs72 + DN(2,1)*clhs123 + DN(2,2)*clhs125 + clhs204;
lhs(9,10)=DN(2,0)*clhs78 + DN(2,1)*clhs126 + DN(2,2)*clhs128;
lhs(9,11)=DN(2,1)*clhs205;
lhs(9,12)=DN(2,0)*clhs84 + DN(2,1)*clhs130 + DN(2,2)*clhs131;
lhs(9,13)=DN(2,0)*clhs94 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs208;
lhs(9,14)=DN(2,0)*clhs100 + DN(2,1)*clhs135 + DN(2,2)*clhs137;
lhs(9,15)=k*(DN(2,1)*N[3] - DN(3,1)*clhs200 - clhs212*clhs58);
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs105 + DN(2,2)*clhs139;
lhs(10,1)=DN(2,0)*clhs25 + DN(2,1)*clhs108 + DN(2,2)*clhs140;
lhs(10,2)=DN(2,0)*clhs30 + DN(2,1)*clhs111 + DN(2,2)*clhs141 + clhs198;
lhs(10,3)=k*(-DN(0,2)*clhs200 + DN(2,2)*N[0] - clhs213*clhs58);
lhs(10,4)=DN(2,0)*clhs39 + DN(2,1)*clhs113 + DN(2,2)*clhs142;
lhs(10,5)=DN(2,0)*clhs51 + DN(2,1)*clhs116 + DN(2,2)*clhs143;
lhs(10,6)=DN(2,0)*clhs56 + DN(2,1)*clhs119 + DN(2,2)*clhs144 + clhs201;
lhs(10,7)=k*(-DN(1,2)*clhs200 + DN(2,2)*N[1] - clhs214*clhs58);
lhs(10,8)=DN(2,0)*clhs63 + DN(2,1)*clhs122 + DN(2,2)*clhs146;
lhs(10,9)=DN(2,0)*clhs75 + DN(2,1)*clhs125 + DN(2,2)*clhs147;
lhs(10,10)=DN(2,0)*clhs80 + DN(2,1)*clhs128 + DN(2,2)*clhs148 + clhs204;
lhs(10,11)=DN(2,2)*clhs205;
lhs(10,12)=DN(2,0)*clhs86 + DN(2,1)*clhs131 + DN(2,2)*clhs150;
lhs(10,13)=DN(2,0)*clhs97 + DN(2,1)*clhs134 + DN(2,2)*clhs151;
lhs(10,14)=DN(2,0)*clhs102 + DN(2,1)*clhs137 + DN(2,2)*clhs152 + clhs208;
lhs(10,15)=k*(DN(2,2)*N[3] - DN(3,2)*clhs200 - clhs215*clhs58);
lhs(11,0)=-DN(2,0)*clhs154 - clhs199;
lhs(11,1)=-DN(2,1)*clhs154 - clhs210;
lhs(11,2)=-DN(2,2)*clhs154 - clhs213;
lhs(11,3)=clhs163;
lhs(11,4)=-DN(2,0)*clhs157 - clhs202;
lhs(11,5)=-DN(2,1)*clhs157 - clhs211;
lhs(11,6)=-DN(2,2)*clhs157 - clhs214;
lhs(11,7)=clhs193;
lhs(11,8)=-DN(2,0)*clhs216;
lhs(11,9)=-DN(2,1)*clhs216;
lhs(11,10)=-DN(2,2)*clhs216;
lhs(11,11)=pow(DN(2,0), 2)*clhs156 + pow(DN(2,1), 2)*clhs156 + pow(DN(2,2), 2)*clhs156 + clhs203;
lhs(11,12)=-DN(2,0)*clhs164 - clhs209;
lhs(11,13)=-DN(2,1)*clhs164 - clhs212;
lhs(11,14)=-DN(2,2)*clhs164 - clhs215;
lhs(11,15)=clhs217;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs221;
lhs(12,1)=DN(3,0)*clhs20 + DN(3,1)*clhs22 + DN(3,2)*clhs25;
lhs(12,2)=DN(3,0)*clhs26 + DN(3,1)*clhs28 + DN(3,2)*clhs30;
lhs(12,3)=k*(-DN(0,0)*clhs223 + DN(3,0)*N[0] - clhs222*clhs58);
lhs(12,4)=DN(3,0)*clhs35 + DN(3,1)*clhs37 + DN(3,2)*clhs39 + clhs224;
lhs(12,5)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs51;
lhs(12,6)=DN(3,0)*clhs52 + DN(3,1)*clhs54 + DN(3,2)*clhs56;
lhs(12,7)=k*(-DN(1,0)*clhs223 + DN(3,0)*N[1] - clhs225*clhs58);
lhs(12,8)=DN(3,0)*clhs59 + DN(3,1)*clhs61 + DN(3,2)*clhs63 + clhs226;
lhs(12,9)=DN(3,0)*clhs70 + DN(3,1)*clhs72 + DN(3,2)*clhs75;
lhs(12,10)=DN(3,0)*clhs76 + DN(3,1)*clhs78 + DN(3,2)*clhs80;
lhs(12,11)=k*(-DN(2,0)*clhs223 + DN(3,0)*N[2] - clhs227*clhs58);
lhs(12,12)=DN(3,0)*clhs82 + DN(3,1)*clhs84 + DN(3,2)*clhs86 + clhs229;
lhs(12,13)=DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs97;
lhs(12,14)=DN(3,0)*clhs98 + DN(3,1)*clhs100 + DN(3,2)*clhs102;
lhs(12,15)=DN(3,0)*clhs230;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs104 + DN(3,2)*clhs105;
lhs(13,1)=DN(3,0)*clhs22 + DN(3,1)*clhs106 + DN(3,2)*clhs108 + clhs221;
lhs(13,2)=DN(3,0)*clhs28 + DN(3,1)*clhs109 + DN(3,2)*clhs111;
lhs(13,3)=k*(-DN(0,1)*clhs223 + DN(3,1)*N[0] - clhs231*clhs58);
lhs(13,4)=DN(3,0)*clhs37 + DN(3,1)*clhs112 + DN(3,2)*clhs113;
lhs(13,5)=DN(3,0)*clhs48 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs224;
lhs(13,6)=DN(3,0)*clhs54 + DN(3,1)*clhs117 + DN(3,2)*clhs119;
lhs(13,7)=k*(-DN(1,1)*clhs223 + DN(3,1)*N[1] - clhs232*clhs58);
lhs(13,8)=DN(3,0)*clhs61 + DN(3,1)*clhs121 + DN(3,2)*clhs122;
lhs(13,9)=DN(3,0)*clhs72 + DN(3,1)*clhs123 + DN(3,2)*clhs125 + clhs226;
lhs(13,10)=DN(3,0)*clhs78 + DN(3,1)*clhs126 + DN(3,2)*clhs128;
lhs(13,11)=k*(-DN(2,1)*clhs223 + DN(3,1)*N[2] - clhs233*clhs58);
lhs(13,12)=DN(3,0)*clhs84 + DN(3,1)*clhs130 + DN(3,2)*clhs131;
lhs(13,13)=DN(3,0)*clhs94 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs229;
lhs(13,14)=DN(3,0)*clhs100 + DN(3,1)*clhs135 + DN(3,2)*clhs137;
lhs(13,15)=DN(3,1)*clhs230;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs105 + DN(3,2)*clhs139;
lhs(14,1)=DN(3,0)*clhs25 + DN(3,1)*clhs108 + DN(3,2)*clhs140;
lhs(14,2)=DN(3,0)*clhs30 + DN(3,1)*clhs111 + DN(3,2)*clhs141 + clhs221;
lhs(14,3)=k*(-DN(0,2)*clhs223 + DN(3,2)*N[0] - clhs234*clhs58);
lhs(14,4)=DN(3,0)*clhs39 + DN(3,1)*clhs113 + DN(3,2)*clhs142;
lhs(14,5)=DN(3,0)*clhs51 + DN(3,1)*clhs116 + DN(3,2)*clhs143;
lhs(14,6)=DN(3,0)*clhs56 + DN(3,1)*clhs119 + DN(3,2)*clhs144 + clhs224;
lhs(14,7)=k*(-DN(1,2)*clhs223 + DN(3,2)*N[1] - clhs235*clhs58);
lhs(14,8)=DN(3,0)*clhs63 + DN(3,1)*clhs122 + DN(3,2)*clhs146;
lhs(14,9)=DN(3,0)*clhs75 + DN(3,1)*clhs125 + DN(3,2)*clhs147;
lhs(14,10)=DN(3,0)*clhs80 + DN(3,1)*clhs128 + DN(3,2)*clhs148 + clhs226;
lhs(14,11)=k*(-DN(2,2)*clhs223 + DN(3,2)*N[2] - clhs236*clhs58);
lhs(14,12)=DN(3,0)*clhs86 + DN(3,1)*clhs131 + DN(3,2)*clhs150;
lhs(14,13)=DN(3,0)*clhs97 + DN(3,1)*clhs134 + DN(3,2)*clhs151;
lhs(14,14)=DN(3,0)*clhs102 + DN(3,1)*clhs137 + DN(3,2)*clhs152 + clhs229;
lhs(14,15)=DN(3,2)*clhs230;
lhs(15,0)=-DN(3,0)*clhs154 - clhs222;
lhs(15,1)=-DN(3,1)*clhs154 - clhs231;
lhs(15,2)=-DN(3,2)*clhs154 - clhs234;
lhs(15,3)=clhs165;
lhs(15,4)=-DN(3,0)*clhs157 - clhs225;
lhs(15,5)=-DN(3,1)*clhs157 - clhs232;
lhs(15,6)=-DN(3,2)*clhs157 - clhs235;
lhs(15,7)=clhs194;
lhs(15,8)=-DN(3,0)*clhs162 - clhs227;
lhs(15,9)=-DN(3,1)*clhs162 - clhs233;
lhs(15,10)=-DN(3,2)*clhs162 - clhs236;
lhs(15,11)=clhs217;
lhs(15,12)=-DN(3,0)*clhs237;
lhs(15,13)=-DN(3,1)*clhs237;
lhs(15,14)=-DN(3,2)*clhs237;
lhs(15,15)=pow(DN(3,0), 2)*clhs156 + pow(DN(3,1), 2)*clhs156 + pow(DN(3,2), 2)*clhs156 + clhs228;


    // Add intermediate results to local system
    noalias(rLHS) += lhs * rData.Weight;
}

template <>
void WeaklyCompressibleNavierStokesVolumetricStrain<WeaklyCompressibleNavierStokesVolumetricStrainData<2,3>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesVolumetricStrainData<2,3>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double k = rData.BulkModulus;
    const double mu = rData.EffectiveViscosity;

    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;
    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3>& v = rData.Velocity;
    const BoundedMatrix<double,2,3>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,3>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = v - vmesh;
    const BoundedMatrix<double,2,3>& f = rData.BodyForce;
    const array_1d<double,3>& eps_vol = rData.Pressure;
    const array_1d<double,3>& eps_vol_n = rData.Pressure_OldStep1;
    const array_1d<double,3>& eps_vol_nn = rData.Pressure_OldStep2;
    const array_1d<double,3>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = k*(N[0]*eps_vol[0] + N[1]*eps_vol[1] + N[2]*eps_vol[2]);
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs8 = 1.0/(rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs9 = crhs1*rho - crhs2 - crhs6 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2]);
const double crhs10 = DN(0,0)*crhs4 + DN(0,1)*crhs5;
const double crhs11 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs12 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs14 = rho*(crhs13*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs15 = crhs11*rho - crhs12 - crhs14 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2]);
const double crhs16 = crhs13 + crhs3;
const double crhs17 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]) + N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]) + N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]);
const double crhs18 = 1.0*crhs8;
const double crhs19 = crhs18*crhs9;
const double crhs20 = crhs15*crhs18;
const double crhs21 = DN(1,0)*crhs4 + DN(1,1)*crhs5;
const double crhs22 = DN(2,0)*crhs4 + DN(2,1)*crhs5;
rhs[0]=-DN(0,0)*crhs0 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1*rho - N[0]*crhs2 - N[0]*crhs6 + 1.0*N[0]*crhs7*crhs8*crhs9*rho + 1.0*crhs10*crhs8*crhs9*rho;
rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs0 - DN(0,1)*stress[1] + N[0]*crhs11*rho - N[0]*crhs12 - N[0]*crhs14 + 1.0*N[0]*crhs15*crhs7*crhs8*rho + 1.0*crhs10*crhs15*crhs8*rho;
rhs[2]=-DN(0,0)*crhs19 - DN(0,1)*crhs20 + N[0]*crhs16 - N[0]*crhs17;
rhs[3]=-DN(1,0)*crhs0 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1*rho - N[1]*crhs2 - N[1]*crhs6 + 1.0*N[1]*crhs7*crhs8*crhs9*rho + 1.0*crhs21*crhs8*crhs9*rho;
rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs0 - DN(1,1)*stress[1] + N[1]*crhs11*rho - N[1]*crhs12 - N[1]*crhs14 + 1.0*N[1]*crhs15*crhs7*crhs8*rho + 1.0*crhs15*crhs21*crhs8*rho;
rhs[5]=-DN(1,0)*crhs19 - DN(1,1)*crhs20 + N[1]*crhs16 - N[1]*crhs17;
rhs[6]=-DN(2,0)*crhs0 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1*rho - N[2]*crhs2 - N[2]*crhs6 + 1.0*N[2]*crhs7*crhs8*crhs9*rho + 1.0*crhs22*crhs8*crhs9*rho;
rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs0 - DN(2,1)*stress[1] + N[2]*crhs11*rho - N[2]*crhs12 - N[2]*crhs14 + 1.0*N[2]*crhs15*crhs7*crhs8*rho + 1.0*crhs15*crhs22*crhs8*rho;
rhs[8]=-DN(2,0)*crhs19 - DN(2,1)*crhs20 + N[2]*crhs16 - N[2]*crhs17;


    noalias(rRHS) += rData.Weight * rhs;
}

template <>
void WeaklyCompressibleNavierStokesVolumetricStrain<WeaklyCompressibleNavierStokesVolumetricStrainData<3,4>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesVolumetricStrainData<3,4>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double k = rData.BulkModulus;
    const double mu = rData.EffectiveViscosity;

    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;
    const double dt = rData.DeltaTime;
    const double h = rData.ElementSize;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,3,4>& v = rData.Velocity;
    const BoundedMatrix<double,3,4>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,3,4>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = v - vmesh;
    const BoundedMatrix<double,3,4>& f = rData.BodyForce;
    const array_1d<double,4>& eps_vol = rData.Pressure;
    const array_1d<double,4>& eps_vol_n = rData.Pressure_OldStep1;
    const array_1d<double,4>& eps_vol_nn = rData.Pressure_OldStep2;
    const array_1d<double,6>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //TODO: Optimize this to directly add to the rRightHandSideVector
    auto& rhs = rData.rhs;

    const double crhs0 = k*(N[0]*eps_vol[0] + N[1]*eps_vol[1] + N[2]*eps_vol[2] + N[3]*eps_vol[3]);
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs7 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs8 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs9 = 1.0/(rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2) + pow(crhs6, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs10 = crhs1*rho - crhs2 - crhs7 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2] + DN(3,0)*eps_vol[3]);
const double crhs11 = DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6;
const double crhs12 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs13 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs14 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs15 = rho*(crhs14*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs16 = crhs12*rho - crhs13 - crhs15 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2] + DN(3,1)*eps_vol[3]);
const double crhs17 = N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2);
const double crhs18 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs19 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs20 = rho*(crhs19*crhs6 + crhs4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs21 = crhs17*rho - crhs18 - crhs20 + k*(DN(0,2)*eps_vol[0] + DN(1,2)*eps_vol[1] + DN(2,2)*eps_vol[2] + DN(3,2)*eps_vol[3]);
const double crhs22 = crhs14 + crhs19 + crhs3;
const double crhs23 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]) + N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]) + N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]) + N[3]*(bdf0*eps_vol[3] + bdf1*eps_vol_n[3] + bdf2*eps_vol_nn[3]);
const double crhs24 = 1.0*crhs9;
const double crhs25 = crhs10*crhs24;
const double crhs26 = crhs16*crhs24;
const double crhs27 = crhs21*crhs24;
const double crhs28 = DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6;
const double crhs29 = DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6;
const double crhs30 = DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6;
rhs[0]=-DN(0,0)*crhs0 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1*rho + 1.0*N[0]*crhs10*crhs8*crhs9*rho - N[0]*crhs2 - N[0]*crhs7 + 1.0*crhs10*crhs11*crhs9*rho;
rhs[1]=-DN(0,0)*stress[3] - DN(0,1)*crhs0 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs12*rho - N[0]*crhs13 - N[0]*crhs15 + 1.0*N[0]*crhs16*crhs8*crhs9*rho + 1.0*crhs11*crhs16*crhs9*rho;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] - DN(0,2)*crhs0 - DN(0,2)*stress[2] + N[0]*crhs17*rho - N[0]*crhs18 - N[0]*crhs20 + 1.0*N[0]*crhs21*crhs8*crhs9*rho + 1.0*crhs11*crhs21*crhs9*rho;
rhs[3]=-DN(0,0)*crhs25 - DN(0,1)*crhs26 - DN(0,2)*crhs27 + N[0]*crhs22 - N[0]*crhs23;
rhs[4]=-DN(1,0)*crhs0 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1*rho + 1.0*N[1]*crhs10*crhs8*crhs9*rho - N[1]*crhs2 - N[1]*crhs7 + 1.0*crhs10*crhs28*crhs9*rho;
rhs[5]=-DN(1,0)*stress[3] - DN(1,1)*crhs0 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs12*rho - N[1]*crhs13 - N[1]*crhs15 + 1.0*N[1]*crhs16*crhs8*crhs9*rho + 1.0*crhs16*crhs28*crhs9*rho;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] - DN(1,2)*crhs0 - DN(1,2)*stress[2] + N[1]*crhs17*rho - N[1]*crhs18 - N[1]*crhs20 + 1.0*N[1]*crhs21*crhs8*crhs9*rho + 1.0*crhs21*crhs28*crhs9*rho;
rhs[7]=-DN(1,0)*crhs25 - DN(1,1)*crhs26 - DN(1,2)*crhs27 + N[1]*crhs22 - N[1]*crhs23;
rhs[8]=-DN(2,0)*crhs0 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1*rho + 1.0*N[2]*crhs10*crhs8*crhs9*rho - N[2]*crhs2 - N[2]*crhs7 + 1.0*crhs10*crhs29*crhs9*rho;
rhs[9]=-DN(2,0)*stress[3] - DN(2,1)*crhs0 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs12*rho - N[2]*crhs13 - N[2]*crhs15 + 1.0*N[2]*crhs16*crhs8*crhs9*rho + 1.0*crhs16*crhs29*crhs9*rho;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] - DN(2,2)*crhs0 - DN(2,2)*stress[2] + N[2]*crhs17*rho - N[2]*crhs18 - N[2]*crhs20 + 1.0*N[2]*crhs21*crhs8*crhs9*rho + 1.0*crhs21*crhs29*crhs9*rho;
rhs[11]=-DN(2,0)*crhs25 - DN(2,1)*crhs26 - DN(2,2)*crhs27 + N[2]*crhs22 - N[2]*crhs23;
rhs[12]=-DN(3,0)*crhs0 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1*rho + 1.0*N[3]*crhs10*crhs8*crhs9*rho - N[3]*crhs2 - N[3]*crhs7 + 1.0*crhs10*crhs30*crhs9*rho;
rhs[13]=-DN(3,0)*stress[3] - DN(3,1)*crhs0 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs12*rho - N[3]*crhs13 - N[3]*crhs15 + 1.0*N[3]*crhs16*crhs8*crhs9*rho + 1.0*crhs16*crhs30*crhs9*rho;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] - DN(3,2)*crhs0 - DN(3,2)*stress[2] + N[3]*crhs17*rho - N[3]*crhs18 - N[3]*crhs20 + 1.0*N[3]*crhs21*crhs8*crhs9*rho + 1.0*crhs21*crhs30*crhs9*rho;
rhs[15]=-DN(3,0)*crhs25 - DN(3,1)*crhs26 - DN(3,2)*crhs27 + N[3]*crhs22 - N[3]*crhs23;


    noalias(rRHS) += rData.Weight * rhs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void WeaklyCompressibleNavierStokesVolumetricStrain<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class WeaklyCompressibleNavierStokesVolumetricStrain< WeaklyCompressibleNavierStokesVolumetricStrainData<2,3> >;
template class WeaklyCompressibleNavierStokesVolumetricStrain< WeaklyCompressibleNavierStokesVolumetricStrainData<3,4> >;

}