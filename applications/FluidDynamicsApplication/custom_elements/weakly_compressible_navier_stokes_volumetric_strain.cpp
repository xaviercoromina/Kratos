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
const double clhs3 = pow(DN(0,0), 2);
const double clhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 = rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 = clhs6*h/stab_c1 + mu;
const double clhs8 = clhs7*k;
const double clhs9 = DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs10 = N[0]*rho;
const double clhs11 = pow(N[0], 2);
const double clhs12 = N[0]*bdf0;
const double clhs13 = clhs12 + clhs9;
const double clhs14 = 1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 = clhs14*pow(rho, 2);
const double clhs16 = clhs15*clhs9;
const double clhs17 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs18 = clhs15*clhs17;
const double clhs19 = N[0]*clhs18;
const double clhs20 = bdf0*clhs11*rho + clhs10*clhs9 + clhs13*clhs16 + clhs13*clhs19;
const double clhs21 = C(0,1)*DN(0,1) + clhs1;
const double clhs22 = C(1,2)*DN(0,1);
const double clhs23 = C(2,2)*DN(0,0) + clhs22;
const double clhs24 = DN(0,0)*clhs8;
const double clhs25 = DN(0,1)*clhs24;
const double clhs26 = clhs14*clhs17;
const double clhs27 = clhs14*rho;
const double clhs28 = clhs27*clhs9;
const double clhs29 = k*(N[0]*bdf0*clhs7 - N[0] - clhs10*clhs26 - clhs28);
const double clhs30 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs31 = C(0,2)*DN(1,0);
const double clhs32 = C(2,2)*DN(1,1) + clhs31;
const double clhs33 = N[1]*clhs12*rho;
const double clhs34 = DN(1,0)*clhs24 + clhs33;
const double clhs35 = DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs36 = N[1]*bdf0;
const double clhs37 = clhs35 + clhs36;
const double clhs38 = clhs10*clhs35 + clhs16*clhs37 + clhs19*clhs37;
const double clhs39 = C(0,1)*DN(1,1) + clhs31;
const double clhs40 = C(1,2)*DN(1,1);
const double clhs41 = C(2,2)*DN(1,0) + clhs40;
const double clhs42 = DN(1,1)*clhs24;
const double clhs43 = DN(0,0)*N[1];
const double clhs44 = DN(1,0)*N[0];
const double clhs45 = clhs17*clhs27;
const double clhs46 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs47 = C(0,2)*DN(2,0);
const double clhs48 = C(2,2)*DN(2,1) + clhs47;
const double clhs49 = N[2]*clhs12*rho;
const double clhs50 = DN(2,0)*clhs24 + clhs49;
const double clhs51 = DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs52 = N[2]*bdf0 + clhs51;
const double clhs53 = clhs10*clhs51 + clhs16*clhs52 + clhs19*clhs52;
const double clhs54 = C(0,1)*DN(2,1) + clhs47;
const double clhs55 = C(1,2)*DN(2,1);
const double clhs56 = C(2,2)*DN(2,0) + clhs55;
const double clhs57 = DN(2,1)*clhs24;
const double clhs58 = DN(0,0)*N[2];
const double clhs59 = DN(2,0)*N[0];
const double clhs60 = C(0,1)*DN(0,0) + clhs22;
const double clhs61 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs62 = pow(DN(0,1), 2);
const double clhs63 = C(0,1)*DN(1,0) + clhs40;
const double clhs64 = DN(0,1)*clhs8;
const double clhs65 = DN(1,0)*clhs64;
const double clhs66 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs67 = DN(1,1)*clhs64 + clhs33;
const double clhs68 = DN(0,1)*N[1];
const double clhs69 = DN(1,1)*N[0];
const double clhs70 = C(0,1)*DN(2,0) + clhs55;
const double clhs71 = DN(2,0)*clhs64;
const double clhs72 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs73 = DN(2,1)*clhs64 + clhs49;
const double clhs74 = DN(0,1)*N[2];
const double clhs75 = DN(2,1)*N[0];
const double clhs76 = clhs13*clhs27;
const double clhs77 = N[0] + clhs76;
const double clhs78 = clhs14*k;
const double clhs79 = clhs27*clhs37;
const double clhs80 = DN(0,0)*clhs78;
const double clhs81 = DN(0,1)*clhs78;
const double clhs82 = -DN(1,0)*clhs80 - DN(1,1)*clhs81 + N[0]*N[1]*bdf0;
const double clhs83 = clhs27*clhs52;
const double clhs84 = -DN(2,0)*clhs80 - DN(2,1)*clhs81 + N[0]*N[2]*bdf0;
const double clhs85 = N[1]*rho;
const double clhs86 = clhs15*clhs35;
const double clhs87 = N[1]*clhs18;
const double clhs88 = clhs13*clhs86 + clhs13*clhs87 + clhs85*clhs9;
const double clhs89 = clhs27*clhs35;
const double clhs90 = pow(DN(1,0), 2);
const double clhs91 = pow(N[1], 2);
const double clhs92 = bdf0*clhs91*rho + clhs35*clhs85 + clhs37*clhs86 + clhs37*clhs87;
const double clhs93 = DN(1,0)*clhs8;
const double clhs94 = DN(1,1)*clhs93;
const double clhs95 = k*(N[1]*bdf0*clhs7 - N[1] - clhs26*clhs85 - clhs89);
const double clhs96 = N[2]*clhs36*rho;
const double clhs97 = DN(2,0)*clhs93 + clhs96;
const double clhs98 = clhs51*clhs85 + clhs52*clhs86 + clhs52*clhs87;
const double clhs99 = DN(2,1)*clhs93;
const double clhs100 = DN(1,0)*N[2];
const double clhs101 = DN(2,0)*N[1];
const double clhs102 = pow(DN(1,1), 2);
const double clhs103 = DN(1,1)*clhs8;
const double clhs104 = DN(2,0)*clhs103;
const double clhs105 = DN(2,1)*clhs103 + clhs96;
const double clhs106 = DN(1,1)*N[2];
const double clhs107 = DN(2,1)*N[1];
const double clhs108 = N[1] + clhs79;
const double clhs109 = -DN(1,0)*DN(2,0)*clhs78 - DN(1,1)*DN(2,1)*clhs78 + N[1]*N[2]*bdf0;
const double clhs110 = N[2]*rho;
const double clhs111 = clhs15*clhs51;
const double clhs112 = N[2]*clhs18;
const double clhs113 = clhs110*clhs9 + clhs111*clhs13 + clhs112*clhs13;
const double clhs114 = clhs27*clhs51;
const double clhs115 = clhs110*clhs35 + clhs111*clhs37 + clhs112*clhs37;
const double clhs116 = pow(DN(2,0), 2);
const double clhs117 = pow(N[2], 2);
const double clhs118 = bdf0*clhs117*rho + clhs110*clhs51 + clhs111*clhs52 + clhs112*clhs52;
const double clhs119 = DN(2,0)*DN(2,1)*clhs8;
const double clhs120 = k*(N[2]*bdf0*clhs7 - N[2] - clhs110*clhs26 - clhs114);
const double clhs121 = pow(DN(2,1), 2);
const double clhs122 = N[2] + clhs83;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs20 + clhs3*clhs8;
lhs(0,1)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs25;
lhs(0,2)=DN(0,0)*clhs29;
lhs(0,3)=DN(0,0)*clhs30 + DN(0,1)*clhs32 + clhs34 + clhs38;
lhs(0,4)=DN(0,0)*clhs39 + DN(0,1)*clhs41 + clhs42;
lhs(0,5)=k*(DN(0,0)*N[1]*bdf0*clhs7 - DN(1,0)*clhs28 - clhs43 - clhs44*clhs45);
lhs(0,6)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + clhs50 + clhs53;
lhs(0,7)=DN(0,0)*clhs54 + DN(0,1)*clhs56 + clhs57;
lhs(0,8)=k*(DN(0,0)*N[2]*bdf0*clhs7 - DN(2,0)*clhs28 - clhs45*clhs59 - clhs58);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs60 + clhs25;
lhs(1,1)=DN(0,0)*clhs23 + DN(0,1)*clhs61 + clhs20 + clhs62*clhs8;
lhs(1,2)=DN(0,1)*clhs29;
lhs(1,3)=DN(0,0)*clhs32 + DN(0,1)*clhs63 + clhs65;
lhs(1,4)=DN(0,0)*clhs41 + DN(0,1)*clhs66 + clhs38 + clhs67;
lhs(1,5)=k*(DN(0,1)*N[1]*bdf0*clhs7 - DN(1,1)*clhs28 - clhs45*clhs69 - clhs68);
lhs(1,6)=DN(0,0)*clhs48 + DN(0,1)*clhs70 + clhs71;
lhs(1,7)=DN(0,0)*clhs56 + DN(0,1)*clhs72 + clhs53 + clhs73;
lhs(1,8)=k*(DN(0,1)*N[2]*bdf0*clhs7 - DN(2,1)*clhs28 - clhs45*clhs75 - clhs74);
lhs(2,0)=DN(0,0)*clhs77;
lhs(2,1)=DN(0,1)*clhs77;
lhs(2,2)=bdf0*clhs11 - clhs3*clhs78 - clhs62*clhs78;
lhs(2,3)=DN(0,0)*clhs79 + clhs44;
lhs(2,4)=DN(0,1)*clhs79 + clhs69;
lhs(2,5)=clhs82;
lhs(2,6)=DN(0,0)*clhs83 + clhs59;
lhs(2,7)=DN(0,1)*clhs83 + clhs75;
lhs(2,8)=clhs84;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs34 + clhs88;
lhs(3,1)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs65;
lhs(3,2)=k*(-DN(0,0)*clhs89 + DN(1,0)*N[0]*bdf0*clhs7 - clhs43*clhs45 - clhs44);
lhs(3,3)=DN(1,0)*clhs30 + DN(1,1)*clhs32 + clhs8*clhs90 + clhs92;
lhs(3,4)=DN(1,0)*clhs39 + DN(1,1)*clhs41 + clhs94;
lhs(3,5)=DN(1,0)*clhs95;
lhs(3,6)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + clhs97 + clhs98;
lhs(3,7)=DN(1,0)*clhs54 + DN(1,1)*clhs56 + clhs99;
lhs(3,8)=k*(DN(1,0)*N[2]*bdf0*clhs7 - DN(2,0)*clhs89 - clhs100 - clhs101*clhs45);
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs60 + clhs42;
lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs61 + clhs67 + clhs88;
lhs(4,2)=k*(-DN(0,1)*clhs89 + DN(1,1)*N[0]*bdf0*clhs7 - clhs45*clhs68 - clhs69);
lhs(4,3)=DN(1,0)*clhs32 + DN(1,1)*clhs63 + clhs94;
lhs(4,4)=DN(1,0)*clhs41 + DN(1,1)*clhs66 + clhs102*clhs8 + clhs92;
lhs(4,5)=DN(1,1)*clhs95;
lhs(4,6)=DN(1,0)*clhs48 + DN(1,1)*clhs70 + clhs104;
lhs(4,7)=DN(1,0)*clhs56 + DN(1,1)*clhs72 + clhs105 + clhs98;
lhs(4,8)=k*(DN(1,1)*N[2]*bdf0*clhs7 - DN(2,1)*clhs89 - clhs106 - clhs107*clhs45);
lhs(5,0)=DN(1,0)*clhs76 + clhs43;
lhs(5,1)=DN(1,1)*clhs76 + clhs68;
lhs(5,2)=clhs82;
lhs(5,3)=DN(1,0)*clhs108;
lhs(5,4)=DN(1,1)*clhs108;
lhs(5,5)=bdf0*clhs91 - clhs102*clhs78 - clhs78*clhs90;
lhs(5,6)=DN(1,0)*clhs83 + clhs101;
lhs(5,7)=DN(1,1)*clhs83 + clhs107;
lhs(5,8)=clhs109;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs113 + clhs50;
lhs(6,1)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs71;
lhs(6,2)=k*(-DN(0,0)*clhs114 + DN(2,0)*N[0]*bdf0*clhs7 - clhs45*clhs58 - clhs59);
lhs(6,3)=DN(2,0)*clhs30 + DN(2,1)*clhs32 + clhs115 + clhs97;
lhs(6,4)=DN(2,0)*clhs39 + DN(2,1)*clhs41 + clhs104;
lhs(6,5)=k*(-DN(1,0)*clhs114 + DN(2,0)*N[1]*bdf0*clhs7 - clhs100*clhs45 - clhs101);
lhs(6,6)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + clhs116*clhs8 + clhs118;
lhs(6,7)=DN(2,0)*clhs54 + DN(2,1)*clhs56 + clhs119;
lhs(6,8)=DN(2,0)*clhs120;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs60 + clhs57;
lhs(7,1)=DN(2,0)*clhs23 + DN(2,1)*clhs61 + clhs113 + clhs73;
lhs(7,2)=k*(-DN(0,1)*clhs114 + DN(2,1)*N[0]*bdf0*clhs7 - clhs45*clhs74 - clhs75);
lhs(7,3)=DN(2,0)*clhs32 + DN(2,1)*clhs63 + clhs99;
lhs(7,4)=DN(2,0)*clhs41 + DN(2,1)*clhs66 + clhs105 + clhs115;
lhs(7,5)=k*(-DN(1,1)*clhs114 + DN(2,1)*N[1]*bdf0*clhs7 - clhs106*clhs45 - clhs107);
lhs(7,6)=DN(2,0)*clhs48 + DN(2,1)*clhs70 + clhs119;
lhs(7,7)=DN(2,0)*clhs56 + DN(2,1)*clhs72 + clhs118 + clhs121*clhs8;
lhs(7,8)=DN(2,1)*clhs120;
lhs(8,0)=DN(2,0)*clhs76 + clhs58;
lhs(8,1)=DN(2,1)*clhs76 + clhs74;
lhs(8,2)=clhs84;
lhs(8,3)=DN(2,0)*clhs79 + clhs100;
lhs(8,4)=DN(2,1)*clhs79 + clhs106;
lhs(8,5)=clhs109;
lhs(8,6)=DN(2,0)*clhs122;
lhs(8,7)=DN(2,1)*clhs122;
lhs(8,8)=bdf0*clhs117 - clhs116*clhs78 - clhs121*clhs78;


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
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 = rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 = clhs9*h/stab_c1 + mu;
const double clhs11 = clhs10*k;
const double clhs12 = DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs13 = N[0]*rho;
const double clhs14 = pow(N[0], 2);
const double clhs15 = N[0]*bdf0;
const double clhs16 = clhs12 + clhs15;
const double clhs17 = 1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = clhs17*pow(rho, 2);
const double clhs19 = clhs12*clhs18;
const double clhs20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs21 = clhs18*clhs20;
const double clhs22 = N[0]*clhs21;
const double clhs23 = bdf0*clhs14*rho + clhs12*clhs13 + clhs16*clhs19 + clhs16*clhs22;
const double clhs24 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs25 = C(1,3)*DN(0,1);
const double clhs26 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs25;
const double clhs27 = C(3,5)*DN(0,0);
const double clhs28 = C(4,5)*DN(0,2);
const double clhs29 = C(1,5)*DN(0,1) + clhs27 + clhs28;
const double clhs30 = DN(0,0)*clhs11;
const double clhs31 = DN(0,1)*clhs30;
const double clhs32 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs33 = C(3,4)*DN(0,1);
const double clhs34 = C(2,3)*DN(0,2) + clhs27 + clhs33;
const double clhs35 = C(2,5)*DN(0,2);
const double clhs36 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs35;
const double clhs37 = DN(0,2)*clhs30;
const double clhs38 = clhs17*clhs20;
const double clhs39 = clhs17*rho;
const double clhs40 = clhs12*clhs39;
const double clhs41 = k*(N[0]*bdf0*clhs10 - N[0] - clhs13*clhs38 - clhs40);
const double clhs42 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs43 = C(0,3)*DN(1,0);
const double clhs44 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs43;
const double clhs45 = C(0,5)*DN(1,0);
const double clhs46 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs45;
const double clhs47 = N[1]*clhs15*rho;
const double clhs48 = DN(1,0)*clhs30 + clhs47;
const double clhs49 = DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs50 = N[1]*bdf0;
const double clhs51 = clhs49 + clhs50;
const double clhs52 = clhs13*clhs49 + clhs19*clhs51 + clhs22*clhs51;
const double clhs53 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs43;
const double clhs54 = C(1,3)*DN(1,1);
const double clhs55 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs54;
const double clhs56 = C(3,5)*DN(1,0);
const double clhs57 = C(4,5)*DN(1,2);
const double clhs58 = C(1,5)*DN(1,1) + clhs56 + clhs57;
const double clhs59 = DN(1,1)*clhs30;
const double clhs60 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs45;
const double clhs61 = C(3,4)*DN(1,1);
const double clhs62 = C(2,3)*DN(1,2) + clhs56 + clhs61;
const double clhs63 = C(2,5)*DN(1,2);
const double clhs64 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs63;
const double clhs65 = DN(1,2)*clhs30;
const double clhs66 = DN(0,0)*N[1];
const double clhs67 = DN(1,0)*N[0];
const double clhs68 = clhs20*clhs39;
const double clhs69 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs70 = C(0,3)*DN(2,0);
const double clhs71 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs70;
const double clhs72 = C(0,5)*DN(2,0);
const double clhs73 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs72;
const double clhs74 = N[2]*clhs15*rho;
const double clhs75 = DN(2,0)*clhs30 + clhs74;
const double clhs76 = DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs77 = N[2]*bdf0;
const double clhs78 = clhs76 + clhs77;
const double clhs79 = clhs13*clhs76 + clhs19*clhs78 + clhs22*clhs78;
const double clhs80 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs70;
const double clhs81 = C(1,3)*DN(2,1);
const double clhs82 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs81;
const double clhs83 = C(3,5)*DN(2,0);
const double clhs84 = C(4,5)*DN(2,2);
const double clhs85 = C(1,5)*DN(2,1) + clhs83 + clhs84;
const double clhs86 = DN(2,1)*clhs30;
const double clhs87 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs72;
const double clhs88 = C(3,4)*DN(2,1);
const double clhs89 = C(2,3)*DN(2,2) + clhs83 + clhs88;
const double clhs90 = C(2,5)*DN(2,2);
const double clhs91 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs90;
const double clhs92 = DN(2,2)*clhs30;
const double clhs93 = DN(0,0)*N[2];
const double clhs94 = DN(2,0)*N[0];
const double clhs95 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs96 = C(0,3)*DN(3,0);
const double clhs97 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs96;
const double clhs98 = C(0,5)*DN(3,0);
const double clhs99 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs98;
const double clhs100 = N[3]*clhs15*rho;
const double clhs101 = DN(3,0)*clhs30 + clhs100;
const double clhs102 = DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs103 = N[3]*bdf0 + clhs102;
const double clhs104 = clhs102*clhs13 + clhs103*clhs19 + clhs103*clhs22;
const double clhs105 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs96;
const double clhs106 = C(1,3)*DN(3,1);
const double clhs107 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs106;
const double clhs108 = C(3,5)*DN(3,0);
const double clhs109 = C(4,5)*DN(3,2);
const double clhs110 = C(1,5)*DN(3,1) + clhs108 + clhs109;
const double clhs111 = DN(3,1)*clhs30;
const double clhs112 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs98;
const double clhs113 = C(3,4)*DN(3,1);
const double clhs114 = C(2,3)*DN(3,2) + clhs108 + clhs113;
const double clhs115 = C(2,5)*DN(3,2);
const double clhs116 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs115;
const double clhs117 = DN(3,2)*clhs30;
const double clhs118 = DN(0,0)*N[3];
const double clhs119 = DN(3,0)*N[0];
const double clhs120 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs25;
const double clhs121 = C(0,4)*DN(0,0) + clhs28 + clhs33;
const double clhs122 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs123 = C(1,4)*DN(0,1);
const double clhs124 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs123;
const double clhs125 = pow(DN(0,1), 2);
const double clhs126 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs123;
const double clhs127 = C(2,4)*DN(0,2);
const double clhs128 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs127;
const double clhs129 = DN(0,1)*clhs11;
const double clhs130 = DN(0,2)*clhs129;
const double clhs131 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs54;
const double clhs132 = C(0,4)*DN(1,0) + clhs57 + clhs61;
const double clhs133 = DN(1,0)*clhs129;
const double clhs134 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs135 = C(1,4)*DN(1,1);
const double clhs136 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs135;
const double clhs137 = DN(1,1)*clhs129;
const double clhs138 = clhs47 + clhs52;
const double clhs139 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs135;
const double clhs140 = C(2,4)*DN(1,2);
const double clhs141 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs140;
const double clhs142 = DN(1,2)*clhs129;
const double clhs143 = DN(0,1)*N[1];
const double clhs144 = DN(1,1)*N[0];
const double clhs145 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs81;
const double clhs146 = C(0,4)*DN(2,0) + clhs84 + clhs88;
const double clhs147 = DN(2,0)*clhs129;
const double clhs148 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs149 = C(1,4)*DN(2,1);
const double clhs150 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs149;
const double clhs151 = DN(2,1)*clhs129;
const double clhs152 = clhs74 + clhs79;
const double clhs153 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs149;
const double clhs154 = C(2,4)*DN(2,2);
const double clhs155 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs154;
const double clhs156 = DN(2,2)*clhs129;
const double clhs157 = DN(0,1)*N[2];
const double clhs158 = DN(2,1)*N[0];
const double clhs159 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs106;
const double clhs160 = C(0,4)*DN(3,0) + clhs109 + clhs113;
const double clhs161 = DN(3,0)*clhs129;
const double clhs162 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs163 = C(1,4)*DN(3,1);
const double clhs164 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs163;
const double clhs165 = DN(3,1)*clhs129;
const double clhs166 = clhs100 + clhs104;
const double clhs167 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs163;
const double clhs168 = C(2,4)*DN(3,2);
const double clhs169 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs168;
const double clhs170 = DN(3,2)*clhs129;
const double clhs171 = DN(0,1)*N[3];
const double clhs172 = DN(3,1)*N[0];
const double clhs173 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
const double clhs174 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs127;
const double clhs175 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs176 = pow(DN(0,2), 2);
const double clhs177 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs63;
const double clhs178 = DN(0,2)*clhs11;
const double clhs179 = DN(1,0)*clhs178;
const double clhs180 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs140;
const double clhs181 = DN(1,1)*clhs178;
const double clhs182 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs183 = DN(1,2)*clhs178;
const double clhs184 = DN(0,2)*N[1];
const double clhs185 = DN(1,2)*N[0];
const double clhs186 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs90;
const double clhs187 = DN(2,0)*clhs178;
const double clhs188 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs154;
const double clhs189 = DN(2,1)*clhs178;
const double clhs190 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs191 = DN(2,2)*clhs178;
const double clhs192 = DN(0,2)*N[2];
const double clhs193 = DN(2,2)*N[0];
const double clhs194 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs115;
const double clhs195 = DN(3,0)*clhs178;
const double clhs196 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs168;
const double clhs197 = DN(3,1)*clhs178;
const double clhs198 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs199 = DN(3,2)*clhs178;
const double clhs200 = DN(0,2)*N[3];
const double clhs201 = DN(3,2)*N[0];
const double clhs202 = clhs16*clhs39;
const double clhs203 = N[0] + clhs202;
const double clhs204 = clhs17*k;
const double clhs205 = clhs39*clhs51;
const double clhs206 = DN(0,0)*clhs204;
const double clhs207 = DN(0,1)*clhs204;
const double clhs208 = DN(0,2)*clhs204;
const double clhs209 = -DN(1,0)*clhs206 - DN(1,1)*clhs207 - DN(1,2)*clhs208 + N[0]*N[1]*bdf0;
const double clhs210 = clhs39*clhs78;
const double clhs211 = -DN(2,0)*clhs206 - DN(2,1)*clhs207 - DN(2,2)*clhs208 + N[0]*N[2]*bdf0;
const double clhs212 = clhs103*clhs39;
const double clhs213 = -DN(3,0)*clhs206 - DN(3,1)*clhs207 - DN(3,2)*clhs208 + N[0]*N[3]*bdf0;
const double clhs214 = N[1]*rho;
const double clhs215 = clhs18*clhs49;
const double clhs216 = N[1]*clhs21;
const double clhs217 = clhs12*clhs214 + clhs16*clhs215 + clhs16*clhs216;
const double clhs218 = clhs39*clhs49;
const double clhs219 = pow(DN(1,0), 2);
const double clhs220 = pow(N[1], 2);
const double clhs221 = bdf0*clhs220*rho + clhs214*clhs49 + clhs215*clhs51 + clhs216*clhs51;
const double clhs222 = DN(1,0)*clhs11;
const double clhs223 = DN(1,1)*clhs222;
const double clhs224 = DN(1,2)*clhs222;
const double clhs225 = k*(N[1]*bdf0*clhs10 - N[1] - clhs214*clhs38 - clhs218);
const double clhs226 = N[2]*clhs50*rho;
const double clhs227 = DN(2,0)*clhs222 + clhs226;
const double clhs228 = clhs214*clhs76 + clhs215*clhs78 + clhs216*clhs78;
const double clhs229 = DN(2,1)*clhs222;
const double clhs230 = DN(2,2)*clhs222;
const double clhs231 = DN(1,0)*N[2];
const double clhs232 = DN(2,0)*N[1];
const double clhs233 = N[3]*clhs50*rho;
const double clhs234 = DN(3,0)*clhs222 + clhs233;
const double clhs235 = clhs102*clhs214 + clhs103*clhs215 + clhs103*clhs216;
const double clhs236 = DN(3,1)*clhs222;
const double clhs237 = DN(3,2)*clhs222;
const double clhs238 = DN(1,0)*N[3];
const double clhs239 = DN(3,0)*N[1];
const double clhs240 = clhs217 + clhs47;
const double clhs241 = pow(DN(1,1), 2);
const double clhs242 = DN(1,1)*clhs11;
const double clhs243 = DN(1,2)*clhs242;
const double clhs244 = DN(2,0)*clhs242;
const double clhs245 = DN(2,1)*clhs242;
const double clhs246 = clhs226 + clhs228;
const double clhs247 = DN(2,2)*clhs242;
const double clhs248 = DN(1,1)*N[2];
const double clhs249 = DN(2,1)*N[1];
const double clhs250 = DN(3,0)*clhs242;
const double clhs251 = DN(3,1)*clhs242;
const double clhs252 = clhs233 + clhs235;
const double clhs253 = DN(3,2)*clhs242;
const double clhs254 = DN(1,1)*N[3];
const double clhs255 = DN(3,1)*N[1];
const double clhs256 = pow(DN(1,2), 2);
const double clhs257 = DN(1,2)*clhs11;
const double clhs258 = DN(2,0)*clhs257;
const double clhs259 = DN(2,1)*clhs257;
const double clhs260 = DN(2,2)*clhs257;
const double clhs261 = DN(1,2)*N[2];
const double clhs262 = DN(2,2)*N[1];
const double clhs263 = DN(3,0)*clhs257;
const double clhs264 = DN(3,1)*clhs257;
const double clhs265 = DN(3,2)*clhs257;
const double clhs266 = DN(1,2)*N[3];
const double clhs267 = DN(3,2)*N[1];
const double clhs268 = N[1] + clhs205;
const double clhs269 = DN(1,0)*clhs204;
const double clhs270 = DN(1,1)*clhs204;
const double clhs271 = DN(1,2)*clhs204;
const double clhs272 = -DN(2,0)*clhs269 - DN(2,1)*clhs270 - DN(2,2)*clhs271 + N[1]*N[2]*bdf0;
const double clhs273 = -DN(3,0)*clhs269 - DN(3,1)*clhs270 - DN(3,2)*clhs271 + N[1]*N[3]*bdf0;
const double clhs274 = N[2]*rho;
const double clhs275 = clhs18*clhs76;
const double clhs276 = N[2]*clhs21;
const double clhs277 = clhs12*clhs274 + clhs16*clhs275 + clhs16*clhs276;
const double clhs278 = clhs39*clhs76;
const double clhs279 = clhs274*clhs49 + clhs275*clhs51 + clhs276*clhs51;
const double clhs280 = pow(DN(2,0), 2);
const double clhs281 = pow(N[2], 2);
const double clhs282 = bdf0*clhs281*rho + clhs274*clhs76 + clhs275*clhs78 + clhs276*clhs78;
const double clhs283 = DN(2,0)*clhs11;
const double clhs284 = DN(2,1)*clhs283;
const double clhs285 = DN(2,2)*clhs283;
const double clhs286 = k*(N[2]*bdf0*clhs10 - N[2] - clhs274*clhs38 - clhs278);
const double clhs287 = N[3]*clhs77*rho;
const double clhs288 = DN(3,0)*clhs283 + clhs287;
const double clhs289 = clhs102*clhs274 + clhs103*clhs275 + clhs103*clhs276;
const double clhs290 = DN(3,1)*clhs283;
const double clhs291 = DN(3,2)*clhs283;
const double clhs292 = DN(2,0)*N[3];
const double clhs293 = DN(3,0)*N[2];
const double clhs294 = clhs277 + clhs74;
const double clhs295 = clhs226 + clhs279;
const double clhs296 = pow(DN(2,1), 2);
const double clhs297 = DN(2,1)*clhs11;
const double clhs298 = DN(2,2)*clhs297;
const double clhs299 = DN(3,0)*clhs297;
const double clhs300 = DN(3,1)*clhs297;
const double clhs301 = clhs287 + clhs289;
const double clhs302 = DN(3,2)*clhs297;
const double clhs303 = DN(2,1)*N[3];
const double clhs304 = DN(3,1)*N[2];
const double clhs305 = pow(DN(2,2), 2);
const double clhs306 = DN(2,2)*clhs11;
const double clhs307 = DN(3,0)*clhs306;
const double clhs308 = DN(3,1)*clhs306;
const double clhs309 = DN(3,2)*clhs306;
const double clhs310 = DN(2,2)*N[3];
const double clhs311 = DN(3,2)*N[2];
const double clhs312 = N[2] + clhs210;
const double clhs313 = -DN(2,0)*DN(3,0)*clhs204 - DN(2,1)*DN(3,1)*clhs204 - DN(2,2)*DN(3,2)*clhs204 + N[2]*N[3]*bdf0;
const double clhs314 = N[3]*rho;
const double clhs315 = clhs102*clhs18;
const double clhs316 = N[3]*clhs21;
const double clhs317 = clhs12*clhs314 + clhs16*clhs315 + clhs16*clhs316;
const double clhs318 = clhs102*clhs39;
const double clhs319 = clhs314*clhs49 + clhs315*clhs51 + clhs316*clhs51;
const double clhs320 = clhs314*clhs76 + clhs315*clhs78 + clhs316*clhs78;
const double clhs321 = pow(DN(3,0), 2);
const double clhs322 = pow(N[3], 2);
const double clhs323 = bdf0*clhs322*rho + clhs102*clhs314 + clhs103*clhs315 + clhs103*clhs316;
const double clhs324 = DN(3,0)*clhs11;
const double clhs325 = DN(3,1)*clhs324;
const double clhs326 = DN(3,2)*clhs324;
const double clhs327 = k*(N[3]*bdf0*clhs10 - N[3] - clhs314*clhs38 - clhs318);
const double clhs328 = clhs100 + clhs317;
const double clhs329 = clhs233 + clhs319;
const double clhs330 = clhs287 + clhs320;
const double clhs331 = pow(DN(3,1), 2);
const double clhs332 = DN(3,1)*DN(3,2)*clhs11;
const double clhs333 = pow(DN(3,2), 2);
const double clhs334 = N[3] + clhs212;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs23;
lhs(0,1)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + DN(0,2)*clhs29 + clhs31;
lhs(0,2)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs37;
lhs(0,3)=DN(0,0)*clhs41;
lhs(0,4)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + DN(0,2)*clhs46 + clhs48 + clhs52;
lhs(0,5)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + DN(0,2)*clhs58 + clhs59;
lhs(0,6)=DN(0,0)*clhs60 + DN(0,1)*clhs62 + DN(0,2)*clhs64 + clhs65;
lhs(0,7)=k*(DN(0,0)*N[1]*bdf0*clhs10 - DN(1,0)*clhs40 - clhs66 - clhs67*clhs68);
lhs(0,8)=DN(0,0)*clhs69 + DN(0,1)*clhs71 + DN(0,2)*clhs73 + clhs75 + clhs79;
lhs(0,9)=DN(0,0)*clhs80 + DN(0,1)*clhs82 + DN(0,2)*clhs85 + clhs86;
lhs(0,10)=DN(0,0)*clhs87 + DN(0,1)*clhs89 + DN(0,2)*clhs91 + clhs92;
lhs(0,11)=k*(DN(0,0)*N[2]*bdf0*clhs10 - DN(2,0)*clhs40 - clhs68*clhs94 - clhs93);
lhs(0,12)=DN(0,0)*clhs95 + DN(0,1)*clhs97 + DN(0,2)*clhs99 + clhs101 + clhs104;
lhs(0,13)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs110 + clhs111;
lhs(0,14)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs117;
lhs(0,15)=k*(DN(0,0)*N[3]*bdf0*clhs10 - DN(3,0)*clhs40 - clhs118 - clhs119*clhs68);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs120 + DN(0,2)*clhs121 + clhs31;
lhs(1,1)=DN(0,0)*clhs26 + DN(0,1)*clhs122 + DN(0,2)*clhs124 + clhs11*clhs125 + clhs23;
lhs(1,2)=DN(0,0)*clhs34 + DN(0,1)*clhs126 + DN(0,2)*clhs128 + clhs130;
lhs(1,3)=DN(0,1)*clhs41;
lhs(1,4)=DN(0,0)*clhs44 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs133;
lhs(1,5)=DN(0,0)*clhs55 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs137 + clhs138;
lhs(1,6)=DN(0,0)*clhs62 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142;
lhs(1,7)=k*(DN(0,1)*N[1]*bdf0*clhs10 - DN(1,1)*clhs40 - clhs143 - clhs144*clhs68);
lhs(1,8)=DN(0,0)*clhs71 + DN(0,1)*clhs145 + DN(0,2)*clhs146 + clhs147;
lhs(1,9)=DN(0,0)*clhs82 + DN(0,1)*clhs148 + DN(0,2)*clhs150 + clhs151 + clhs152;
lhs(1,10)=DN(0,0)*clhs89 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs156;
lhs(1,11)=k*(DN(0,1)*N[2]*bdf0*clhs10 - DN(2,1)*clhs40 - clhs157 - clhs158*clhs68);
lhs(1,12)=DN(0,0)*clhs97 + DN(0,1)*clhs159 + DN(0,2)*clhs160 + clhs161;
lhs(1,13)=DN(0,0)*clhs107 + DN(0,1)*clhs162 + DN(0,2)*clhs164 + clhs165 + clhs166;
lhs(1,14)=DN(0,0)*clhs114 + DN(0,1)*clhs167 + DN(0,2)*clhs169 + clhs170;
lhs(1,15)=k*(DN(0,1)*N[3]*bdf0*clhs10 - DN(3,1)*clhs40 - clhs171 - clhs172*clhs68);
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs121 + DN(0,2)*clhs173 + clhs37;
lhs(2,1)=DN(0,0)*clhs29 + DN(0,1)*clhs124 + DN(0,2)*clhs174 + clhs130;
lhs(2,2)=DN(0,0)*clhs36 + DN(0,1)*clhs128 + DN(0,2)*clhs175 + clhs11*clhs176 + clhs23;
lhs(2,3)=DN(0,2)*clhs41;
lhs(2,4)=DN(0,0)*clhs46 + DN(0,1)*clhs132 + DN(0,2)*clhs177 + clhs179;
lhs(2,5)=DN(0,0)*clhs58 + DN(0,1)*clhs136 + DN(0,2)*clhs180 + clhs181;
lhs(2,6)=DN(0,0)*clhs64 + DN(0,1)*clhs141 + DN(0,2)*clhs182 + clhs138 + clhs183;
lhs(2,7)=k*(DN(0,2)*N[1]*bdf0*clhs10 - DN(1,2)*clhs40 - clhs184 - clhs185*clhs68);
lhs(2,8)=DN(0,0)*clhs73 + DN(0,1)*clhs146 + DN(0,2)*clhs186 + clhs187;
lhs(2,9)=DN(0,0)*clhs85 + DN(0,1)*clhs150 + DN(0,2)*clhs188 + clhs189;
lhs(2,10)=DN(0,0)*clhs91 + DN(0,1)*clhs155 + DN(0,2)*clhs190 + clhs152 + clhs191;
lhs(2,11)=k*(DN(0,2)*N[2]*bdf0*clhs10 - DN(2,2)*clhs40 - clhs192 - clhs193*clhs68);
lhs(2,12)=DN(0,0)*clhs99 + DN(0,1)*clhs160 + DN(0,2)*clhs194 + clhs195;
lhs(2,13)=DN(0,0)*clhs110 + DN(0,1)*clhs164 + DN(0,2)*clhs196 + clhs197;
lhs(2,14)=DN(0,0)*clhs116 + DN(0,1)*clhs169 + DN(0,2)*clhs198 + clhs166 + clhs199;
lhs(2,15)=k*(DN(0,2)*N[3]*bdf0*clhs10 - DN(3,2)*clhs40 - clhs200 - clhs201*clhs68);
lhs(3,0)=DN(0,0)*clhs203;
lhs(3,1)=DN(0,1)*clhs203;
lhs(3,2)=DN(0,2)*clhs203;
lhs(3,3)=bdf0*clhs14 - clhs125*clhs204 - clhs176*clhs204 - clhs204*clhs5;
lhs(3,4)=DN(0,0)*clhs205 + clhs67;
lhs(3,5)=DN(0,1)*clhs205 + clhs144;
lhs(3,6)=DN(0,2)*clhs205 + clhs185;
lhs(3,7)=clhs209;
lhs(3,8)=DN(0,0)*clhs210 + clhs94;
lhs(3,9)=DN(0,1)*clhs210 + clhs158;
lhs(3,10)=DN(0,2)*clhs210 + clhs193;
lhs(3,11)=clhs211;
lhs(3,12)=DN(0,0)*clhs212 + clhs119;
lhs(3,13)=DN(0,1)*clhs212 + clhs172;
lhs(3,14)=DN(0,2)*clhs212 + clhs201;
lhs(3,15)=clhs213;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs217 + clhs48;
lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + DN(1,2)*clhs29 + clhs133;
lhs(4,2)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs179;
lhs(4,3)=k*(-DN(0,0)*clhs218 + DN(1,0)*N[0]*bdf0*clhs10 - clhs66*clhs68 - clhs67);
lhs(4,4)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + DN(1,2)*clhs46 + clhs11*clhs219 + clhs221;
lhs(4,5)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + DN(1,2)*clhs58 + clhs223;
lhs(4,6)=DN(1,0)*clhs60 + DN(1,1)*clhs62 + DN(1,2)*clhs64 + clhs224;
lhs(4,7)=DN(1,0)*clhs225;
lhs(4,8)=DN(1,0)*clhs69 + DN(1,1)*clhs71 + DN(1,2)*clhs73 + clhs227 + clhs228;
lhs(4,9)=DN(1,0)*clhs80 + DN(1,1)*clhs82 + DN(1,2)*clhs85 + clhs229;
lhs(4,10)=DN(1,0)*clhs87 + DN(1,1)*clhs89 + DN(1,2)*clhs91 + clhs230;
lhs(4,11)=k*(DN(1,0)*N[2]*bdf0*clhs10 - DN(2,0)*clhs218 - clhs231 - clhs232*clhs68);
lhs(4,12)=DN(1,0)*clhs95 + DN(1,1)*clhs97 + DN(1,2)*clhs99 + clhs234 + clhs235;
lhs(4,13)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs110 + clhs236;
lhs(4,14)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs237;
lhs(4,15)=k*(DN(1,0)*N[3]*bdf0*clhs10 - DN(3,0)*clhs218 - clhs238 - clhs239*clhs68);
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs120 + DN(1,2)*clhs121 + clhs59;
lhs(5,1)=DN(1,0)*clhs26 + DN(1,1)*clhs122 + DN(1,2)*clhs124 + clhs137 + clhs240;
lhs(5,2)=DN(1,0)*clhs34 + DN(1,1)*clhs126 + DN(1,2)*clhs128 + clhs181;
lhs(5,3)=k*(-DN(0,1)*clhs218 + DN(1,1)*N[0]*bdf0*clhs10 - clhs143*clhs68 - clhs144);
lhs(5,4)=DN(1,0)*clhs44 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs223;
lhs(5,5)=DN(1,0)*clhs55 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs11*clhs241 + clhs221;
lhs(5,6)=DN(1,0)*clhs62 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs243;
lhs(5,7)=DN(1,1)*clhs225;
lhs(5,8)=DN(1,0)*clhs71 + DN(1,1)*clhs145 + DN(1,2)*clhs146 + clhs244;
lhs(5,9)=DN(1,0)*clhs82 + DN(1,1)*clhs148 + DN(1,2)*clhs150 + clhs245 + clhs246;
lhs(5,10)=DN(1,0)*clhs89 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs247;
lhs(5,11)=k*(DN(1,1)*N[2]*bdf0*clhs10 - DN(2,1)*clhs218 - clhs248 - clhs249*clhs68);
lhs(5,12)=DN(1,0)*clhs97 + DN(1,1)*clhs159 + DN(1,2)*clhs160 + clhs250;
lhs(5,13)=DN(1,0)*clhs107 + DN(1,1)*clhs162 + DN(1,2)*clhs164 + clhs251 + clhs252;
lhs(5,14)=DN(1,0)*clhs114 + DN(1,1)*clhs167 + DN(1,2)*clhs169 + clhs253;
lhs(5,15)=k*(DN(1,1)*N[3]*bdf0*clhs10 - DN(3,1)*clhs218 - clhs254 - clhs255*clhs68);
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs121 + DN(1,2)*clhs173 + clhs65;
lhs(6,1)=DN(1,0)*clhs29 + DN(1,1)*clhs124 + DN(1,2)*clhs174 + clhs142;
lhs(6,2)=DN(1,0)*clhs36 + DN(1,1)*clhs128 + DN(1,2)*clhs175 + clhs183 + clhs240;
lhs(6,3)=k*(-DN(0,2)*clhs218 + DN(1,2)*N[0]*bdf0*clhs10 - clhs184*clhs68 - clhs185);
lhs(6,4)=DN(1,0)*clhs46 + DN(1,1)*clhs132 + DN(1,2)*clhs177 + clhs224;
lhs(6,5)=DN(1,0)*clhs58 + DN(1,1)*clhs136 + DN(1,2)*clhs180 + clhs243;
lhs(6,6)=DN(1,0)*clhs64 + DN(1,1)*clhs141 + DN(1,2)*clhs182 + clhs11*clhs256 + clhs221;
lhs(6,7)=DN(1,2)*clhs225;
lhs(6,8)=DN(1,0)*clhs73 + DN(1,1)*clhs146 + DN(1,2)*clhs186 + clhs258;
lhs(6,9)=DN(1,0)*clhs85 + DN(1,1)*clhs150 + DN(1,2)*clhs188 + clhs259;
lhs(6,10)=DN(1,0)*clhs91 + DN(1,1)*clhs155 + DN(1,2)*clhs190 + clhs246 + clhs260;
lhs(6,11)=k*(DN(1,2)*N[2]*bdf0*clhs10 - DN(2,2)*clhs218 - clhs261 - clhs262*clhs68);
lhs(6,12)=DN(1,0)*clhs99 + DN(1,1)*clhs160 + DN(1,2)*clhs194 + clhs263;
lhs(6,13)=DN(1,0)*clhs110 + DN(1,1)*clhs164 + DN(1,2)*clhs196 + clhs264;
lhs(6,14)=DN(1,0)*clhs116 + DN(1,1)*clhs169 + DN(1,2)*clhs198 + clhs252 + clhs265;
lhs(6,15)=k*(DN(1,2)*N[3]*bdf0*clhs10 - DN(3,2)*clhs218 - clhs266 - clhs267*clhs68);
lhs(7,0)=DN(1,0)*clhs202 + clhs66;
lhs(7,1)=DN(1,1)*clhs202 + clhs143;
lhs(7,2)=DN(1,2)*clhs202 + clhs184;
lhs(7,3)=clhs209;
lhs(7,4)=DN(1,0)*clhs268;
lhs(7,5)=DN(1,1)*clhs268;
lhs(7,6)=DN(1,2)*clhs268;
lhs(7,7)=bdf0*clhs220 - clhs204*clhs219 - clhs204*clhs241 - clhs204*clhs256;
lhs(7,8)=DN(1,0)*clhs210 + clhs232;
lhs(7,9)=DN(1,1)*clhs210 + clhs249;
lhs(7,10)=DN(1,2)*clhs210 + clhs262;
lhs(7,11)=clhs272;
lhs(7,12)=DN(1,0)*clhs212 + clhs239;
lhs(7,13)=DN(1,1)*clhs212 + clhs255;
lhs(7,14)=DN(1,2)*clhs212 + clhs267;
lhs(7,15)=clhs273;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs277 + clhs75;
lhs(8,1)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + DN(2,2)*clhs29 + clhs147;
lhs(8,2)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs187;
lhs(8,3)=k*(-DN(0,0)*clhs278 + DN(2,0)*N[0]*bdf0*clhs10 - clhs68*clhs93 - clhs94);
lhs(8,4)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + DN(2,2)*clhs46 + clhs227 + clhs279;
lhs(8,5)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + DN(2,2)*clhs58 + clhs244;
lhs(8,6)=DN(2,0)*clhs60 + DN(2,1)*clhs62 + DN(2,2)*clhs64 + clhs258;
lhs(8,7)=k*(-DN(1,0)*clhs278 + DN(2,0)*N[1]*bdf0*clhs10 - clhs231*clhs68 - clhs232);
lhs(8,8)=DN(2,0)*clhs69 + DN(2,1)*clhs71 + DN(2,2)*clhs73 + clhs11*clhs280 + clhs282;
lhs(8,9)=DN(2,0)*clhs80 + DN(2,1)*clhs82 + DN(2,2)*clhs85 + clhs284;
lhs(8,10)=DN(2,0)*clhs87 + DN(2,1)*clhs89 + DN(2,2)*clhs91 + clhs285;
lhs(8,11)=DN(2,0)*clhs286;
lhs(8,12)=DN(2,0)*clhs95 + DN(2,1)*clhs97 + DN(2,2)*clhs99 + clhs288 + clhs289;
lhs(8,13)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs110 + clhs290;
lhs(8,14)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs291;
lhs(8,15)=k*(DN(2,0)*N[3]*bdf0*clhs10 - DN(3,0)*clhs278 - clhs292 - clhs293*clhs68);
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs120 + DN(2,2)*clhs121 + clhs86;
lhs(9,1)=DN(2,0)*clhs26 + DN(2,1)*clhs122 + DN(2,2)*clhs124 + clhs151 + clhs294;
lhs(9,2)=DN(2,0)*clhs34 + DN(2,1)*clhs126 + DN(2,2)*clhs128 + clhs189;
lhs(9,3)=k*(-DN(0,1)*clhs278 + DN(2,1)*N[0]*bdf0*clhs10 - clhs157*clhs68 - clhs158);
lhs(9,4)=DN(2,0)*clhs44 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs229;
lhs(9,5)=DN(2,0)*clhs55 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs245 + clhs295;
lhs(9,6)=DN(2,0)*clhs62 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs259;
lhs(9,7)=k*(-DN(1,1)*clhs278 + DN(2,1)*N[1]*bdf0*clhs10 - clhs248*clhs68 - clhs249);
lhs(9,8)=DN(2,0)*clhs71 + DN(2,1)*clhs145 + DN(2,2)*clhs146 + clhs284;
lhs(9,9)=DN(2,0)*clhs82 + DN(2,1)*clhs148 + DN(2,2)*clhs150 + clhs11*clhs296 + clhs282;
lhs(9,10)=DN(2,0)*clhs89 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs298;
lhs(9,11)=DN(2,1)*clhs286;
lhs(9,12)=DN(2,0)*clhs97 + DN(2,1)*clhs159 + DN(2,2)*clhs160 + clhs299;
lhs(9,13)=DN(2,0)*clhs107 + DN(2,1)*clhs162 + DN(2,2)*clhs164 + clhs300 + clhs301;
lhs(9,14)=DN(2,0)*clhs114 + DN(2,1)*clhs167 + DN(2,2)*clhs169 + clhs302;
lhs(9,15)=k*(DN(2,1)*N[3]*bdf0*clhs10 - DN(3,1)*clhs278 - clhs303 - clhs304*clhs68);
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs121 + DN(2,2)*clhs173 + clhs92;
lhs(10,1)=DN(2,0)*clhs29 + DN(2,1)*clhs124 + DN(2,2)*clhs174 + clhs156;
lhs(10,2)=DN(2,0)*clhs36 + DN(2,1)*clhs128 + DN(2,2)*clhs175 + clhs191 + clhs294;
lhs(10,3)=k*(-DN(0,2)*clhs278 + DN(2,2)*N[0]*bdf0*clhs10 - clhs192*clhs68 - clhs193);
lhs(10,4)=DN(2,0)*clhs46 + DN(2,1)*clhs132 + DN(2,2)*clhs177 + clhs230;
lhs(10,5)=DN(2,0)*clhs58 + DN(2,1)*clhs136 + DN(2,2)*clhs180 + clhs247;
lhs(10,6)=DN(2,0)*clhs64 + DN(2,1)*clhs141 + DN(2,2)*clhs182 + clhs260 + clhs295;
lhs(10,7)=k*(-DN(1,2)*clhs278 + DN(2,2)*N[1]*bdf0*clhs10 - clhs261*clhs68 - clhs262);
lhs(10,8)=DN(2,0)*clhs73 + DN(2,1)*clhs146 + DN(2,2)*clhs186 + clhs285;
lhs(10,9)=DN(2,0)*clhs85 + DN(2,1)*clhs150 + DN(2,2)*clhs188 + clhs298;
lhs(10,10)=DN(2,0)*clhs91 + DN(2,1)*clhs155 + DN(2,2)*clhs190 + clhs11*clhs305 + clhs282;
lhs(10,11)=DN(2,2)*clhs286;
lhs(10,12)=DN(2,0)*clhs99 + DN(2,1)*clhs160 + DN(2,2)*clhs194 + clhs307;
lhs(10,13)=DN(2,0)*clhs110 + DN(2,1)*clhs164 + DN(2,2)*clhs196 + clhs308;
lhs(10,14)=DN(2,0)*clhs116 + DN(2,1)*clhs169 + DN(2,2)*clhs198 + clhs301 + clhs309;
lhs(10,15)=k*(DN(2,2)*N[3]*bdf0*clhs10 - DN(3,2)*clhs278 - clhs310 - clhs311*clhs68);
lhs(11,0)=DN(2,0)*clhs202 + clhs93;
lhs(11,1)=DN(2,1)*clhs202 + clhs157;
lhs(11,2)=DN(2,2)*clhs202 + clhs192;
lhs(11,3)=clhs211;
lhs(11,4)=DN(2,0)*clhs205 + clhs231;
lhs(11,5)=DN(2,1)*clhs205 + clhs248;
lhs(11,6)=DN(2,2)*clhs205 + clhs261;
lhs(11,7)=clhs272;
lhs(11,8)=DN(2,0)*clhs312;
lhs(11,9)=DN(2,1)*clhs312;
lhs(11,10)=DN(2,2)*clhs312;
lhs(11,11)=bdf0*clhs281 - clhs204*clhs280 - clhs204*clhs296 - clhs204*clhs305;
lhs(11,12)=DN(2,0)*clhs212 + clhs293;
lhs(11,13)=DN(2,1)*clhs212 + clhs304;
lhs(11,14)=DN(2,2)*clhs212 + clhs311;
lhs(11,15)=clhs313;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs101 + clhs317;
lhs(12,1)=DN(3,0)*clhs24 + DN(3,1)*clhs26 + DN(3,2)*clhs29 + clhs161;
lhs(12,2)=DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs195;
lhs(12,3)=k*(-DN(0,0)*clhs318 + DN(3,0)*N[0]*bdf0*clhs10 - clhs118*clhs68 - clhs119);
lhs(12,4)=DN(3,0)*clhs42 + DN(3,1)*clhs44 + DN(3,2)*clhs46 + clhs234 + clhs319;
lhs(12,5)=DN(3,0)*clhs53 + DN(3,1)*clhs55 + DN(3,2)*clhs58 + clhs250;
lhs(12,6)=DN(3,0)*clhs60 + DN(3,1)*clhs62 + DN(3,2)*clhs64 + clhs263;
lhs(12,7)=k*(-DN(1,0)*clhs318 + DN(3,0)*N[1]*bdf0*clhs10 - clhs238*clhs68 - clhs239);
lhs(12,8)=DN(3,0)*clhs69 + DN(3,1)*clhs71 + DN(3,2)*clhs73 + clhs288 + clhs320;
lhs(12,9)=DN(3,0)*clhs80 + DN(3,1)*clhs82 + DN(3,2)*clhs85 + clhs299;
lhs(12,10)=DN(3,0)*clhs87 + DN(3,1)*clhs89 + DN(3,2)*clhs91 + clhs307;
lhs(12,11)=k*(-DN(2,0)*clhs318 + DN(3,0)*N[2]*bdf0*clhs10 - clhs292*clhs68 - clhs293);
lhs(12,12)=DN(3,0)*clhs95 + DN(3,1)*clhs97 + DN(3,2)*clhs99 + clhs11*clhs321 + clhs323;
lhs(12,13)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs110 + clhs325;
lhs(12,14)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs326;
lhs(12,15)=DN(3,0)*clhs327;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs120 + DN(3,2)*clhs121 + clhs111;
lhs(13,1)=DN(3,0)*clhs26 + DN(3,1)*clhs122 + DN(3,2)*clhs124 + clhs165 + clhs328;
lhs(13,2)=DN(3,0)*clhs34 + DN(3,1)*clhs126 + DN(3,2)*clhs128 + clhs197;
lhs(13,3)=k*(-DN(0,1)*clhs318 + DN(3,1)*N[0]*bdf0*clhs10 - clhs171*clhs68 - clhs172);
lhs(13,4)=DN(3,0)*clhs44 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs236;
lhs(13,5)=DN(3,0)*clhs55 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs251 + clhs329;
lhs(13,6)=DN(3,0)*clhs62 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs264;
lhs(13,7)=k*(-DN(1,1)*clhs318 + DN(3,1)*N[1]*bdf0*clhs10 - clhs254*clhs68 - clhs255);
lhs(13,8)=DN(3,0)*clhs71 + DN(3,1)*clhs145 + DN(3,2)*clhs146 + clhs290;
lhs(13,9)=DN(3,0)*clhs82 + DN(3,1)*clhs148 + DN(3,2)*clhs150 + clhs300 + clhs330;
lhs(13,10)=DN(3,0)*clhs89 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs308;
lhs(13,11)=k*(-DN(2,1)*clhs318 + DN(3,1)*N[2]*bdf0*clhs10 - clhs303*clhs68 - clhs304);
lhs(13,12)=DN(3,0)*clhs97 + DN(3,1)*clhs159 + DN(3,2)*clhs160 + clhs325;
lhs(13,13)=DN(3,0)*clhs107 + DN(3,1)*clhs162 + DN(3,2)*clhs164 + clhs11*clhs331 + clhs323;
lhs(13,14)=DN(3,0)*clhs114 + DN(3,1)*clhs167 + DN(3,2)*clhs169 + clhs332;
lhs(13,15)=DN(3,1)*clhs327;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs121 + DN(3,2)*clhs173 + clhs117;
lhs(14,1)=DN(3,0)*clhs29 + DN(3,1)*clhs124 + DN(3,2)*clhs174 + clhs170;
lhs(14,2)=DN(3,0)*clhs36 + DN(3,1)*clhs128 + DN(3,2)*clhs175 + clhs199 + clhs328;
lhs(14,3)=k*(-DN(0,2)*clhs318 + DN(3,2)*N[0]*bdf0*clhs10 - clhs200*clhs68 - clhs201);
lhs(14,4)=DN(3,0)*clhs46 + DN(3,1)*clhs132 + DN(3,2)*clhs177 + clhs237;
lhs(14,5)=DN(3,0)*clhs58 + DN(3,1)*clhs136 + DN(3,2)*clhs180 + clhs253;
lhs(14,6)=DN(3,0)*clhs64 + DN(3,1)*clhs141 + DN(3,2)*clhs182 + clhs265 + clhs329;
lhs(14,7)=k*(-DN(1,2)*clhs318 + DN(3,2)*N[1]*bdf0*clhs10 - clhs266*clhs68 - clhs267);
lhs(14,8)=DN(3,0)*clhs73 + DN(3,1)*clhs146 + DN(3,2)*clhs186 + clhs291;
lhs(14,9)=DN(3,0)*clhs85 + DN(3,1)*clhs150 + DN(3,2)*clhs188 + clhs302;
lhs(14,10)=DN(3,0)*clhs91 + DN(3,1)*clhs155 + DN(3,2)*clhs190 + clhs309 + clhs330;
lhs(14,11)=k*(-DN(2,2)*clhs318 + DN(3,2)*N[2]*bdf0*clhs10 - clhs310*clhs68 - clhs311);
lhs(14,12)=DN(3,0)*clhs99 + DN(3,1)*clhs160 + DN(3,2)*clhs194 + clhs326;
lhs(14,13)=DN(3,0)*clhs110 + DN(3,1)*clhs164 + DN(3,2)*clhs196 + clhs332;
lhs(14,14)=DN(3,0)*clhs116 + DN(3,1)*clhs169 + DN(3,2)*clhs198 + clhs11*clhs333 + clhs323;
lhs(14,15)=DN(3,2)*clhs327;
lhs(15,0)=DN(3,0)*clhs202 + clhs118;
lhs(15,1)=DN(3,1)*clhs202 + clhs171;
lhs(15,2)=DN(3,2)*clhs202 + clhs200;
lhs(15,3)=clhs213;
lhs(15,4)=DN(3,0)*clhs205 + clhs238;
lhs(15,5)=DN(3,1)*clhs205 + clhs254;
lhs(15,6)=DN(3,2)*clhs205 + clhs266;
lhs(15,7)=clhs273;
lhs(15,8)=DN(3,0)*clhs210 + clhs292;
lhs(15,9)=DN(3,1)*clhs210 + clhs303;
lhs(15,10)=DN(3,2)*clhs210 + clhs310;
lhs(15,11)=clhs313;
lhs(15,12)=DN(3,0)*clhs334;
lhs(15,13)=DN(3,1)*clhs334;
lhs(15,14)=DN(3,2)*clhs334;
lhs(15,15)=bdf0*clhs322 - clhs204*clhs321 - clhs204*clhs331 - clhs204*clhs333;


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

    const double crhs0 = N[0]*eps_vol[0] + N[1]*eps_vol[1] + N[2]*eps_vol[2];
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0);
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs8 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]) + N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]) + N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]);
const double crhs9 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs10 = crhs3 + crhs9;
const double crhs11 = k*(crhs10 + crhs8)*(crhs7*h/stab_c1 + mu);
const double crhs12 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs13 = 1.0/(crhs7/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs14 = crhs1*rho - crhs2 - crhs6 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2]);
const double crhs15 = DN(0,0)*crhs4 + DN(0,1)*crhs5;
const double crhs16 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs17 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs18 = rho*(crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs5*crhs9);
const double crhs19 = crhs16*rho - crhs17 - crhs18 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2]);
const double crhs20 = 1.0*crhs13;
const double crhs21 = crhs14*crhs20;
const double crhs22 = crhs19*crhs20;
const double crhs23 = DN(1,0)*crhs4 + DN(1,1)*crhs5;
const double crhs24 = DN(2,0)*crhs4 + DN(2,1)*crhs5;
rhs[0]=DN(0,0)*crhs0*k - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1*rho + 1.0*N[0]*crhs12*crhs13*crhs14*rho - N[0]*crhs2 - N[0]*crhs6 + 1.0*crhs13*crhs14*crhs15*rho;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0*k - DN(0,1)*crhs11 - DN(0,1)*stress[1] + 1.0*N[0]*crhs12*crhs13*crhs19*rho + N[0]*crhs16*rho - N[0]*crhs17 - N[0]*crhs18 + 1.0*crhs13*crhs15*crhs19*rho;
rhs[2]=DN(0,0)*crhs21 + DN(0,1)*crhs22 - N[0]*crhs10 - N[0]*crhs8;
rhs[3]=DN(1,0)*crhs0*k - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1*rho + 1.0*N[1]*crhs12*crhs13*crhs14*rho - N[1]*crhs2 - N[1]*crhs6 + 1.0*crhs13*crhs14*crhs23*rho;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0*k - DN(1,1)*crhs11 - DN(1,1)*stress[1] + 1.0*N[1]*crhs12*crhs13*crhs19*rho + N[1]*crhs16*rho - N[1]*crhs17 - N[1]*crhs18 + 1.0*crhs13*crhs19*crhs23*rho;
rhs[5]=DN(1,0)*crhs21 + DN(1,1)*crhs22 - N[1]*crhs10 - N[1]*crhs8;
rhs[6]=DN(2,0)*crhs0*k - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1*rho + 1.0*N[2]*crhs12*crhs13*crhs14*rho - N[2]*crhs2 - N[2]*crhs6 + 1.0*crhs13*crhs14*crhs24*rho;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0*k - DN(2,1)*crhs11 - DN(2,1)*stress[1] + 1.0*N[2]*crhs12*crhs13*crhs19*rho + N[2]*crhs16*rho - N[2]*crhs17 - N[2]*crhs18 + 1.0*crhs13*crhs19*crhs24*rho;
rhs[8]=DN(2,0)*crhs21 + DN(2,1)*crhs22 - N[2]*crhs10 - N[2]*crhs8;


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

    const double crhs0 = N[0]*eps_vol[0] + N[1]*eps_vol[1] + N[2]*eps_vol[2] + N[3]*eps_vol[3];
const double crhs1 = N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0);
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs7 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs8 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2) + pow(crhs6, 2));
const double crhs9 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]) + N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]) + N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]) + N[3]*(bdf0*eps_vol[3] + bdf1*eps_vol_n[3] + bdf2*eps_vol_nn[3]);
const double crhs10 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs11 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs12 = crhs10 + crhs11 + crhs3;
const double crhs13 = k*(crhs12 + crhs9)*(crhs8*h/stab_c1 + mu);
const double crhs14 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs15 = 1.0/(crhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 = crhs1*rho - crhs2 - crhs7 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2] + DN(3,0)*eps_vol[3]);
const double crhs17 = DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6;
const double crhs18 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs19 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs20 = rho*(crhs10*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs21 = crhs18*rho - crhs19 - crhs20 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2] + DN(3,1)*eps_vol[3]);
const double crhs22 = N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2);
const double crhs23 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs24 = rho*(crhs11*crhs6 + crhs4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs25 = crhs22*rho - crhs23 - crhs24 + k*(DN(0,2)*eps_vol[0] + DN(1,2)*eps_vol[1] + DN(2,2)*eps_vol[2] + DN(3,2)*eps_vol[3]);
const double crhs26 = 1.0*crhs15;
const double crhs27 = crhs16*crhs26;
const double crhs28 = crhs21*crhs26;
const double crhs29 = crhs25*crhs26;
const double crhs30 = DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6;
const double crhs31 = DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6;
const double crhs32 = DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6;
rhs[0]=DN(0,0)*crhs0*k - DN(0,0)*crhs13 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1*rho + 1.0*N[0]*crhs14*crhs15*crhs16*rho - N[0]*crhs2 - N[0]*crhs7 + 1.0*crhs15*crhs16*crhs17*rho;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0*k - DN(0,1)*crhs13 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + 1.0*N[0]*crhs14*crhs15*crhs21*rho + N[0]*crhs18*rho - N[0]*crhs19 - N[0]*crhs20 + 1.0*crhs15*crhs17*crhs21*rho;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0*k - DN(0,2)*crhs13 - DN(0,2)*stress[2] + 1.0*N[0]*crhs14*crhs15*crhs25*rho + N[0]*crhs22*rho - N[0]*crhs23 - N[0]*crhs24 + 1.0*crhs15*crhs17*crhs25*rho;
rhs[3]=DN(0,0)*crhs27 + DN(0,1)*crhs28 + DN(0,2)*crhs29 - N[0]*crhs12 - N[0]*crhs9;
rhs[4]=DN(1,0)*crhs0*k - DN(1,0)*crhs13 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1*rho + 1.0*N[1]*crhs14*crhs15*crhs16*rho - N[1]*crhs2 - N[1]*crhs7 + 1.0*crhs15*crhs16*crhs30*rho;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0*k - DN(1,1)*crhs13 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + 1.0*N[1]*crhs14*crhs15*crhs21*rho + N[1]*crhs18*rho - N[1]*crhs19 - N[1]*crhs20 + 1.0*crhs15*crhs21*crhs30*rho;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0*k - DN(1,2)*crhs13 - DN(1,2)*stress[2] + 1.0*N[1]*crhs14*crhs15*crhs25*rho + N[1]*crhs22*rho - N[1]*crhs23 - N[1]*crhs24 + 1.0*crhs15*crhs25*crhs30*rho;
rhs[7]=DN(1,0)*crhs27 + DN(1,1)*crhs28 + DN(1,2)*crhs29 - N[1]*crhs12 - N[1]*crhs9;
rhs[8]=DN(2,0)*crhs0*k - DN(2,0)*crhs13 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1*rho + 1.0*N[2]*crhs14*crhs15*crhs16*rho - N[2]*crhs2 - N[2]*crhs7 + 1.0*crhs15*crhs16*crhs31*rho;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0*k - DN(2,1)*crhs13 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + 1.0*N[2]*crhs14*crhs15*crhs21*rho + N[2]*crhs18*rho - N[2]*crhs19 - N[2]*crhs20 + 1.0*crhs15*crhs21*crhs31*rho;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0*k - DN(2,2)*crhs13 - DN(2,2)*stress[2] + 1.0*N[2]*crhs14*crhs15*crhs25*rho + N[2]*crhs22*rho - N[2]*crhs23 - N[2]*crhs24 + 1.0*crhs15*crhs25*crhs31*rho;
rhs[11]=DN(2,0)*crhs27 + DN(2,1)*crhs28 + DN(2,2)*crhs29 - N[2]*crhs12 - N[2]*crhs9;
rhs[12]=DN(3,0)*crhs0*k - DN(3,0)*crhs13 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1*rho + 1.0*N[3]*crhs14*crhs15*crhs16*rho - N[3]*crhs2 - N[3]*crhs7 + 1.0*crhs15*crhs16*crhs32*rho;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0*k - DN(3,1)*crhs13 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + 1.0*N[3]*crhs14*crhs15*crhs21*rho + N[3]*crhs18*rho - N[3]*crhs19 - N[3]*crhs20 + 1.0*crhs15*crhs21*crhs32*rho;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0*k - DN(3,2)*crhs13 - DN(3,2)*stress[2] + 1.0*N[3]*crhs14*crhs15*crhs25*rho + N[3]*crhs22*rho - N[3]*crhs23 - N[3]*crhs24 + 1.0*crhs15*crhs25*crhs32*rho;
rhs[15]=DN(3,0)*crhs27 + DN(3,1)*crhs28 + DN(3,2)*crhs29 - N[3]*crhs12 - N[3]*crhs9;


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