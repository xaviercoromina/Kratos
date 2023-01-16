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
const double clhs11 = pow(N[0], 2)*bdf0;
const double clhs12 = N[0]*bdf0;
const double clhs13 = clhs12 + clhs9;
const double clhs14 = 1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 = clhs14*pow(rho, 2);
const double clhs16 = clhs15*clhs9;
const double clhs17 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs18 = clhs15*clhs17;
const double clhs19 = N[0]*clhs18;
const double clhs20 = clhs10*clhs9 + clhs11*rho + clhs13*clhs16 + clhs13*clhs19;
const double clhs21 = C(0,1)*DN(0,1) + clhs1;
const double clhs22 = C(1,2)*DN(0,1);
const double clhs23 = C(2,2)*DN(0,0) + clhs22;
const double clhs24 = DN(0,0)*clhs8;
const double clhs25 = DN(0,1)*clhs24;
const double clhs26 = clhs12*clhs7;
const double clhs27 = clhs14*clhs17;
const double clhs28 = clhs14*rho;
const double clhs29 = clhs28*clhs9;
const double clhs30 = k*(N[0] + clhs10*clhs27 + clhs26 + clhs29);
const double clhs31 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs32 = C(0,2)*DN(1,0);
const double clhs33 = C(2,2)*DN(1,1) + clhs32;
const double clhs34 = N[1]*clhs12;
const double clhs35 = clhs34*rho;
const double clhs36 = DN(1,0)*clhs24 + clhs35;
const double clhs37 = DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs38 = N[1]*bdf0;
const double clhs39 = clhs37 + clhs38;
const double clhs40 = clhs10*clhs37 + clhs16*clhs39 + clhs19*clhs39;
const double clhs41 = C(0,1)*DN(1,1) + clhs32;
const double clhs42 = C(1,2)*DN(1,1);
const double clhs43 = C(2,2)*DN(1,0) + clhs42;
const double clhs44 = DN(1,1)*clhs24;
const double clhs45 = DN(0,0)*N[1];
const double clhs46 = clhs38*clhs7;
const double clhs47 = DN(1,0)*N[0];
const double clhs48 = clhs17*clhs28;
const double clhs49 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs50 = C(0,2)*DN(2,0);
const double clhs51 = C(2,2)*DN(2,1) + clhs50;
const double clhs52 = N[2]*clhs12;
const double clhs53 = clhs52*rho;
const double clhs54 = DN(2,0)*clhs24 + clhs53;
const double clhs55 = DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs56 = N[2]*bdf0;
const double clhs57 = clhs55 + clhs56;
const double clhs58 = clhs10*clhs55 + clhs16*clhs57 + clhs19*clhs57;
const double clhs59 = C(0,1)*DN(2,1) + clhs50;
const double clhs60 = C(1,2)*DN(2,1);
const double clhs61 = C(2,2)*DN(2,0) + clhs60;
const double clhs62 = DN(2,1)*clhs24;
const double clhs63 = DN(0,0)*N[2];
const double clhs64 = clhs56*clhs7;
const double clhs65 = DN(2,0)*N[0];
const double clhs66 = C(0,1)*DN(0,0) + clhs22;
const double clhs67 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs68 = pow(DN(0,1), 2);
const double clhs69 = C(0,1)*DN(1,0) + clhs42;
const double clhs70 = DN(0,1)*clhs8;
const double clhs71 = DN(1,0)*clhs70;
const double clhs72 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs73 = DN(1,1)*clhs70 + clhs35;
const double clhs74 = DN(0,1)*N[1];
const double clhs75 = DN(1,1)*N[0];
const double clhs76 = C(0,1)*DN(2,0) + clhs60;
const double clhs77 = DN(2,0)*clhs70;
const double clhs78 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs79 = DN(2,1)*clhs70 + clhs53;
const double clhs80 = DN(0,1)*N[2];
const double clhs81 = DN(2,1)*N[0];
const double clhs82 = clhs13*clhs28;
const double clhs83 = N[0] + clhs82;
const double clhs84 = clhs14*k;
const double clhs85 = clhs28*clhs39;
const double clhs86 = DN(0,0)*clhs84;
const double clhs87 = DN(0,1)*clhs84;
const double clhs88 = -DN(1,0)*clhs86 - DN(1,1)*clhs87 - clhs34;
const double clhs89 = clhs28*clhs57;
const double clhs90 = -DN(2,0)*clhs86 - DN(2,1)*clhs87 - clhs52;
const double clhs91 = N[1]*rho;
const double clhs92 = clhs15*clhs37;
const double clhs93 = N[1]*clhs18;
const double clhs94 = clhs13*clhs92 + clhs13*clhs93 + clhs9*clhs91;
const double clhs95 = clhs28*clhs37;
const double clhs96 = pow(DN(1,0), 2);
const double clhs97 = pow(N[1], 2)*bdf0;
const double clhs98 = clhs37*clhs91 + clhs39*clhs92 + clhs39*clhs93 + clhs97*rho;
const double clhs99 = DN(1,0)*clhs8;
const double clhs100 = DN(1,1)*clhs99;
const double clhs101 = k*(N[1] + clhs27*clhs91 + clhs46 + clhs95);
const double clhs102 = N[2]*clhs38;
const double clhs103 = clhs102*rho;
const double clhs104 = DN(2,0)*clhs99 + clhs103;
const double clhs105 = clhs55*clhs91 + clhs57*clhs92 + clhs57*clhs93;
const double clhs106 = DN(2,1)*clhs99;
const double clhs107 = DN(1,0)*N[2];
const double clhs108 = DN(2,0)*N[1];
const double clhs109 = pow(DN(1,1), 2);
const double clhs110 = DN(1,1)*clhs8;
const double clhs111 = DN(2,0)*clhs110;
const double clhs112 = DN(2,1)*clhs110 + clhs103;
const double clhs113 = DN(1,1)*N[2];
const double clhs114 = DN(2,1)*N[1];
const double clhs115 = N[1] + clhs85;
const double clhs116 = -DN(1,0)*DN(2,0)*clhs84 - DN(1,1)*DN(2,1)*clhs84 - clhs102;
const double clhs117 = N[2]*rho;
const double clhs118 = clhs15*clhs55;
const double clhs119 = N[2]*clhs18;
const double clhs120 = clhs117*clhs9 + clhs118*clhs13 + clhs119*clhs13;
const double clhs121 = clhs28*clhs55;
const double clhs122 = clhs117*clhs37 + clhs118*clhs39 + clhs119*clhs39;
const double clhs123 = pow(DN(2,0), 2);
const double clhs124 = pow(N[2], 2)*bdf0;
const double clhs125 = clhs117*clhs55 + clhs118*clhs57 + clhs119*clhs57 + clhs124*rho;
const double clhs126 = DN(2,0)*DN(2,1)*clhs8;
const double clhs127 = k*(N[2] + clhs117*clhs27 + clhs121 + clhs64);
const double clhs128 = pow(DN(2,1), 2);
const double clhs129 = N[2] + clhs89;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs20 + clhs3*clhs8;
lhs(0,1)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs25;
lhs(0,2)=-DN(0,0)*clhs30;
lhs(0,3)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + clhs36 + clhs40;
lhs(0,4)=DN(0,0)*clhs41 + DN(0,1)*clhs43 + clhs44;
lhs(0,5)=-k*(DN(0,0)*clhs46 + DN(1,0)*clhs29 + clhs45 + clhs47*clhs48);
lhs(0,6)=DN(0,0)*clhs49 + DN(0,1)*clhs51 + clhs54 + clhs58;
lhs(0,7)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + clhs62;
lhs(0,8)=-k*(DN(0,0)*clhs64 + DN(2,0)*clhs29 + clhs48*clhs65 + clhs63);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs66 + clhs25;
lhs(1,1)=DN(0,0)*clhs23 + DN(0,1)*clhs67 + clhs20 + clhs68*clhs8;
lhs(1,2)=-DN(0,1)*clhs30;
lhs(1,3)=DN(0,0)*clhs33 + DN(0,1)*clhs69 + clhs71;
lhs(1,4)=DN(0,0)*clhs43 + DN(0,1)*clhs72 + clhs40 + clhs73;
lhs(1,5)=-k*(DN(0,1)*clhs46 + DN(1,1)*clhs29 + clhs48*clhs75 + clhs74);
lhs(1,6)=DN(0,0)*clhs51 + DN(0,1)*clhs76 + clhs77;
lhs(1,7)=DN(0,0)*clhs61 + DN(0,1)*clhs78 + clhs58 + clhs79;
lhs(1,8)=-k*(DN(0,1)*clhs64 + DN(2,1)*clhs29 + clhs48*clhs81 + clhs80);
lhs(2,0)=DN(0,0)*clhs83;
lhs(2,1)=DN(0,1)*clhs83;
lhs(2,2)=-clhs11 - clhs3*clhs84 - clhs68*clhs84;
lhs(2,3)=DN(0,0)*clhs85 + clhs47;
lhs(2,4)=DN(0,1)*clhs85 + clhs75;
lhs(2,5)=clhs88;
lhs(2,6)=DN(0,0)*clhs89 + clhs65;
lhs(2,7)=DN(0,1)*clhs89 + clhs81;
lhs(2,8)=clhs90;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs36 + clhs94;
lhs(3,1)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs71;
lhs(3,2)=-k*(DN(0,0)*clhs95 + DN(1,0)*clhs26 + clhs45*clhs48 + clhs47);
lhs(3,3)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + clhs8*clhs96 + clhs98;
lhs(3,4)=DN(1,0)*clhs41 + DN(1,1)*clhs43 + clhs100;
lhs(3,5)=-DN(1,0)*clhs101;
lhs(3,6)=DN(1,0)*clhs49 + DN(1,1)*clhs51 + clhs104 + clhs105;
lhs(3,7)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + clhs106;
lhs(3,8)=-k*(DN(1,0)*clhs64 + DN(2,0)*clhs95 + clhs107 + clhs108*clhs48);
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs66 + clhs44;
lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs67 + clhs73 + clhs94;
lhs(4,2)=-k*(DN(0,1)*clhs95 + DN(1,1)*clhs26 + clhs48*clhs74 + clhs75);
lhs(4,3)=DN(1,0)*clhs33 + DN(1,1)*clhs69 + clhs100;
lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs72 + clhs109*clhs8 + clhs98;
lhs(4,5)=-DN(1,1)*clhs101;
lhs(4,6)=DN(1,0)*clhs51 + DN(1,1)*clhs76 + clhs111;
lhs(4,7)=DN(1,0)*clhs61 + DN(1,1)*clhs78 + clhs105 + clhs112;
lhs(4,8)=-k*(DN(1,1)*clhs64 + DN(2,1)*clhs95 + clhs113 + clhs114*clhs48);
lhs(5,0)=DN(1,0)*clhs82 + clhs45;
lhs(5,1)=DN(1,1)*clhs82 + clhs74;
lhs(5,2)=clhs88;
lhs(5,3)=DN(1,0)*clhs115;
lhs(5,4)=DN(1,1)*clhs115;
lhs(5,5)=-clhs109*clhs84 - clhs84*clhs96 - clhs97;
lhs(5,6)=DN(1,0)*clhs89 + clhs108;
lhs(5,7)=DN(1,1)*clhs89 + clhs114;
lhs(5,8)=clhs116;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs120 + clhs54;
lhs(6,1)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs77;
lhs(6,2)=-k*(DN(0,0)*clhs121 + DN(2,0)*clhs26 + clhs48*clhs63 + clhs65);
lhs(6,3)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + clhs104 + clhs122;
lhs(6,4)=DN(2,0)*clhs41 + DN(2,1)*clhs43 + clhs111;
lhs(6,5)=-k*(DN(1,0)*clhs121 + DN(2,0)*clhs46 + clhs107*clhs48 + clhs108);
lhs(6,6)=DN(2,0)*clhs49 + DN(2,1)*clhs51 + clhs123*clhs8 + clhs125;
lhs(6,7)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + clhs126;
lhs(6,8)=-DN(2,0)*clhs127;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs66 + clhs62;
lhs(7,1)=DN(2,0)*clhs23 + DN(2,1)*clhs67 + clhs120 + clhs79;
lhs(7,2)=-k*(DN(0,1)*clhs121 + DN(2,1)*clhs26 + clhs48*clhs80 + clhs81);
lhs(7,3)=DN(2,0)*clhs33 + DN(2,1)*clhs69 + clhs106;
lhs(7,4)=DN(2,0)*clhs43 + DN(2,1)*clhs72 + clhs112 + clhs122;
lhs(7,5)=-k*(DN(1,1)*clhs121 + DN(2,1)*clhs46 + clhs113*clhs48 + clhs114);
lhs(7,6)=DN(2,0)*clhs51 + DN(2,1)*clhs76 + clhs126;
lhs(7,7)=DN(2,0)*clhs61 + DN(2,1)*clhs78 + clhs125 + clhs128*clhs8;
lhs(7,8)=-DN(2,1)*clhs127;
lhs(8,0)=DN(2,0)*clhs82 + clhs63;
lhs(8,1)=DN(2,1)*clhs82 + clhs80;
lhs(8,2)=clhs90;
lhs(8,3)=DN(2,0)*clhs85 + clhs107;
lhs(8,4)=DN(2,1)*clhs85 + clhs113;
lhs(8,5)=clhs116;
lhs(8,6)=DN(2,0)*clhs129;
lhs(8,7)=DN(2,1)*clhs129;
lhs(8,8)=-clhs123*clhs84 - clhs124 - clhs128*clhs84;


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
const double clhs14 = pow(N[0], 2)*bdf0;
const double clhs15 = N[0]*bdf0;
const double clhs16 = clhs12 + clhs15;
const double clhs17 = 1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = clhs17*pow(rho, 2);
const double clhs19 = clhs12*clhs18;
const double clhs20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs21 = clhs18*clhs20;
const double clhs22 = N[0]*clhs21;
const double clhs23 = clhs12*clhs13 + clhs14*rho + clhs16*clhs19 + clhs16*clhs22;
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
const double clhs38 = clhs10*clhs15;
const double clhs39 = clhs17*clhs20;
const double clhs40 = clhs17*rho;
const double clhs41 = clhs12*clhs40;
const double clhs42 = k*(N[0] + clhs13*clhs39 + clhs38 + clhs41);
const double clhs43 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs44 = C(0,3)*DN(1,0);
const double clhs45 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs44;
const double clhs46 = C(0,5)*DN(1,0);
const double clhs47 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs46;
const double clhs48 = N[1]*clhs15;
const double clhs49 = clhs48*rho;
const double clhs50 = DN(1,0)*clhs30 + clhs49;
const double clhs51 = DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs52 = N[1]*bdf0;
const double clhs53 = clhs51 + clhs52;
const double clhs54 = clhs13*clhs51 + clhs19*clhs53 + clhs22*clhs53;
const double clhs55 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs44;
const double clhs56 = C(1,3)*DN(1,1);
const double clhs57 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs56;
const double clhs58 = C(3,5)*DN(1,0);
const double clhs59 = C(4,5)*DN(1,2);
const double clhs60 = C(1,5)*DN(1,1) + clhs58 + clhs59;
const double clhs61 = DN(1,1)*clhs30;
const double clhs62 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs46;
const double clhs63 = C(3,4)*DN(1,1);
const double clhs64 = C(2,3)*DN(1,2) + clhs58 + clhs63;
const double clhs65 = C(2,5)*DN(1,2);
const double clhs66 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs65;
const double clhs67 = DN(1,2)*clhs30;
const double clhs68 = DN(0,0)*N[1];
const double clhs69 = clhs10*clhs52;
const double clhs70 = DN(1,0)*N[0];
const double clhs71 = clhs20*clhs40;
const double clhs72 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs73 = C(0,3)*DN(2,0);
const double clhs74 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs73;
const double clhs75 = C(0,5)*DN(2,0);
const double clhs76 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs75;
const double clhs77 = N[2]*clhs15;
const double clhs78 = clhs77*rho;
const double clhs79 = DN(2,0)*clhs30 + clhs78;
const double clhs80 = DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs81 = N[2]*bdf0;
const double clhs82 = clhs80 + clhs81;
const double clhs83 = clhs13*clhs80 + clhs19*clhs82 + clhs22*clhs82;
const double clhs84 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs73;
const double clhs85 = C(1,3)*DN(2,1);
const double clhs86 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs85;
const double clhs87 = C(3,5)*DN(2,0);
const double clhs88 = C(4,5)*DN(2,2);
const double clhs89 = C(1,5)*DN(2,1) + clhs87 + clhs88;
const double clhs90 = DN(2,1)*clhs30;
const double clhs91 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs75;
const double clhs92 = C(3,4)*DN(2,1);
const double clhs93 = C(2,3)*DN(2,2) + clhs87 + clhs92;
const double clhs94 = C(2,5)*DN(2,2);
const double clhs95 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs94;
const double clhs96 = DN(2,2)*clhs30;
const double clhs97 = DN(0,0)*N[2];
const double clhs98 = clhs10*clhs81;
const double clhs99 = DN(2,0)*N[0];
const double clhs100 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs101 = C(0,3)*DN(3,0);
const double clhs102 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs101;
const double clhs103 = C(0,5)*DN(3,0);
const double clhs104 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs103;
const double clhs105 = N[3]*clhs15;
const double clhs106 = clhs105*rho;
const double clhs107 = DN(3,0)*clhs30 + clhs106;
const double clhs108 = DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs109 = N[3]*bdf0;
const double clhs110 = clhs108 + clhs109;
const double clhs111 = clhs108*clhs13 + clhs110*clhs19 + clhs110*clhs22;
const double clhs112 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs101;
const double clhs113 = C(1,3)*DN(3,1);
const double clhs114 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs113;
const double clhs115 = C(3,5)*DN(3,0);
const double clhs116 = C(4,5)*DN(3,2);
const double clhs117 = C(1,5)*DN(3,1) + clhs115 + clhs116;
const double clhs118 = DN(3,1)*clhs30;
const double clhs119 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs103;
const double clhs120 = C(3,4)*DN(3,1);
const double clhs121 = C(2,3)*DN(3,2) + clhs115 + clhs120;
const double clhs122 = C(2,5)*DN(3,2);
const double clhs123 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs122;
const double clhs124 = DN(3,2)*clhs30;
const double clhs125 = DN(0,0)*N[3];
const double clhs126 = clhs10*clhs109;
const double clhs127 = DN(3,0)*N[0];
const double clhs128 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs25;
const double clhs129 = C(0,4)*DN(0,0) + clhs28 + clhs33;
const double clhs130 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs131 = C(1,4)*DN(0,1);
const double clhs132 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs131;
const double clhs133 = pow(DN(0,1), 2);
const double clhs134 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs131;
const double clhs135 = C(2,4)*DN(0,2);
const double clhs136 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs135;
const double clhs137 = DN(0,1)*clhs11;
const double clhs138 = DN(0,2)*clhs137;
const double clhs139 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs56;
const double clhs140 = C(0,4)*DN(1,0) + clhs59 + clhs63;
const double clhs141 = DN(1,0)*clhs137;
const double clhs142 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs143 = C(1,4)*DN(1,1);
const double clhs144 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs143;
const double clhs145 = DN(1,1)*clhs137;
const double clhs146 = clhs49 + clhs54;
const double clhs147 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs143;
const double clhs148 = C(2,4)*DN(1,2);
const double clhs149 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs148;
const double clhs150 = DN(1,2)*clhs137;
const double clhs151 = DN(0,1)*N[1];
const double clhs152 = DN(1,1)*N[0];
const double clhs153 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs85;
const double clhs154 = C(0,4)*DN(2,0) + clhs88 + clhs92;
const double clhs155 = DN(2,0)*clhs137;
const double clhs156 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs157 = C(1,4)*DN(2,1);
const double clhs158 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs157;
const double clhs159 = DN(2,1)*clhs137;
const double clhs160 = clhs78 + clhs83;
const double clhs161 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs157;
const double clhs162 = C(2,4)*DN(2,2);
const double clhs163 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs162;
const double clhs164 = DN(2,2)*clhs137;
const double clhs165 = DN(0,1)*N[2];
const double clhs166 = DN(2,1)*N[0];
const double clhs167 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs113;
const double clhs168 = C(0,4)*DN(3,0) + clhs116 + clhs120;
const double clhs169 = DN(3,0)*clhs137;
const double clhs170 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs171 = C(1,4)*DN(3,1);
const double clhs172 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs171;
const double clhs173 = DN(3,1)*clhs137;
const double clhs174 = clhs106 + clhs111;
const double clhs175 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs171;
const double clhs176 = C(2,4)*DN(3,2);
const double clhs177 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs176;
const double clhs178 = DN(3,2)*clhs137;
const double clhs179 = DN(0,1)*N[3];
const double clhs180 = DN(3,1)*N[0];
const double clhs181 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
const double clhs182 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs135;
const double clhs183 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs184 = pow(DN(0,2), 2);
const double clhs185 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs65;
const double clhs186 = DN(0,2)*clhs11;
const double clhs187 = DN(1,0)*clhs186;
const double clhs188 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs148;
const double clhs189 = DN(1,1)*clhs186;
const double clhs190 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs191 = DN(1,2)*clhs186;
const double clhs192 = DN(0,2)*N[1];
const double clhs193 = DN(1,2)*N[0];
const double clhs194 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs94;
const double clhs195 = DN(2,0)*clhs186;
const double clhs196 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs162;
const double clhs197 = DN(2,1)*clhs186;
const double clhs198 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs199 = DN(2,2)*clhs186;
const double clhs200 = DN(0,2)*N[2];
const double clhs201 = DN(2,2)*N[0];
const double clhs202 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs122;
const double clhs203 = DN(3,0)*clhs186;
const double clhs204 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs176;
const double clhs205 = DN(3,1)*clhs186;
const double clhs206 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs207 = DN(3,2)*clhs186;
const double clhs208 = DN(0,2)*N[3];
const double clhs209 = DN(3,2)*N[0];
const double clhs210 = clhs16*clhs40;
const double clhs211 = N[0] + clhs210;
const double clhs212 = clhs17*k;
const double clhs213 = clhs40*clhs53;
const double clhs214 = DN(0,0)*clhs212;
const double clhs215 = DN(0,1)*clhs212;
const double clhs216 = DN(0,2)*clhs212;
const double clhs217 = -DN(1,0)*clhs214 - DN(1,1)*clhs215 - DN(1,2)*clhs216 - clhs48;
const double clhs218 = clhs40*clhs82;
const double clhs219 = -DN(2,0)*clhs214 - DN(2,1)*clhs215 - DN(2,2)*clhs216 - clhs77;
const double clhs220 = clhs110*clhs40;
const double clhs221 = -DN(3,0)*clhs214 - DN(3,1)*clhs215 - DN(3,2)*clhs216 - clhs105;
const double clhs222 = N[1]*rho;
const double clhs223 = clhs18*clhs51;
const double clhs224 = N[1]*clhs21;
const double clhs225 = clhs12*clhs222 + clhs16*clhs223 + clhs16*clhs224;
const double clhs226 = clhs40*clhs51;
const double clhs227 = pow(DN(1,0), 2);
const double clhs228 = pow(N[1], 2)*bdf0;
const double clhs229 = clhs222*clhs51 + clhs223*clhs53 + clhs224*clhs53 + clhs228*rho;
const double clhs230 = DN(1,0)*clhs11;
const double clhs231 = DN(1,1)*clhs230;
const double clhs232 = DN(1,2)*clhs230;
const double clhs233 = k*(N[1] + clhs222*clhs39 + clhs226 + clhs69);
const double clhs234 = N[2]*clhs52;
const double clhs235 = clhs234*rho;
const double clhs236 = DN(2,0)*clhs230 + clhs235;
const double clhs237 = clhs222*clhs80 + clhs223*clhs82 + clhs224*clhs82;
const double clhs238 = DN(2,1)*clhs230;
const double clhs239 = DN(2,2)*clhs230;
const double clhs240 = DN(1,0)*N[2];
const double clhs241 = DN(2,0)*N[1];
const double clhs242 = N[3]*clhs52;
const double clhs243 = clhs242*rho;
const double clhs244 = DN(3,0)*clhs230 + clhs243;
const double clhs245 = clhs108*clhs222 + clhs110*clhs223 + clhs110*clhs224;
const double clhs246 = DN(3,1)*clhs230;
const double clhs247 = DN(3,2)*clhs230;
const double clhs248 = DN(1,0)*N[3];
const double clhs249 = DN(3,0)*N[1];
const double clhs250 = clhs225 + clhs49;
const double clhs251 = pow(DN(1,1), 2);
const double clhs252 = DN(1,1)*clhs11;
const double clhs253 = DN(1,2)*clhs252;
const double clhs254 = DN(2,0)*clhs252;
const double clhs255 = DN(2,1)*clhs252;
const double clhs256 = clhs235 + clhs237;
const double clhs257 = DN(2,2)*clhs252;
const double clhs258 = DN(1,1)*N[2];
const double clhs259 = DN(2,1)*N[1];
const double clhs260 = DN(3,0)*clhs252;
const double clhs261 = DN(3,1)*clhs252;
const double clhs262 = clhs243 + clhs245;
const double clhs263 = DN(3,2)*clhs252;
const double clhs264 = DN(1,1)*N[3];
const double clhs265 = DN(3,1)*N[1];
const double clhs266 = pow(DN(1,2), 2);
const double clhs267 = DN(1,2)*clhs11;
const double clhs268 = DN(2,0)*clhs267;
const double clhs269 = DN(2,1)*clhs267;
const double clhs270 = DN(2,2)*clhs267;
const double clhs271 = DN(1,2)*N[2];
const double clhs272 = DN(2,2)*N[1];
const double clhs273 = DN(3,0)*clhs267;
const double clhs274 = DN(3,1)*clhs267;
const double clhs275 = DN(3,2)*clhs267;
const double clhs276 = DN(1,2)*N[3];
const double clhs277 = DN(3,2)*N[1];
const double clhs278 = N[1] + clhs213;
const double clhs279 = DN(1,0)*clhs212;
const double clhs280 = DN(1,1)*clhs212;
const double clhs281 = DN(1,2)*clhs212;
const double clhs282 = -DN(2,0)*clhs279 - DN(2,1)*clhs280 - DN(2,2)*clhs281 - clhs234;
const double clhs283 = -DN(3,0)*clhs279 - DN(3,1)*clhs280 - DN(3,2)*clhs281 - clhs242;
const double clhs284 = N[2]*rho;
const double clhs285 = clhs18*clhs80;
const double clhs286 = N[2]*clhs21;
const double clhs287 = clhs12*clhs284 + clhs16*clhs285 + clhs16*clhs286;
const double clhs288 = clhs40*clhs80;
const double clhs289 = clhs284*clhs51 + clhs285*clhs53 + clhs286*clhs53;
const double clhs290 = pow(DN(2,0), 2);
const double clhs291 = pow(N[2], 2)*bdf0;
const double clhs292 = clhs284*clhs80 + clhs285*clhs82 + clhs286*clhs82 + clhs291*rho;
const double clhs293 = DN(2,0)*clhs11;
const double clhs294 = DN(2,1)*clhs293;
const double clhs295 = DN(2,2)*clhs293;
const double clhs296 = k*(N[2] + clhs284*clhs39 + clhs288 + clhs98);
const double clhs297 = N[3]*clhs81;
const double clhs298 = clhs297*rho;
const double clhs299 = DN(3,0)*clhs293 + clhs298;
const double clhs300 = clhs108*clhs284 + clhs110*clhs285 + clhs110*clhs286;
const double clhs301 = DN(3,1)*clhs293;
const double clhs302 = DN(3,2)*clhs293;
const double clhs303 = DN(2,0)*N[3];
const double clhs304 = DN(3,0)*N[2];
const double clhs305 = clhs287 + clhs78;
const double clhs306 = clhs235 + clhs289;
const double clhs307 = pow(DN(2,1), 2);
const double clhs308 = DN(2,1)*clhs11;
const double clhs309 = DN(2,2)*clhs308;
const double clhs310 = DN(3,0)*clhs308;
const double clhs311 = DN(3,1)*clhs308;
const double clhs312 = clhs298 + clhs300;
const double clhs313 = DN(3,2)*clhs308;
const double clhs314 = DN(2,1)*N[3];
const double clhs315 = DN(3,1)*N[2];
const double clhs316 = pow(DN(2,2), 2);
const double clhs317 = DN(2,2)*clhs11;
const double clhs318 = DN(3,0)*clhs317;
const double clhs319 = DN(3,1)*clhs317;
const double clhs320 = DN(3,2)*clhs317;
const double clhs321 = DN(2,2)*N[3];
const double clhs322 = DN(3,2)*N[2];
const double clhs323 = N[2] + clhs218;
const double clhs324 = -DN(2,0)*DN(3,0)*clhs212 - DN(2,1)*DN(3,1)*clhs212 - DN(2,2)*DN(3,2)*clhs212 - clhs297;
const double clhs325 = N[3]*rho;
const double clhs326 = clhs108*clhs18;
const double clhs327 = N[3]*clhs21;
const double clhs328 = clhs12*clhs325 + clhs16*clhs326 + clhs16*clhs327;
const double clhs329 = clhs108*clhs40;
const double clhs330 = clhs325*clhs51 + clhs326*clhs53 + clhs327*clhs53;
const double clhs331 = clhs325*clhs80 + clhs326*clhs82 + clhs327*clhs82;
const double clhs332 = pow(DN(3,0), 2);
const double clhs333 = pow(N[3], 2)*bdf0;
const double clhs334 = clhs108*clhs325 + clhs110*clhs326 + clhs110*clhs327 + clhs333*rho;
const double clhs335 = DN(3,0)*clhs11;
const double clhs336 = DN(3,1)*clhs335;
const double clhs337 = DN(3,2)*clhs335;
const double clhs338 = k*(N[3] + clhs126 + clhs325*clhs39 + clhs329);
const double clhs339 = clhs106 + clhs328;
const double clhs340 = clhs243 + clhs330;
const double clhs341 = clhs298 + clhs331;
const double clhs342 = pow(DN(3,1), 2);
const double clhs343 = DN(3,1)*DN(3,2)*clhs11;
const double clhs344 = pow(DN(3,2), 2);
const double clhs345 = N[3] + clhs220;
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs11*clhs5 + clhs23;
lhs(0,1)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + DN(0,2)*clhs29 + clhs31;
lhs(0,2)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs37;
lhs(0,3)=-DN(0,0)*clhs42;
lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + DN(0,2)*clhs47 + clhs50 + clhs54;
lhs(0,5)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs61;
lhs(0,6)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs67;
lhs(0,7)=-k*(DN(0,0)*clhs69 + DN(1,0)*clhs41 + clhs68 + clhs70*clhs71);
lhs(0,8)=DN(0,0)*clhs72 + DN(0,1)*clhs74 + DN(0,2)*clhs76 + clhs79 + clhs83;
lhs(0,9)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs89 + clhs90;
lhs(0,10)=DN(0,0)*clhs91 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs96;
lhs(0,11)=-k*(DN(0,0)*clhs98 + DN(2,0)*clhs41 + clhs71*clhs99 + clhs97);
lhs(0,12)=DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs107 + clhs111;
lhs(0,13)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs117 + clhs118;
lhs(0,14)=DN(0,0)*clhs119 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs124;
lhs(0,15)=-k*(DN(0,0)*clhs126 + DN(3,0)*clhs41 + clhs125 + clhs127*clhs71);
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs128 + DN(0,2)*clhs129 + clhs31;
lhs(1,1)=DN(0,0)*clhs26 + DN(0,1)*clhs130 + DN(0,2)*clhs132 + clhs11*clhs133 + clhs23;
lhs(1,2)=DN(0,0)*clhs34 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs138;
lhs(1,3)=-DN(0,1)*clhs42;
lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs139 + DN(0,2)*clhs140 + clhs141;
lhs(1,5)=DN(0,0)*clhs57 + DN(0,1)*clhs142 + DN(0,2)*clhs144 + clhs145 + clhs146;
lhs(1,6)=DN(0,0)*clhs64 + DN(0,1)*clhs147 + DN(0,2)*clhs149 + clhs150;
lhs(1,7)=-k*(DN(0,1)*clhs69 + DN(1,1)*clhs41 + clhs151 + clhs152*clhs71);
lhs(1,8)=DN(0,0)*clhs74 + DN(0,1)*clhs153 + DN(0,2)*clhs154 + clhs155;
lhs(1,9)=DN(0,0)*clhs86 + DN(0,1)*clhs156 + DN(0,2)*clhs158 + clhs159 + clhs160;
lhs(1,10)=DN(0,0)*clhs93 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs164;
lhs(1,11)=-k*(DN(0,1)*clhs98 + DN(2,1)*clhs41 + clhs165 + clhs166*clhs71);
lhs(1,12)=DN(0,0)*clhs102 + DN(0,1)*clhs167 + DN(0,2)*clhs168 + clhs169;
lhs(1,13)=DN(0,0)*clhs114 + DN(0,1)*clhs170 + DN(0,2)*clhs172 + clhs173 + clhs174;
lhs(1,14)=DN(0,0)*clhs121 + DN(0,1)*clhs175 + DN(0,2)*clhs177 + clhs178;
lhs(1,15)=-k*(DN(0,1)*clhs126 + DN(3,1)*clhs41 + clhs179 + clhs180*clhs71);
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs129 + DN(0,2)*clhs181 + clhs37;
lhs(2,1)=DN(0,0)*clhs29 + DN(0,1)*clhs132 + DN(0,2)*clhs182 + clhs138;
lhs(2,2)=DN(0,0)*clhs36 + DN(0,1)*clhs136 + DN(0,2)*clhs183 + clhs11*clhs184 + clhs23;
lhs(2,3)=-DN(0,2)*clhs42;
lhs(2,4)=DN(0,0)*clhs47 + DN(0,1)*clhs140 + DN(0,2)*clhs185 + clhs187;
lhs(2,5)=DN(0,0)*clhs60 + DN(0,1)*clhs144 + DN(0,2)*clhs188 + clhs189;
lhs(2,6)=DN(0,0)*clhs66 + DN(0,1)*clhs149 + DN(0,2)*clhs190 + clhs146 + clhs191;
lhs(2,7)=-k*(DN(0,2)*clhs69 + DN(1,2)*clhs41 + clhs192 + clhs193*clhs71);
lhs(2,8)=DN(0,0)*clhs76 + DN(0,1)*clhs154 + DN(0,2)*clhs194 + clhs195;
lhs(2,9)=DN(0,0)*clhs89 + DN(0,1)*clhs158 + DN(0,2)*clhs196 + clhs197;
lhs(2,10)=DN(0,0)*clhs95 + DN(0,1)*clhs163 + DN(0,2)*clhs198 + clhs160 + clhs199;
lhs(2,11)=-k*(DN(0,2)*clhs98 + DN(2,2)*clhs41 + clhs200 + clhs201*clhs71);
lhs(2,12)=DN(0,0)*clhs104 + DN(0,1)*clhs168 + DN(0,2)*clhs202 + clhs203;
lhs(2,13)=DN(0,0)*clhs117 + DN(0,1)*clhs172 + DN(0,2)*clhs204 + clhs205;
lhs(2,14)=DN(0,0)*clhs123 + DN(0,1)*clhs177 + DN(0,2)*clhs206 + clhs174 + clhs207;
lhs(2,15)=-k*(DN(0,2)*clhs126 + DN(3,2)*clhs41 + clhs208 + clhs209*clhs71);
lhs(3,0)=DN(0,0)*clhs211;
lhs(3,1)=DN(0,1)*clhs211;
lhs(3,2)=DN(0,2)*clhs211;
lhs(3,3)=-clhs133*clhs212 - clhs14 - clhs184*clhs212 - clhs212*clhs5;
lhs(3,4)=DN(0,0)*clhs213 + clhs70;
lhs(3,5)=DN(0,1)*clhs213 + clhs152;
lhs(3,6)=DN(0,2)*clhs213 + clhs193;
lhs(3,7)=clhs217;
lhs(3,8)=DN(0,0)*clhs218 + clhs99;
lhs(3,9)=DN(0,1)*clhs218 + clhs166;
lhs(3,10)=DN(0,2)*clhs218 + clhs201;
lhs(3,11)=clhs219;
lhs(3,12)=DN(0,0)*clhs220 + clhs127;
lhs(3,13)=DN(0,1)*clhs220 + clhs180;
lhs(3,14)=DN(0,2)*clhs220 + clhs209;
lhs(3,15)=clhs221;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs225 + clhs50;
lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + DN(1,2)*clhs29 + clhs141;
lhs(4,2)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs187;
lhs(4,3)=-k*(DN(0,0)*clhs226 + DN(1,0)*clhs38 + clhs68*clhs71 + clhs70);
lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + DN(1,2)*clhs47 + clhs11*clhs227 + clhs229;
lhs(4,5)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs231;
lhs(4,6)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs232;
lhs(4,7)=-DN(1,0)*clhs233;
lhs(4,8)=DN(1,0)*clhs72 + DN(1,1)*clhs74 + DN(1,2)*clhs76 + clhs236 + clhs237;
lhs(4,9)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs89 + clhs238;
lhs(4,10)=DN(1,0)*clhs91 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs239;
lhs(4,11)=-k*(DN(1,0)*clhs98 + DN(2,0)*clhs226 + clhs240 + clhs241*clhs71);
lhs(4,12)=DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs244 + clhs245;
lhs(4,13)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs117 + clhs246;
lhs(4,14)=DN(1,0)*clhs119 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs247;
lhs(4,15)=-k*(DN(1,0)*clhs126 + DN(3,0)*clhs226 + clhs248 + clhs249*clhs71);
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs128 + DN(1,2)*clhs129 + clhs61;
lhs(5,1)=DN(1,0)*clhs26 + DN(1,1)*clhs130 + DN(1,2)*clhs132 + clhs145 + clhs250;
lhs(5,2)=DN(1,0)*clhs34 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs189;
lhs(5,3)=-k*(DN(0,1)*clhs226 + DN(1,1)*clhs38 + clhs151*clhs71 + clhs152);
lhs(5,4)=DN(1,0)*clhs45 + DN(1,1)*clhs139 + DN(1,2)*clhs140 + clhs231;
lhs(5,5)=DN(1,0)*clhs57 + DN(1,1)*clhs142 + DN(1,2)*clhs144 + clhs11*clhs251 + clhs229;
lhs(5,6)=DN(1,0)*clhs64 + DN(1,1)*clhs147 + DN(1,2)*clhs149 + clhs253;
lhs(5,7)=-DN(1,1)*clhs233;
lhs(5,8)=DN(1,0)*clhs74 + DN(1,1)*clhs153 + DN(1,2)*clhs154 + clhs254;
lhs(5,9)=DN(1,0)*clhs86 + DN(1,1)*clhs156 + DN(1,2)*clhs158 + clhs255 + clhs256;
lhs(5,10)=DN(1,0)*clhs93 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs257;
lhs(5,11)=-k*(DN(1,1)*clhs98 + DN(2,1)*clhs226 + clhs258 + clhs259*clhs71);
lhs(5,12)=DN(1,0)*clhs102 + DN(1,1)*clhs167 + DN(1,2)*clhs168 + clhs260;
lhs(5,13)=DN(1,0)*clhs114 + DN(1,1)*clhs170 + DN(1,2)*clhs172 + clhs261 + clhs262;
lhs(5,14)=DN(1,0)*clhs121 + DN(1,1)*clhs175 + DN(1,2)*clhs177 + clhs263;
lhs(5,15)=-k*(DN(1,1)*clhs126 + DN(3,1)*clhs226 + clhs264 + clhs265*clhs71);
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs129 + DN(1,2)*clhs181 + clhs67;
lhs(6,1)=DN(1,0)*clhs29 + DN(1,1)*clhs132 + DN(1,2)*clhs182 + clhs150;
lhs(6,2)=DN(1,0)*clhs36 + DN(1,1)*clhs136 + DN(1,2)*clhs183 + clhs191 + clhs250;
lhs(6,3)=-k*(DN(0,2)*clhs226 + DN(1,2)*clhs38 + clhs192*clhs71 + clhs193);
lhs(6,4)=DN(1,0)*clhs47 + DN(1,1)*clhs140 + DN(1,2)*clhs185 + clhs232;
lhs(6,5)=DN(1,0)*clhs60 + DN(1,1)*clhs144 + DN(1,2)*clhs188 + clhs253;
lhs(6,6)=DN(1,0)*clhs66 + DN(1,1)*clhs149 + DN(1,2)*clhs190 + clhs11*clhs266 + clhs229;
lhs(6,7)=-DN(1,2)*clhs233;
lhs(6,8)=DN(1,0)*clhs76 + DN(1,1)*clhs154 + DN(1,2)*clhs194 + clhs268;
lhs(6,9)=DN(1,0)*clhs89 + DN(1,1)*clhs158 + DN(1,2)*clhs196 + clhs269;
lhs(6,10)=DN(1,0)*clhs95 + DN(1,1)*clhs163 + DN(1,2)*clhs198 + clhs256 + clhs270;
lhs(6,11)=-k*(DN(1,2)*clhs98 + DN(2,2)*clhs226 + clhs271 + clhs272*clhs71);
lhs(6,12)=DN(1,0)*clhs104 + DN(1,1)*clhs168 + DN(1,2)*clhs202 + clhs273;
lhs(6,13)=DN(1,0)*clhs117 + DN(1,1)*clhs172 + DN(1,2)*clhs204 + clhs274;
lhs(6,14)=DN(1,0)*clhs123 + DN(1,1)*clhs177 + DN(1,2)*clhs206 + clhs262 + clhs275;
lhs(6,15)=-k*(DN(1,2)*clhs126 + DN(3,2)*clhs226 + clhs276 + clhs277*clhs71);
lhs(7,0)=DN(1,0)*clhs210 + clhs68;
lhs(7,1)=DN(1,1)*clhs210 + clhs151;
lhs(7,2)=DN(1,2)*clhs210 + clhs192;
lhs(7,3)=clhs217;
lhs(7,4)=DN(1,0)*clhs278;
lhs(7,5)=DN(1,1)*clhs278;
lhs(7,6)=DN(1,2)*clhs278;
lhs(7,7)=-clhs212*clhs227 - clhs212*clhs251 - clhs212*clhs266 - clhs228;
lhs(7,8)=DN(1,0)*clhs218 + clhs241;
lhs(7,9)=DN(1,1)*clhs218 + clhs259;
lhs(7,10)=DN(1,2)*clhs218 + clhs272;
lhs(7,11)=clhs282;
lhs(7,12)=DN(1,0)*clhs220 + clhs249;
lhs(7,13)=DN(1,1)*clhs220 + clhs265;
lhs(7,14)=DN(1,2)*clhs220 + clhs277;
lhs(7,15)=clhs283;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs287 + clhs79;
lhs(8,1)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + DN(2,2)*clhs29 + clhs155;
lhs(8,2)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs195;
lhs(8,3)=-k*(DN(0,0)*clhs288 + DN(2,0)*clhs38 + clhs71*clhs97 + clhs99);
lhs(8,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + DN(2,2)*clhs47 + clhs236 + clhs289;
lhs(8,5)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs254;
lhs(8,6)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs268;
lhs(8,7)=-k*(DN(1,0)*clhs288 + DN(2,0)*clhs69 + clhs240*clhs71 + clhs241);
lhs(8,8)=DN(2,0)*clhs72 + DN(2,1)*clhs74 + DN(2,2)*clhs76 + clhs11*clhs290 + clhs292;
lhs(8,9)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs89 + clhs294;
lhs(8,10)=DN(2,0)*clhs91 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs295;
lhs(8,11)=-DN(2,0)*clhs296;
lhs(8,12)=DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs299 + clhs300;
lhs(8,13)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs117 + clhs301;
lhs(8,14)=DN(2,0)*clhs119 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs302;
lhs(8,15)=-k*(DN(2,0)*clhs126 + DN(3,0)*clhs288 + clhs303 + clhs304*clhs71);
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs128 + DN(2,2)*clhs129 + clhs90;
lhs(9,1)=DN(2,0)*clhs26 + DN(2,1)*clhs130 + DN(2,2)*clhs132 + clhs159 + clhs305;
lhs(9,2)=DN(2,0)*clhs34 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs197;
lhs(9,3)=-k*(DN(0,1)*clhs288 + DN(2,1)*clhs38 + clhs165*clhs71 + clhs166);
lhs(9,4)=DN(2,0)*clhs45 + DN(2,1)*clhs139 + DN(2,2)*clhs140 + clhs238;
lhs(9,5)=DN(2,0)*clhs57 + DN(2,1)*clhs142 + DN(2,2)*clhs144 + clhs255 + clhs306;
lhs(9,6)=DN(2,0)*clhs64 + DN(2,1)*clhs147 + DN(2,2)*clhs149 + clhs269;
lhs(9,7)=-k*(DN(1,1)*clhs288 + DN(2,1)*clhs69 + clhs258*clhs71 + clhs259);
lhs(9,8)=DN(2,0)*clhs74 + DN(2,1)*clhs153 + DN(2,2)*clhs154 + clhs294;
lhs(9,9)=DN(2,0)*clhs86 + DN(2,1)*clhs156 + DN(2,2)*clhs158 + clhs11*clhs307 + clhs292;
lhs(9,10)=DN(2,0)*clhs93 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs309;
lhs(9,11)=-DN(2,1)*clhs296;
lhs(9,12)=DN(2,0)*clhs102 + DN(2,1)*clhs167 + DN(2,2)*clhs168 + clhs310;
lhs(9,13)=DN(2,0)*clhs114 + DN(2,1)*clhs170 + DN(2,2)*clhs172 + clhs311 + clhs312;
lhs(9,14)=DN(2,0)*clhs121 + DN(2,1)*clhs175 + DN(2,2)*clhs177 + clhs313;
lhs(9,15)=-k*(DN(2,1)*clhs126 + DN(3,1)*clhs288 + clhs314 + clhs315*clhs71);
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs129 + DN(2,2)*clhs181 + clhs96;
lhs(10,1)=DN(2,0)*clhs29 + DN(2,1)*clhs132 + DN(2,2)*clhs182 + clhs164;
lhs(10,2)=DN(2,0)*clhs36 + DN(2,1)*clhs136 + DN(2,2)*clhs183 + clhs199 + clhs305;
lhs(10,3)=-k*(DN(0,2)*clhs288 + DN(2,2)*clhs38 + clhs200*clhs71 + clhs201);
lhs(10,4)=DN(2,0)*clhs47 + DN(2,1)*clhs140 + DN(2,2)*clhs185 + clhs239;
lhs(10,5)=DN(2,0)*clhs60 + DN(2,1)*clhs144 + DN(2,2)*clhs188 + clhs257;
lhs(10,6)=DN(2,0)*clhs66 + DN(2,1)*clhs149 + DN(2,2)*clhs190 + clhs270 + clhs306;
lhs(10,7)=-k*(DN(1,2)*clhs288 + DN(2,2)*clhs69 + clhs271*clhs71 + clhs272);
lhs(10,8)=DN(2,0)*clhs76 + DN(2,1)*clhs154 + DN(2,2)*clhs194 + clhs295;
lhs(10,9)=DN(2,0)*clhs89 + DN(2,1)*clhs158 + DN(2,2)*clhs196 + clhs309;
lhs(10,10)=DN(2,0)*clhs95 + DN(2,1)*clhs163 + DN(2,2)*clhs198 + clhs11*clhs316 + clhs292;
lhs(10,11)=-DN(2,2)*clhs296;
lhs(10,12)=DN(2,0)*clhs104 + DN(2,1)*clhs168 + DN(2,2)*clhs202 + clhs318;
lhs(10,13)=DN(2,0)*clhs117 + DN(2,1)*clhs172 + DN(2,2)*clhs204 + clhs319;
lhs(10,14)=DN(2,0)*clhs123 + DN(2,1)*clhs177 + DN(2,2)*clhs206 + clhs312 + clhs320;
lhs(10,15)=-k*(DN(2,2)*clhs126 + DN(3,2)*clhs288 + clhs321 + clhs322*clhs71);
lhs(11,0)=DN(2,0)*clhs210 + clhs97;
lhs(11,1)=DN(2,1)*clhs210 + clhs165;
lhs(11,2)=DN(2,2)*clhs210 + clhs200;
lhs(11,3)=clhs219;
lhs(11,4)=DN(2,0)*clhs213 + clhs240;
lhs(11,5)=DN(2,1)*clhs213 + clhs258;
lhs(11,6)=DN(2,2)*clhs213 + clhs271;
lhs(11,7)=clhs282;
lhs(11,8)=DN(2,0)*clhs323;
lhs(11,9)=DN(2,1)*clhs323;
lhs(11,10)=DN(2,2)*clhs323;
lhs(11,11)=-clhs212*clhs290 - clhs212*clhs307 - clhs212*clhs316 - clhs291;
lhs(11,12)=DN(2,0)*clhs220 + clhs304;
lhs(11,13)=DN(2,1)*clhs220 + clhs315;
lhs(11,14)=DN(2,2)*clhs220 + clhs322;
lhs(11,15)=clhs324;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs107 + clhs328;
lhs(12,1)=DN(3,0)*clhs24 + DN(3,1)*clhs26 + DN(3,2)*clhs29 + clhs169;
lhs(12,2)=DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs203;
lhs(12,3)=-k*(DN(0,0)*clhs329 + DN(3,0)*clhs38 + clhs125*clhs71 + clhs127);
lhs(12,4)=DN(3,0)*clhs43 + DN(3,1)*clhs45 + DN(3,2)*clhs47 + clhs244 + clhs330;
lhs(12,5)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs260;
lhs(12,6)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs273;
lhs(12,7)=-k*(DN(1,0)*clhs329 + DN(3,0)*clhs69 + clhs248*clhs71 + clhs249);
lhs(12,8)=DN(3,0)*clhs72 + DN(3,1)*clhs74 + DN(3,2)*clhs76 + clhs299 + clhs331;
lhs(12,9)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs89 + clhs310;
lhs(12,10)=DN(3,0)*clhs91 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs318;
lhs(12,11)=-k*(DN(2,0)*clhs329 + DN(3,0)*clhs98 + clhs303*clhs71 + clhs304);
lhs(12,12)=DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs11*clhs332 + clhs334;
lhs(12,13)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs117 + clhs336;
lhs(12,14)=DN(3,0)*clhs119 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs337;
lhs(12,15)=-DN(3,0)*clhs338;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs128 + DN(3,2)*clhs129 + clhs118;
lhs(13,1)=DN(3,0)*clhs26 + DN(3,1)*clhs130 + DN(3,2)*clhs132 + clhs173 + clhs339;
lhs(13,2)=DN(3,0)*clhs34 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs205;
lhs(13,3)=-k*(DN(0,1)*clhs329 + DN(3,1)*clhs38 + clhs179*clhs71 + clhs180);
lhs(13,4)=DN(3,0)*clhs45 + DN(3,1)*clhs139 + DN(3,2)*clhs140 + clhs246;
lhs(13,5)=DN(3,0)*clhs57 + DN(3,1)*clhs142 + DN(3,2)*clhs144 + clhs261 + clhs340;
lhs(13,6)=DN(3,0)*clhs64 + DN(3,1)*clhs147 + DN(3,2)*clhs149 + clhs274;
lhs(13,7)=-k*(DN(1,1)*clhs329 + DN(3,1)*clhs69 + clhs264*clhs71 + clhs265);
lhs(13,8)=DN(3,0)*clhs74 + DN(3,1)*clhs153 + DN(3,2)*clhs154 + clhs301;
lhs(13,9)=DN(3,0)*clhs86 + DN(3,1)*clhs156 + DN(3,2)*clhs158 + clhs311 + clhs341;
lhs(13,10)=DN(3,0)*clhs93 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs319;
lhs(13,11)=-k*(DN(2,1)*clhs329 + DN(3,1)*clhs98 + clhs314*clhs71 + clhs315);
lhs(13,12)=DN(3,0)*clhs102 + DN(3,1)*clhs167 + DN(3,2)*clhs168 + clhs336;
lhs(13,13)=DN(3,0)*clhs114 + DN(3,1)*clhs170 + DN(3,2)*clhs172 + clhs11*clhs342 + clhs334;
lhs(13,14)=DN(3,0)*clhs121 + DN(3,1)*clhs175 + DN(3,2)*clhs177 + clhs343;
lhs(13,15)=-DN(3,1)*clhs338;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs129 + DN(3,2)*clhs181 + clhs124;
lhs(14,1)=DN(3,0)*clhs29 + DN(3,1)*clhs132 + DN(3,2)*clhs182 + clhs178;
lhs(14,2)=DN(3,0)*clhs36 + DN(3,1)*clhs136 + DN(3,2)*clhs183 + clhs207 + clhs339;
lhs(14,3)=-k*(DN(0,2)*clhs329 + DN(3,2)*clhs38 + clhs208*clhs71 + clhs209);
lhs(14,4)=DN(3,0)*clhs47 + DN(3,1)*clhs140 + DN(3,2)*clhs185 + clhs247;
lhs(14,5)=DN(3,0)*clhs60 + DN(3,1)*clhs144 + DN(3,2)*clhs188 + clhs263;
lhs(14,6)=DN(3,0)*clhs66 + DN(3,1)*clhs149 + DN(3,2)*clhs190 + clhs275 + clhs340;
lhs(14,7)=-k*(DN(1,2)*clhs329 + DN(3,2)*clhs69 + clhs276*clhs71 + clhs277);
lhs(14,8)=DN(3,0)*clhs76 + DN(3,1)*clhs154 + DN(3,2)*clhs194 + clhs302;
lhs(14,9)=DN(3,0)*clhs89 + DN(3,1)*clhs158 + DN(3,2)*clhs196 + clhs313;
lhs(14,10)=DN(3,0)*clhs95 + DN(3,1)*clhs163 + DN(3,2)*clhs198 + clhs320 + clhs341;
lhs(14,11)=-k*(DN(2,2)*clhs329 + DN(3,2)*clhs98 + clhs321*clhs71 + clhs322);
lhs(14,12)=DN(3,0)*clhs104 + DN(3,1)*clhs168 + DN(3,2)*clhs202 + clhs337;
lhs(14,13)=DN(3,0)*clhs117 + DN(3,1)*clhs172 + DN(3,2)*clhs204 + clhs343;
lhs(14,14)=DN(3,0)*clhs123 + DN(3,1)*clhs177 + DN(3,2)*clhs206 + clhs11*clhs344 + clhs334;
lhs(14,15)=-DN(3,2)*clhs338;
lhs(15,0)=DN(3,0)*clhs210 + clhs125;
lhs(15,1)=DN(3,1)*clhs210 + clhs179;
lhs(15,2)=DN(3,2)*clhs210 + clhs208;
lhs(15,3)=clhs221;
lhs(15,4)=DN(3,0)*clhs213 + clhs248;
lhs(15,5)=DN(3,1)*clhs213 + clhs264;
lhs(15,6)=DN(3,2)*clhs213 + clhs276;
lhs(15,7)=clhs283;
lhs(15,8)=DN(3,0)*clhs218 + clhs303;
lhs(15,9)=DN(3,1)*clhs218 + clhs314;
lhs(15,10)=DN(3,2)*clhs218 + clhs321;
lhs(15,11)=clhs324;
lhs(15,12)=DN(3,0)*clhs345;
lhs(15,13)=DN(3,1)*clhs345;
lhs(15,14)=DN(3,2)*clhs345;
lhs(15,15)=-clhs212*clhs332 - clhs212*clhs342 - clhs212*clhs344 - clhs333;


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
const double crhs8 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]);
const double crhs9 = N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]);
const double crhs10 = N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]);
const double crhs11 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs12 = crhs11 + crhs3;
const double crhs13 = k*(crhs7*h/stab_c1 + mu)*(-crhs10 + crhs12 - crhs8 - crhs9);
const double crhs14 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs15 = 1.0/(crhs7/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 = crhs1*rho - crhs2 - crhs6 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2]);
const double crhs17 = DN(0,0)*crhs4 + DN(0,1)*crhs5;
const double crhs18 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1);
const double crhs19 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs20 = rho*(crhs11*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs21 = crhs18*rho - crhs19 - crhs20 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2]);
const double crhs22 = crhs10 + crhs8 + crhs9;
const double crhs23 = 1.0*crhs15;
const double crhs24 = crhs16*crhs23;
const double crhs25 = crhs21*crhs23;
const double crhs26 = DN(1,0)*crhs4 + DN(1,1)*crhs5;
const double crhs27 = DN(2,0)*crhs4 + DN(2,1)*crhs5;
rhs[0]=DN(0,0)*crhs0*k - DN(0,0)*crhs13 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1*rho + 1.0*N[0]*crhs14*crhs15*crhs16*rho - N[0]*crhs2 - N[0]*crhs6 + 1.0*crhs15*crhs16*crhs17*rho;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0*k - DN(0,1)*crhs13 - DN(0,1)*stress[1] + 1.0*N[0]*crhs14*crhs15*crhs21*rho + N[0]*crhs18*rho - N[0]*crhs19 - N[0]*crhs20 + 1.0*crhs15*crhs17*crhs21*rho;
rhs[2]=DN(0,0)*crhs24 + DN(0,1)*crhs25 - N[0]*crhs12 + N[0]*crhs22;
rhs[3]=DN(1,0)*crhs0*k - DN(1,0)*crhs13 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1*rho + 1.0*N[1]*crhs14*crhs15*crhs16*rho - N[1]*crhs2 - N[1]*crhs6 + 1.0*crhs15*crhs16*crhs26*rho;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0*k - DN(1,1)*crhs13 - DN(1,1)*stress[1] + 1.0*N[1]*crhs14*crhs15*crhs21*rho + N[1]*crhs18*rho - N[1]*crhs19 - N[1]*crhs20 + 1.0*crhs15*crhs21*crhs26*rho;
rhs[5]=DN(1,0)*crhs24 + DN(1,1)*crhs25 - N[1]*crhs12 + N[1]*crhs22;
rhs[6]=DN(2,0)*crhs0*k - DN(2,0)*crhs13 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1*rho + 1.0*N[2]*crhs14*crhs15*crhs16*rho - N[2]*crhs2 - N[2]*crhs6 + 1.0*crhs15*crhs16*crhs27*rho;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0*k - DN(2,1)*crhs13 - DN(2,1)*stress[1] + 1.0*N[2]*crhs14*crhs15*crhs21*rho + N[2]*crhs18*rho - N[2]*crhs19 - N[2]*crhs20 + 1.0*crhs15*crhs21*crhs27*rho;
rhs[8]=DN(2,0)*crhs24 + DN(2,1)*crhs25 - N[2]*crhs12 + N[2]*crhs22;


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
const double crhs9 = N[0]*(bdf0*eps_vol[0] + bdf1*eps_vol_n[0] + bdf2*eps_vol_nn[0]);
const double crhs10 = N[1]*(bdf0*eps_vol[1] + bdf1*eps_vol_n[1] + bdf2*eps_vol_nn[1]);
const double crhs11 = N[2]*(bdf0*eps_vol[2] + bdf1*eps_vol_n[2] + bdf2*eps_vol_nn[2]);
const double crhs12 = N[3]*(bdf0*eps_vol[3] + bdf1*eps_vol_n[3] + bdf2*eps_vol_nn[3]);
const double crhs13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs14 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs15 = crhs13 + crhs14 + crhs3;
const double crhs16 = k*(crhs8*h/stab_c1 + mu)*(-crhs10 - crhs11 - crhs12 + crhs15 - crhs9);
const double crhs17 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs18 = 1.0/(crhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs19 = crhs1*rho - crhs2 - crhs7 + k*(DN(0,0)*eps_vol[0] + DN(1,0)*eps_vol[1] + DN(2,0)*eps_vol[2] + DN(3,0)*eps_vol[3]);
const double crhs20 = DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6;
const double crhs21 = N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1);
const double crhs22 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs23 = rho*(crhs13*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs24 = crhs21*rho - crhs22 - crhs23 + k*(DN(0,1)*eps_vol[0] + DN(1,1)*eps_vol[1] + DN(2,1)*eps_vol[2] + DN(3,1)*eps_vol[3]);
const double crhs25 = N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2);
const double crhs26 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs27 = rho*(crhs14*crhs6 + crhs4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs28 = crhs25*rho - crhs26 - crhs27 + k*(DN(0,2)*eps_vol[0] + DN(1,2)*eps_vol[1] + DN(2,2)*eps_vol[2] + DN(3,2)*eps_vol[3]);
const double crhs29 = crhs10 + crhs11 + crhs12 + crhs9;
const double crhs30 = 1.0*crhs18;
const double crhs31 = crhs19*crhs30;
const double crhs32 = crhs24*crhs30;
const double crhs33 = crhs28*crhs30;
const double crhs34 = DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6;
const double crhs35 = DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6;
const double crhs36 = DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6;
rhs[0]=DN(0,0)*crhs0*k - DN(0,0)*crhs16 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1*rho + 1.0*N[0]*crhs17*crhs18*crhs19*rho - N[0]*crhs2 - N[0]*crhs7 + 1.0*crhs18*crhs19*crhs20*rho;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0*k - DN(0,1)*crhs16 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + 1.0*N[0]*crhs17*crhs18*crhs24*rho + N[0]*crhs21*rho - N[0]*crhs22 - N[0]*crhs23 + 1.0*crhs18*crhs20*crhs24*rho;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0*k - DN(0,2)*crhs16 - DN(0,2)*stress[2] + 1.0*N[0]*crhs17*crhs18*crhs28*rho + N[0]*crhs25*rho - N[0]*crhs26 - N[0]*crhs27 + 1.0*crhs18*crhs20*crhs28*rho;
rhs[3]=DN(0,0)*crhs31 + DN(0,1)*crhs32 + DN(0,2)*crhs33 - N[0]*crhs15 + N[0]*crhs29;
rhs[4]=DN(1,0)*crhs0*k - DN(1,0)*crhs16 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1*rho + 1.0*N[1]*crhs17*crhs18*crhs19*rho - N[1]*crhs2 - N[1]*crhs7 + 1.0*crhs18*crhs19*crhs34*rho;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0*k - DN(1,1)*crhs16 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + 1.0*N[1]*crhs17*crhs18*crhs24*rho + N[1]*crhs21*rho - N[1]*crhs22 - N[1]*crhs23 + 1.0*crhs18*crhs24*crhs34*rho;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0*k - DN(1,2)*crhs16 - DN(1,2)*stress[2] + 1.0*N[1]*crhs17*crhs18*crhs28*rho + N[1]*crhs25*rho - N[1]*crhs26 - N[1]*crhs27 + 1.0*crhs18*crhs28*crhs34*rho;
rhs[7]=DN(1,0)*crhs31 + DN(1,1)*crhs32 + DN(1,2)*crhs33 - N[1]*crhs15 + N[1]*crhs29;
rhs[8]=DN(2,0)*crhs0*k - DN(2,0)*crhs16 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1*rho + 1.0*N[2]*crhs17*crhs18*crhs19*rho - N[2]*crhs2 - N[2]*crhs7 + 1.0*crhs18*crhs19*crhs35*rho;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0*k - DN(2,1)*crhs16 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + 1.0*N[2]*crhs17*crhs18*crhs24*rho + N[2]*crhs21*rho - N[2]*crhs22 - N[2]*crhs23 + 1.0*crhs18*crhs24*crhs35*rho;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0*k - DN(2,2)*crhs16 - DN(2,2)*stress[2] + 1.0*N[2]*crhs17*crhs18*crhs28*rho + N[2]*crhs25*rho - N[2]*crhs26 - N[2]*crhs27 + 1.0*crhs18*crhs28*crhs35*rho;
rhs[11]=DN(2,0)*crhs31 + DN(2,1)*crhs32 + DN(2,2)*crhs33 - N[2]*crhs15 + N[2]*crhs29;
rhs[12]=DN(3,0)*crhs0*k - DN(3,0)*crhs16 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1*rho + 1.0*N[3]*crhs17*crhs18*crhs19*rho - N[3]*crhs2 - N[3]*crhs7 + 1.0*crhs18*crhs19*crhs36*rho;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0*k - DN(3,1)*crhs16 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + 1.0*N[3]*crhs17*crhs18*crhs24*rho + N[3]*crhs21*rho - N[3]*crhs22 - N[3]*crhs23 + 1.0*crhs18*crhs24*crhs36*rho;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0*k - DN(3,2)*crhs16 - DN(3,2)*stress[2] + 1.0*N[3]*crhs17*crhs18*crhs28*rho + N[3]*crhs25*rho - N[3]*crhs26 - N[3]*crhs27 + 1.0*crhs18*crhs28*crhs36*rho;
rhs[15]=DN(3,0)*crhs31 + DN(3,1)*crhs32 + DN(3,2)*crhs33 - N[3]*crhs15 + N[3]*crhs29;


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