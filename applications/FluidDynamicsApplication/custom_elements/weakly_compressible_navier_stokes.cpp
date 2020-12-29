#include "weakly_compressible_navier_stokes.h"
#include "custom_utilities/weakly_compressible_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::WeaklyCompressibleNavierStokes(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
WeaklyCompressibleNavierStokes<TElementData>::~WeaklyCompressibleNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokes<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokes>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer WeaklyCompressibleNavierStokes<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<WeaklyCompressibleNavierStokes>(NewId, pGeom, pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int WeaklyCompressibleNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
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

template <class TElementData>
std::string WeaklyCompressibleNavierStokes<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "WeaklyCompressibleNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::PrintInfo(std::ostream& rOStream) const
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
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData,
    MatrixType& rLHS,
    VectorType& rRHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData,
    MatrixType& rLHS)
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData,
    VectorType& rRHS)
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void WeaklyCompressibleNavierStokes<TElementData>::AddBoundaryTraction(
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
void WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> >::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesData<2,3>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const array_1d<double,3>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

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

    // LHS Gauss point contribution
    const double weight = rData.Weight;

    const double crLHS0 =     pow(weight, 2);
const double crLHS1 =     C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS2 =     C(0,2)*DN(0,0);
const double crLHS3 =     C(2,2)*DN(0,1) + crLHS2;
const double crLHS4 =     pow(DN(0,0), 2);
const double crLHS5 =     N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crLHS6 =     N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS7 =     N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS8 =     crLHS5*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2));
const double crLHS9 =     crLHS8*h/stab_c1 + mu;
const double crLHS10 =     pow(N[0], 2);
const double crLHS11 =     bdf0*crLHS5;
const double crLHS12 =     N[0]*crLHS5;
const double crLHS13 =     DN(0,0)*crLHS6 + DN(0,1)*crLHS7;
const double crLHS14 =     N[0]*bdf0 + crLHS13;
const double crLHS15 =     pow(crLHS5, 2);
const double crLHS16 =     DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS17 =     1.0/(crLHS5*dyn_tau/dt + crLHS8/h + mu*stab_c1/pow(h, 2));
const double crLHS18 =     1.0*N[0]*crLHS15*crLHS16*crLHS17;
const double crLHS19 =     1.0*crLHS13*crLHS15*crLHS17;
const double crLHS20 =     crLHS10*crLHS11 + crLHS12*crLHS13 + crLHS14*crLHS18 + crLHS14*crLHS19;
const double crLHS21 =     C(0,1)*DN(0,1) + crLHS2;
const double crLHS22 =     C(1,2)*DN(0,1);
const double crLHS23 =     C(2,2)*DN(0,0) + crLHS22;
const double crLHS24 =     DN(0,0)*crLHS9;
const double crLHS25 =     DN(0,1)*crLHS24;
const double crLHS26 =     pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], -2);
const double crLHS27 =     1.0/crLHS5;
const double crLHS28 =     N[0]*bdf0*crLHS26*crLHS27;
const double crLHS29 =     1.0*crLHS16*crLHS17;
const double crLHS30 =     1.0*crLHS17*crLHS5;
const double crLHS31 =     crLHS13*crLHS30;
const double crLHS32 =     crLHS0*(-N[0] + crLHS12*crLHS29 + crLHS28*crLHS9 + crLHS31);
const double crLHS33 =     C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS34 =     C(0,2)*DN(1,0);
const double crLHS35 =     C(2,2)*DN(1,1) + crLHS34;
const double crLHS36 =     DN(1,0)*crLHS24;
const double crLHS37 =     N[0]*bdf0*crLHS5;
const double crLHS38 =     N[1]*crLHS37;
const double crLHS39 =     DN(1,0)*crLHS6 + DN(1,1)*crLHS7;
const double crLHS40 =     N[1]*bdf0;
const double crLHS41 =     crLHS39 + crLHS40;
const double crLHS42 =     crLHS12*crLHS39 + crLHS18*crLHS41 + crLHS19*crLHS41 + crLHS38;
const double crLHS43 =     C(0,1)*DN(1,1) + crLHS34;
const double crLHS44 =     C(1,2)*DN(1,1);
const double crLHS45 =     C(2,2)*DN(1,0) + crLHS44;
const double crLHS46 =     DN(1,1)*crLHS24;
const double crLHS47 =     DN(0,0)*N[1];
const double crLHS48 =     bdf0*crLHS26*crLHS27*crLHS9;
const double crLHS49 =     DN(1,0)*N[0];
const double crLHS50 =     1.0*crLHS16*crLHS17*crLHS5;
const double crLHS51 =     1.0*DN(1,0)*crLHS17*crLHS5;
const double crLHS52 =     C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS53 =     C(0,2)*DN(2,0);
const double crLHS54 =     C(2,2)*DN(2,1) + crLHS53;
const double crLHS55 =     DN(2,0)*crLHS24;
const double crLHS56 =     N[2]*crLHS37;
const double crLHS57 =     DN(2,0)*crLHS6 + DN(2,1)*crLHS7;
const double crLHS58 =     N[2]*bdf0;
const double crLHS59 =     crLHS57 + crLHS58;
const double crLHS60 =     crLHS12*crLHS57 + crLHS18*crLHS59 + crLHS19*crLHS59 + crLHS56;
const double crLHS61 =     C(0,1)*DN(2,1) + crLHS53;
const double crLHS62 =     C(1,2)*DN(2,1);
const double crLHS63 =     C(2,2)*DN(2,0) + crLHS62;
const double crLHS64 =     DN(2,1)*crLHS24;
const double crLHS65 =     DN(0,0)*N[2];
const double crLHS66 =     DN(2,0)*N[0];
const double crLHS67 =     C(0,1)*DN(0,0) + crLHS22;
const double crLHS68 =     C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS69 =     pow(DN(0,1), 2);
const double crLHS70 =     C(0,1)*DN(1,0) + crLHS44;
const double crLHS71 =     DN(0,1)*crLHS9;
const double crLHS72 =     DN(1,0)*crLHS71;
const double crLHS73 =     C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS74 =     DN(1,1)*crLHS71;
const double crLHS75 =     DN(0,1)*N[1];
const double crLHS76 =     DN(1,1)*N[0];
const double crLHS77 =     1.0*DN(1,1)*crLHS17*crLHS5;
const double crLHS78 =     C(0,1)*DN(2,0) + crLHS62;
const double crLHS79 =     DN(2,0)*crLHS71;
const double crLHS80 =     C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS81 =     DN(2,1)*crLHS71;
const double crLHS82 =     DN(0,1)*N[2];
const double crLHS83 =     DN(2,1)*N[0];
const double crLHS84 =     crLHS14*crLHS30;
const double crLHS85 =     crLHS0*(N[0] + crLHS84);
const double crLHS86 =     bdf0*crLHS26*crLHS27;
const double crLHS87 =     1.0*crLHS17;
const double crLHS88 =     1.0*DN(0,0)*crLHS17*crLHS5;
const double crLHS89 =     1.0*DN(0,1)*crLHS17*crLHS5;
const double crLHS90 =     1.0*DN(0,0)*crLHS17;
const double crLHS91 =     1.0*DN(0,1)*crLHS17;
const double crLHS92 =     crLHS0*(DN(1,0)*crLHS90 + DN(1,1)*crLHS91 + N[1]*crLHS28);
const double crLHS93 =     crLHS0*(DN(2,0)*crLHS90 + DN(2,1)*crLHS91 + N[2]*crLHS28);
const double crLHS94 =     N[1]*crLHS5;
const double crLHS95 =     1.0*N[1]*crLHS15*crLHS16*crLHS17;
const double crLHS96 =     1.0*crLHS15*crLHS17*crLHS39;
const double crLHS97 =     crLHS13*crLHS94 + crLHS14*crLHS95 + crLHS14*crLHS96 + crLHS38;
const double crLHS98 =     pow(DN(1,0), 2);
const double crLHS99 =     pow(N[1], 2);
const double crLHS100 =     crLHS11*crLHS99 + crLHS39*crLHS94 + crLHS41*crLHS95 + crLHS41*crLHS96;
const double crLHS101 =     DN(1,0)*crLHS9;
const double crLHS102 =     DN(1,1)*crLHS101;
const double crLHS103 =     crLHS26*crLHS27*crLHS9;
const double crLHS104 =     crLHS30*crLHS39;
const double crLHS105 =     crLHS0*(-N[1] + crLHS103*crLHS40 + crLHS104 + crLHS29*crLHS94);
const double crLHS106 =     DN(2,0)*crLHS101;
const double crLHS107 =     N[1]*N[2]*bdf0;
const double crLHS108 =     crLHS107*crLHS5;
const double crLHS109 =     crLHS108 + crLHS57*crLHS94 + crLHS59*crLHS95 + crLHS59*crLHS96;
const double crLHS110 =     DN(2,1)*crLHS101;
const double crLHS111 =     DN(1,0)*N[2];
const double crLHS112 =     DN(2,0)*N[1];
const double crLHS113 =     pow(DN(1,1), 2);
const double crLHS114 =     DN(1,1)*crLHS9;
const double crLHS115 =     DN(2,0)*crLHS114;
const double crLHS116 =     DN(2,1)*crLHS114;
const double crLHS117 =     DN(1,1)*N[2];
const double crLHS118 =     DN(2,1)*N[1];
const double crLHS119 =     crLHS30*crLHS41;
const double crLHS120 =     crLHS0*(N[1] + crLHS119);
const double crLHS121 =     crLHS0*(1.0*DN(1,0)*DN(2,0)*crLHS17 + 1.0*DN(1,1)*DN(2,1)*crLHS17 + crLHS107*crLHS26*crLHS27);
const double crLHS122 =     N[2]*crLHS5;
const double crLHS123 =     1.0*N[2]*crLHS15*crLHS16*crLHS17;
const double crLHS124 =     1.0*crLHS15*crLHS17*crLHS57;
const double crLHS125 =     crLHS122*crLHS13 + crLHS123*crLHS14 + crLHS124*crLHS14 + crLHS56;
const double crLHS126 =     crLHS108 + crLHS122*crLHS39 + crLHS123*crLHS41 + crLHS124*crLHS41;
const double crLHS127 =     pow(DN(2,0), 2);
const double crLHS128 =     pow(N[2], 2);
const double crLHS129 =     crLHS11*crLHS128 + crLHS122*crLHS57 + crLHS123*crLHS59 + crLHS124*crLHS59;
const double crLHS130 =     DN(2,0)*DN(2,1)*crLHS9;
const double crLHS131 =     crLHS0*(-N[2] + crLHS103*crLHS58 + crLHS122*crLHS29 + crLHS30*crLHS57);
const double crLHS132 =     pow(DN(2,1), 2);
const double crLHS133 =     crLHS0*(N[2] + crLHS30*crLHS59);
    rLHS(0,0)+=crLHS0*(DN(0,0)*crLHS1 + DN(0,1)*crLHS3 + crLHS20 + crLHS4*crLHS9);
    rLHS(0,1)+=crLHS0*(DN(0,0)*crLHS21 + DN(0,1)*crLHS23 + crLHS25);
    rLHS(0,2)+=DN(0,0)*crLHS32;
    rLHS(0,3)+=crLHS0*(DN(0,0)*crLHS33 + DN(0,1)*crLHS35 + crLHS36 + crLHS42);
    rLHS(0,4)+=crLHS0*(DN(0,0)*crLHS43 + DN(0,1)*crLHS45 + crLHS46);
    rLHS(0,5)+=crLHS0*(crLHS13*crLHS51 + crLHS47*crLHS48 - crLHS47 + crLHS49*crLHS50);
    rLHS(0,6)+=crLHS0*(DN(0,0)*crLHS52 + DN(0,1)*crLHS54 + crLHS55 + crLHS60);
    rLHS(0,7)+=crLHS0*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + crLHS64);
    rLHS(0,8)+=crLHS0*(DN(2,0)*crLHS31 + crLHS48*crLHS65 + crLHS50*crLHS66 - crLHS65);
    rLHS(1,0)+=crLHS0*(DN(0,0)*crLHS3 + DN(0,1)*crLHS67 + crLHS25);
    rLHS(1,1)+=crLHS0*(DN(0,0)*crLHS23 + DN(0,1)*crLHS68 + crLHS20 + crLHS69*crLHS9);
    rLHS(1,2)+=DN(0,1)*crLHS32;
    rLHS(1,3)+=crLHS0*(DN(0,0)*crLHS35 + DN(0,1)*crLHS70 + crLHS72);
    rLHS(1,4)+=crLHS0*(DN(0,0)*crLHS45 + DN(0,1)*crLHS73 + crLHS42 + crLHS74);
    rLHS(1,5)+=crLHS0*(crLHS13*crLHS77 + crLHS48*crLHS75 + crLHS50*crLHS76 - crLHS75);
    rLHS(1,6)+=crLHS0*(DN(0,0)*crLHS54 + DN(0,1)*crLHS78 + crLHS79);
    rLHS(1,7)+=crLHS0*(DN(0,0)*crLHS63 + DN(0,1)*crLHS80 + crLHS60 + crLHS81);
    rLHS(1,8)+=crLHS0*(DN(2,1)*crLHS31 + crLHS48*crLHS82 + crLHS50*crLHS83 - crLHS82);
    rLHS(2,0)+=DN(0,0)*crLHS85;
    rLHS(2,1)+=DN(0,1)*crLHS85;
    rLHS(2,2)+=crLHS0*(crLHS10*crLHS86 + crLHS4*crLHS87 + crLHS69*crLHS87);
    rLHS(2,3)+=crLHS0*(crLHS41*crLHS88 + crLHS49);
    rLHS(2,4)+=crLHS0*(crLHS41*crLHS89 + crLHS76);
    rLHS(2,5)+=crLHS92;
    rLHS(2,6)+=crLHS0*(crLHS59*crLHS88 + crLHS66);
    rLHS(2,7)+=crLHS0*(crLHS59*crLHS89 + crLHS83);
    rLHS(2,8)+=crLHS93;
    rLHS(3,0)+=crLHS0*(DN(1,0)*crLHS1 + DN(1,1)*crLHS3 + crLHS36 + crLHS97);
    rLHS(3,1)+=crLHS0*(DN(1,0)*crLHS21 + DN(1,1)*crLHS23 + crLHS72);
    rLHS(3,2)+=crLHS0*(crLHS39*crLHS88 + crLHS47*crLHS50 + crLHS48*crLHS49 - crLHS49);
    rLHS(3,3)+=crLHS0*(DN(1,0)*crLHS33 + DN(1,1)*crLHS35 + crLHS100 + crLHS9*crLHS98);
    rLHS(3,4)+=crLHS0*(DN(1,0)*crLHS43 + DN(1,1)*crLHS45 + crLHS102);
    rLHS(3,5)+=DN(1,0)*crLHS105;
    rLHS(3,6)+=crLHS0*(DN(1,0)*crLHS52 + DN(1,1)*crLHS54 + crLHS106 + crLHS109);
    rLHS(3,7)+=crLHS0*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + crLHS110);
    rLHS(3,8)+=crLHS0*(DN(2,0)*crLHS104 + crLHS111*crLHS48 - crLHS111 + crLHS112*crLHS50);
    rLHS(4,0)+=crLHS0*(DN(1,0)*crLHS3 + DN(1,1)*crLHS67 + crLHS46);
    rLHS(4,1)+=crLHS0*(DN(1,0)*crLHS23 + DN(1,1)*crLHS68 + crLHS74 + crLHS97);
    rLHS(4,2)+=crLHS0*(crLHS39*crLHS89 + crLHS48*crLHS76 + crLHS50*crLHS75 - crLHS76);
    rLHS(4,3)+=crLHS0*(DN(1,0)*crLHS35 + DN(1,1)*crLHS70 + crLHS102);
    rLHS(4,4)+=crLHS0*(DN(1,0)*crLHS45 + DN(1,1)*crLHS73 + crLHS100 + crLHS113*crLHS9);
    rLHS(4,5)+=DN(1,1)*crLHS105;
    rLHS(4,6)+=crLHS0*(DN(1,0)*crLHS54 + DN(1,1)*crLHS78 + crLHS115);
    rLHS(4,7)+=crLHS0*(DN(1,0)*crLHS63 + DN(1,1)*crLHS80 + crLHS109 + crLHS116);
    rLHS(4,8)+=crLHS0*(DN(2,1)*crLHS104 + crLHS117*crLHS48 - crLHS117 + crLHS118*crLHS50);
    rLHS(5,0)+=crLHS0*(crLHS14*crLHS51 + crLHS47);
    rLHS(5,1)+=crLHS0*(crLHS14*crLHS77 + crLHS75);
    rLHS(5,2)+=crLHS92;
    rLHS(5,3)+=DN(1,0)*crLHS120;
    rLHS(5,4)+=DN(1,1)*crLHS120;
    rLHS(5,5)+=crLHS0*(crLHS113*crLHS87 + crLHS86*crLHS99 + crLHS87*crLHS98);
    rLHS(5,6)+=crLHS0*(crLHS112 + crLHS51*crLHS59);
    rLHS(5,7)+=crLHS0*(crLHS118 + crLHS59*crLHS77);
    rLHS(5,8)+=crLHS121;
    rLHS(6,0)+=crLHS0*(DN(2,0)*crLHS1 + DN(2,1)*crLHS3 + crLHS125 + crLHS55);
    rLHS(6,1)+=crLHS0*(DN(2,0)*crLHS21 + DN(2,1)*crLHS23 + crLHS79);
    rLHS(6,2)+=crLHS0*(crLHS48*crLHS66 + crLHS50*crLHS65 + crLHS57*crLHS88 - crLHS66);
    rLHS(6,3)+=crLHS0*(DN(2,0)*crLHS33 + DN(2,1)*crLHS35 + crLHS106 + crLHS126);
    rLHS(6,4)+=crLHS0*(DN(2,0)*crLHS43 + DN(2,1)*crLHS45 + crLHS115);
    rLHS(6,5)+=crLHS0*(crLHS111*crLHS50 + crLHS112*crLHS48 - crLHS112 + crLHS51*crLHS57);
    rLHS(6,6)+=crLHS0*(DN(2,0)*crLHS52 + DN(2,1)*crLHS54 + crLHS127*crLHS9 + crLHS129);
    rLHS(6,7)+=crLHS0*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + crLHS130);
    rLHS(6,8)+=DN(2,0)*crLHS131;
    rLHS(7,0)+=crLHS0*(DN(2,0)*crLHS3 + DN(2,1)*crLHS67 + crLHS64);
    rLHS(7,1)+=crLHS0*(DN(2,0)*crLHS23 + DN(2,1)*crLHS68 + crLHS125 + crLHS81);
    rLHS(7,2)+=crLHS0*(crLHS48*crLHS83 + crLHS50*crLHS82 + crLHS57*crLHS89 - crLHS83);
    rLHS(7,3)+=crLHS0*(DN(2,0)*crLHS35 + DN(2,1)*crLHS70 + crLHS110);
    rLHS(7,4)+=crLHS0*(DN(2,0)*crLHS45 + DN(2,1)*crLHS73 + crLHS116 + crLHS126);
    rLHS(7,5)+=crLHS0*(crLHS117*crLHS50 + crLHS118*crLHS48 - crLHS118 + crLHS57*crLHS77);
    rLHS(7,6)+=crLHS0*(DN(2,0)*crLHS54 + DN(2,1)*crLHS78 + crLHS130);
    rLHS(7,7)+=crLHS0*(DN(2,0)*crLHS63 + DN(2,1)*crLHS80 + crLHS129 + crLHS132*crLHS9);
    rLHS(7,8)+=DN(2,1)*crLHS131;
    rLHS(8,0)+=crLHS0*(DN(2,0)*crLHS84 + crLHS65);
    rLHS(8,1)+=crLHS0*(DN(2,1)*crLHS84 + crLHS82);
    rLHS(8,2)+=crLHS93;
    rLHS(8,3)+=crLHS0*(DN(2,0)*crLHS119 + crLHS111);
    rLHS(8,4)+=crLHS0*(DN(2,1)*crLHS119 + crLHS117);
    rLHS(8,5)+=crLHS121;
    rLHS(8,6)+=DN(2,0)*crLHS133;
    rLHS(8,7)+=DN(2,1)*crLHS133;
    rLHS(8,8)+=crLHS0*(crLHS127*crLHS87 + crLHS128*crLHS86 + crLHS132*crLHS87);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointLHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    MatrixType& rLHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const array_1d<double,4>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

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

    // LHS Gauss point contribution
    const double weight = rData.Weight;

    const double crLHS0 =     pow(weight, 2);
const double crLHS1 =     C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS2 =     C(0,3)*DN(0,0);
const double crLHS3 =     C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS2;
const double crLHS4 =     C(0,5)*DN(0,0);
const double crLHS5 =     C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS4;
const double crLHS6 =     pow(DN(0,0), 2);
const double crLHS7 =     N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crLHS8 =     N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 =     N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS10 =     N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS11 =     crLHS7*stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2));
const double crLHS12 =     crLHS11*h/stab_c1 + mu;
const double crLHS13 =     pow(N[0], 2);
const double crLHS14 =     bdf0*crLHS7;
const double crLHS15 =     N[0]*crLHS7;
const double crLHS16 =     DN(0,0)*crLHS8 + DN(0,1)*crLHS9 + DN(0,2)*crLHS10;
const double crLHS17 =     N[0]*bdf0 + crLHS16;
const double crLHS18 =     pow(crLHS7, 2);
const double crLHS19 =     DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS20 =     1.0/(crLHS11/h + crLHS7*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double crLHS21 =     1.0*N[0]*crLHS18*crLHS19*crLHS20;
const double crLHS22 =     1.0*crLHS16*crLHS18*crLHS20;
const double crLHS23 =     crLHS13*crLHS14 + crLHS15*crLHS16 + crLHS17*crLHS21 + crLHS17*crLHS22;
const double crLHS24 =     C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS2;
const double crLHS25 =     C(1,3)*DN(0,1);
const double crLHS26 =     C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS25;
const double crLHS27 =     C(3,5)*DN(0,0);
const double crLHS28 =     C(4,5)*DN(0,2);
const double crLHS29 =     C(1,5)*DN(0,1) + crLHS27 + crLHS28;
const double crLHS30 =     DN(0,0)*crLHS12;
const double crLHS31 =     DN(0,1)*crLHS30;
const double crLHS32 =     C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS4;
const double crLHS33 =     C(3,4)*DN(0,1);
const double crLHS34 =     C(2,3)*DN(0,2) + crLHS27 + crLHS33;
const double crLHS35 =     C(2,5)*DN(0,2);
const double crLHS36 =     C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS35;
const double crLHS37 =     DN(0,2)*crLHS30;
const double crLHS38 =     pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], -2);
const double crLHS39 =     1.0/crLHS7;
const double crLHS40 =     N[0]*bdf0*crLHS38*crLHS39;
const double crLHS41 =     1.0*crLHS19*crLHS20;
const double crLHS42 =     1.0*crLHS20*crLHS7;
const double crLHS43 =     crLHS16*crLHS42;
const double crLHS44 =     crLHS0*(-N[0] + crLHS12*crLHS40 + crLHS15*crLHS41 + crLHS43);
const double crLHS45 =     C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS46 =     C(0,3)*DN(1,0);
const double crLHS47 =     C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS46;
const double crLHS48 =     C(0,5)*DN(1,0);
const double crLHS49 =     C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS48;
const double crLHS50 =     DN(1,0)*crLHS30;
const double crLHS51 =     N[0]*bdf0*crLHS7;
const double crLHS52 =     N[1]*crLHS51;
const double crLHS53 =     DN(1,0)*crLHS8 + DN(1,1)*crLHS9 + DN(1,2)*crLHS10;
const double crLHS54 =     N[1]*bdf0 + crLHS53;
const double crLHS55 =     crLHS15*crLHS53 + crLHS21*crLHS54 + crLHS22*crLHS54 + crLHS52;
const double crLHS56 =     C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS46;
const double crLHS57 =     C(1,3)*DN(1,1);
const double crLHS58 =     C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS57;
const double crLHS59 =     C(3,5)*DN(1,0);
const double crLHS60 =     C(4,5)*DN(1,2);
const double crLHS61 =     C(1,5)*DN(1,1) + crLHS59 + crLHS60;
const double crLHS62 =     DN(1,1)*crLHS30;
const double crLHS63 =     C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS48;
const double crLHS64 =     C(3,4)*DN(1,1);
const double crLHS65 =     C(2,3)*DN(1,2) + crLHS59 + crLHS64;
const double crLHS66 =     C(2,5)*DN(1,2);
const double crLHS67 =     C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS66;
const double crLHS68 =     DN(1,2)*crLHS30;
const double crLHS69 =     DN(0,0)*N[1];
const double crLHS70 =     bdf0*crLHS12*crLHS38*crLHS39;
const double crLHS71 =     DN(1,0)*N[0];
const double crLHS72 =     1.0*crLHS19*crLHS20*crLHS7;
const double crLHS73 =     1.0*DN(1,0)*crLHS20*crLHS7;
const double crLHS74 =     C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS75 =     C(0,3)*DN(2,0);
const double crLHS76 =     C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS75;
const double crLHS77 =     C(0,5)*DN(2,0);
const double crLHS78 =     C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS77;
const double crLHS79 =     DN(2,0)*crLHS30;
const double crLHS80 =     N[2]*crLHS51;
const double crLHS81 =     DN(2,0)*crLHS8 + DN(2,1)*crLHS9 + DN(2,2)*crLHS10;
const double crLHS82 =     N[2]*bdf0;
const double crLHS83 =     crLHS81 + crLHS82;
const double crLHS84 =     crLHS15*crLHS81 + crLHS21*crLHS83 + crLHS22*crLHS83 + crLHS80;
const double crLHS85 =     C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS75;
const double crLHS86 =     C(1,3)*DN(2,1);
const double crLHS87 =     C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS86;
const double crLHS88 =     C(3,5)*DN(2,0);
const double crLHS89 =     C(4,5)*DN(2,2);
const double crLHS90 =     C(1,5)*DN(2,1) + crLHS88 + crLHS89;
const double crLHS91 =     DN(2,1)*crLHS30;
const double crLHS92 =     C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS77;
const double crLHS93 =     C(3,4)*DN(2,1);
const double crLHS94 =     C(2,3)*DN(2,2) + crLHS88 + crLHS93;
const double crLHS95 =     C(2,5)*DN(2,2);
const double crLHS96 =     C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS95;
const double crLHS97 =     DN(2,2)*crLHS30;
const double crLHS98 =     DN(0,0)*N[2];
const double crLHS99 =     DN(2,0)*N[0];
const double crLHS100 =     1.0*DN(2,0)*crLHS20*crLHS7;
const double crLHS101 =     C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS102 =     C(0,3)*DN(3,0);
const double crLHS103 =     C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS102;
const double crLHS104 =     C(0,5)*DN(3,0);
const double crLHS105 =     C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS104;
const double crLHS106 =     DN(3,0)*crLHS30;
const double crLHS107 =     N[3]*crLHS51;
const double crLHS108 =     DN(3,0)*crLHS8 + DN(3,1)*crLHS9 + DN(3,2)*crLHS10;
const double crLHS109 =     N[3]*bdf0;
const double crLHS110 =     crLHS108 + crLHS109;
const double crLHS111 =     crLHS107 + crLHS108*crLHS15 + crLHS110*crLHS21 + crLHS110*crLHS22;
const double crLHS112 =     C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS102;
const double crLHS113 =     C(1,3)*DN(3,1);
const double crLHS114 =     C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS113;
const double crLHS115 =     C(3,5)*DN(3,0);
const double crLHS116 =     C(4,5)*DN(3,2);
const double crLHS117 =     C(1,5)*DN(3,1) + crLHS115 + crLHS116;
const double crLHS118 =     DN(3,1)*crLHS30;
const double crLHS119 =     C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS104;
const double crLHS120 =     C(3,4)*DN(3,1);
const double crLHS121 =     C(2,3)*DN(3,2) + crLHS115 + crLHS120;
const double crLHS122 =     C(2,5)*DN(3,2);
const double crLHS123 =     C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS122;
const double crLHS124 =     DN(3,2)*crLHS30;
const double crLHS125 =     DN(0,0)*N[3];
const double crLHS126 =     DN(3,0)*N[0];
const double crLHS127 =     C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS25;
const double crLHS128 =     C(0,4)*DN(0,0) + crLHS28 + crLHS33;
const double crLHS129 =     C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS130 =     C(1,4)*DN(0,1);
const double crLHS131 =     C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS130;
const double crLHS132 =     pow(DN(0,1), 2);
const double crLHS133 =     C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS130;
const double crLHS134 =     C(2,4)*DN(0,2);
const double crLHS135 =     C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS134;
const double crLHS136 =     DN(0,1)*crLHS12;
const double crLHS137 =     DN(0,2)*crLHS136;
const double crLHS138 =     C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS57;
const double crLHS139 =     C(0,4)*DN(1,0) + crLHS60 + crLHS64;
const double crLHS140 =     DN(1,0)*crLHS136;
const double crLHS141 =     C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS142 =     C(1,4)*DN(1,1);
const double crLHS143 =     C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS142;
const double crLHS144 =     DN(1,1)*crLHS136;
const double crLHS145 =     C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS142;
const double crLHS146 =     C(2,4)*DN(1,2);
const double crLHS147 =     C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS146;
const double crLHS148 =     DN(1,2)*crLHS136;
const double crLHS149 =     DN(0,1)*N[1];
const double crLHS150 =     DN(1,1)*N[0];
const double crLHS151 =     1.0*DN(1,1)*crLHS20*crLHS7;
const double crLHS152 =     C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS86;
const double crLHS153 =     C(0,4)*DN(2,0) + crLHS89 + crLHS93;
const double crLHS154 =     DN(2,0)*crLHS136;
const double crLHS155 =     C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS156 =     C(1,4)*DN(2,1);
const double crLHS157 =     C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS156;
const double crLHS158 =     DN(2,1)*crLHS136;
const double crLHS159 =     C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS156;
const double crLHS160 =     C(2,4)*DN(2,2);
const double crLHS161 =     C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS160;
const double crLHS162 =     DN(2,2)*crLHS136;
const double crLHS163 =     DN(0,1)*N[2];
const double crLHS164 =     DN(2,1)*N[0];
const double crLHS165 =     1.0*DN(2,1)*crLHS20*crLHS7;
const double crLHS166 =     C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS113;
const double crLHS167 =     C(0,4)*DN(3,0) + crLHS116 + crLHS120;
const double crLHS168 =     DN(3,0)*crLHS136;
const double crLHS169 =     C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS170 =     C(1,4)*DN(3,1);
const double crLHS171 =     C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS170;
const double crLHS172 =     DN(3,1)*crLHS136;
const double crLHS173 =     C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS170;
const double crLHS174 =     C(2,4)*DN(3,2);
const double crLHS175 =     C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS174;
const double crLHS176 =     DN(3,2)*crLHS136;
const double crLHS177 =     DN(0,1)*N[3];
const double crLHS178 =     DN(3,1)*N[0];
const double crLHS179 =     C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS35;
const double crLHS180 =     C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS134;
const double crLHS181 =     C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS182 =     pow(DN(0,2), 2);
const double crLHS183 =     C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS66;
const double crLHS184 =     DN(0,2)*crLHS12;
const double crLHS185 =     DN(1,0)*crLHS184;
const double crLHS186 =     C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS146;
const double crLHS187 =     DN(1,1)*crLHS184;
const double crLHS188 =     C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS189 =     DN(1,2)*crLHS184;
const double crLHS190 =     DN(0,2)*N[1];
const double crLHS191 =     DN(1,2)*N[0];
const double crLHS192 =     1.0*DN(1,2)*crLHS20*crLHS7;
const double crLHS193 =     C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS95;
const double crLHS194 =     DN(2,0)*crLHS184;
const double crLHS195 =     C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS160;
const double crLHS196 =     DN(2,1)*crLHS184;
const double crLHS197 =     C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS198 =     DN(2,2)*crLHS184;
const double crLHS199 =     DN(0,2)*N[2];
const double crLHS200 =     DN(2,2)*N[0];
const double crLHS201 =     1.0*DN(2,2)*crLHS20*crLHS7;
const double crLHS202 =     C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS122;
const double crLHS203 =     DN(3,0)*crLHS184;
const double crLHS204 =     C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS174;
const double crLHS205 =     DN(3,1)*crLHS184;
const double crLHS206 =     C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS207 =     DN(3,2)*crLHS184;
const double crLHS208 =     DN(0,2)*N[3];
const double crLHS209 =     DN(3,2)*N[0];
const double crLHS210 =     crLHS17*crLHS42;
const double crLHS211 =     crLHS0*(N[0] + crLHS210);
const double crLHS212 =     bdf0*crLHS38*crLHS39;
const double crLHS213 =     1.0*crLHS20;
const double crLHS214 =     1.0*DN(0,0)*crLHS20*crLHS7;
const double crLHS215 =     1.0*DN(0,1)*crLHS20*crLHS7;
const double crLHS216 =     1.0*DN(0,2)*crLHS20*crLHS7;
const double crLHS217 =     1.0*DN(0,0)*crLHS20;
const double crLHS218 =     1.0*DN(0,1)*crLHS20;
const double crLHS219 =     1.0*DN(0,2)*crLHS20;
const double crLHS220 =     crLHS0*(DN(1,0)*crLHS217 + DN(1,1)*crLHS218 + DN(1,2)*crLHS219 + N[1]*crLHS40);
const double crLHS221 =     crLHS0*(DN(2,0)*crLHS217 + DN(2,1)*crLHS218 + DN(2,2)*crLHS219 + N[2]*crLHS40);
const double crLHS222 =     crLHS0*(DN(3,0)*crLHS217 + DN(3,1)*crLHS218 + DN(3,2)*crLHS219 + N[3]*crLHS40);
const double crLHS223 =     N[1]*crLHS7;
const double crLHS224 =     1.0*N[1]*crLHS18*crLHS19*crLHS20;
const double crLHS225 =     1.0*crLHS18*crLHS20*crLHS53;
const double crLHS226 =     crLHS16*crLHS223 + crLHS17*crLHS224 + crLHS17*crLHS225 + crLHS52;
const double crLHS227 =     pow(DN(1,0), 2);
const double crLHS228 =     pow(N[1], 2);
const double crLHS229 =     crLHS14*crLHS228 + crLHS223*crLHS53 + crLHS224*crLHS54 + crLHS225*crLHS54;
const double crLHS230 =     DN(1,0)*crLHS12;
const double crLHS231 =     DN(1,1)*crLHS230;
const double crLHS232 =     DN(1,2)*crLHS230;
const double crLHS233 =     N[1]*bdf0*crLHS38*crLHS39;
const double crLHS234 =     crLHS42*crLHS53;
const double crLHS235 =     crLHS0*(-N[1] + crLHS12*crLHS233 + crLHS223*crLHS41 + crLHS234);
const double crLHS236 =     DN(2,0)*crLHS230;
const double crLHS237 =     N[1]*bdf0*crLHS7;
const double crLHS238 =     N[2]*crLHS237;
const double crLHS239 =     crLHS223*crLHS81 + crLHS224*crLHS83 + crLHS225*crLHS83 + crLHS238;
const double crLHS240 =     DN(2,1)*crLHS230;
const double crLHS241 =     DN(2,2)*crLHS230;
const double crLHS242 =     DN(1,0)*N[2];
const double crLHS243 =     DN(2,0)*N[1];
const double crLHS244 =     DN(3,0)*crLHS230;
const double crLHS245 =     N[3]*crLHS237;
const double crLHS246 =     crLHS108*crLHS223 + crLHS110*crLHS224 + crLHS110*crLHS225 + crLHS245;
const double crLHS247 =     DN(3,1)*crLHS230;
const double crLHS248 =     DN(3,2)*crLHS230;
const double crLHS249 =     DN(1,0)*N[3];
const double crLHS250 =     DN(3,0)*N[1];
const double crLHS251 =     pow(DN(1,1), 2);
const double crLHS252 =     DN(1,1)*crLHS12;
const double crLHS253 =     DN(1,2)*crLHS252;
const double crLHS254 =     DN(2,0)*crLHS252;
const double crLHS255 =     DN(2,1)*crLHS252;
const double crLHS256 =     DN(2,2)*crLHS252;
const double crLHS257 =     DN(1,1)*N[2];
const double crLHS258 =     DN(2,1)*N[1];
const double crLHS259 =     DN(3,0)*crLHS252;
const double crLHS260 =     DN(3,1)*crLHS252;
const double crLHS261 =     DN(3,2)*crLHS252;
const double crLHS262 =     DN(1,1)*N[3];
const double crLHS263 =     DN(3,1)*N[1];
const double crLHS264 =     pow(DN(1,2), 2);
const double crLHS265 =     DN(1,2)*crLHS12;
const double crLHS266 =     DN(2,0)*crLHS265;
const double crLHS267 =     DN(2,1)*crLHS265;
const double crLHS268 =     DN(2,2)*crLHS265;
const double crLHS269 =     DN(1,2)*N[2];
const double crLHS270 =     DN(2,2)*N[1];
const double crLHS271 =     DN(3,0)*crLHS265;
const double crLHS272 =     DN(3,1)*crLHS265;
const double crLHS273 =     DN(3,2)*crLHS265;
const double crLHS274 =     DN(1,2)*N[3];
const double crLHS275 =     DN(3,2)*N[1];
const double crLHS276 =     crLHS42*crLHS54;
const double crLHS277 =     crLHS0*(N[1] + crLHS276);
const double crLHS278 =     1.0*DN(1,0)*crLHS20;
const double crLHS279 =     1.0*DN(1,1)*crLHS20;
const double crLHS280 =     1.0*DN(1,2)*crLHS20;
const double crLHS281 =     crLHS0*(DN(2,0)*crLHS278 + DN(2,1)*crLHS279 + DN(2,2)*crLHS280 + N[2]*crLHS233);
const double crLHS282 =     crLHS0*(DN(3,0)*crLHS278 + DN(3,1)*crLHS279 + DN(3,2)*crLHS280 + N[3]*crLHS233);
const double crLHS283 =     N[2]*crLHS7;
const double crLHS284 =     1.0*N[2]*crLHS18*crLHS19*crLHS20;
const double crLHS285 =     1.0*crLHS18*crLHS20*crLHS81;
const double crLHS286 =     crLHS16*crLHS283 + crLHS17*crLHS284 + crLHS17*crLHS285 + crLHS80;
const double crLHS287 =     crLHS238 + crLHS283*crLHS53 + crLHS284*crLHS54 + crLHS285*crLHS54;
const double crLHS288 =     pow(DN(2,0), 2);
const double crLHS289 =     pow(N[2], 2);
const double crLHS290 =     crLHS14*crLHS289 + crLHS283*crLHS81 + crLHS284*crLHS83 + crLHS285*crLHS83;
const double crLHS291 =     DN(2,0)*crLHS12;
const double crLHS292 =     DN(2,1)*crLHS291;
const double crLHS293 =     DN(2,2)*crLHS291;
const double crLHS294 =     crLHS12*crLHS38*crLHS39;
const double crLHS295 =     crLHS42*crLHS81;
const double crLHS296 =     crLHS0*(-N[2] + crLHS283*crLHS41 + crLHS294*crLHS82 + crLHS295);
const double crLHS297 =     DN(3,0)*crLHS291;
const double crLHS298 =     N[2]*N[3]*bdf0;
const double crLHS299 =     crLHS298*crLHS7;
const double crLHS300 =     crLHS108*crLHS283 + crLHS110*crLHS284 + crLHS110*crLHS285 + crLHS299;
const double crLHS301 =     DN(3,1)*crLHS291;
const double crLHS302 =     DN(3,2)*crLHS291;
const double crLHS303 =     DN(2,0)*N[3];
const double crLHS304 =     DN(3,0)*N[2];
const double crLHS305 =     pow(DN(2,1), 2);
const double crLHS306 =     DN(2,1)*crLHS12;
const double crLHS307 =     DN(2,2)*crLHS306;
const double crLHS308 =     DN(3,0)*crLHS306;
const double crLHS309 =     DN(3,1)*crLHS306;
const double crLHS310 =     DN(3,2)*crLHS306;
const double crLHS311 =     DN(2,1)*N[3];
const double crLHS312 =     DN(3,1)*N[2];
const double crLHS313 =     pow(DN(2,2), 2);
const double crLHS314 =     DN(2,2)*crLHS12;
const double crLHS315 =     DN(3,0)*crLHS314;
const double crLHS316 =     DN(3,1)*crLHS314;
const double crLHS317 =     DN(3,2)*crLHS314;
const double crLHS318 =     DN(2,2)*N[3];
const double crLHS319 =     DN(3,2)*N[2];
const double crLHS320 =     crLHS42*crLHS83;
const double crLHS321 =     crLHS0*(N[2] + crLHS320);
const double crLHS322 =     crLHS0*(1.0*DN(2,0)*DN(3,0)*crLHS20 + 1.0*DN(2,1)*DN(3,1)*crLHS20 + 1.0*DN(2,2)*DN(3,2)*crLHS20 + crLHS298*crLHS38*crLHS39);
const double crLHS323 =     N[3]*crLHS7;
const double crLHS324 =     1.0*N[3]*crLHS18*crLHS19*crLHS20;
const double crLHS325 =     1.0*crLHS108*crLHS18*crLHS20;
const double crLHS326 =     crLHS107 + crLHS16*crLHS323 + crLHS17*crLHS324 + crLHS17*crLHS325;
const double crLHS327 =     crLHS245 + crLHS323*crLHS53 + crLHS324*crLHS54 + crLHS325*crLHS54;
const double crLHS328 =     crLHS299 + crLHS323*crLHS81 + crLHS324*crLHS83 + crLHS325*crLHS83;
const double crLHS329 =     pow(DN(3,0), 2);
const double crLHS330 =     pow(N[3], 2);
const double crLHS331 =     crLHS108*crLHS323 + crLHS110*crLHS324 + crLHS110*crLHS325 + crLHS14*crLHS330;
const double crLHS332 =     DN(3,0)*crLHS12;
const double crLHS333 =     DN(3,1)*crLHS332;
const double crLHS334 =     DN(3,2)*crLHS332;
const double crLHS335 =     crLHS0*(-N[3] + crLHS108*crLHS42 + crLHS109*crLHS294 + crLHS323*crLHS41);
const double crLHS336 =     pow(DN(3,1), 2);
const double crLHS337 =     DN(3,1)*DN(3,2)*crLHS12;
const double crLHS338 =     pow(DN(3,2), 2);
const double crLHS339 =     crLHS0*(N[3] + crLHS110*crLHS42);
    rLHS(0,0)+=crLHS0*(DN(0,0)*crLHS1 + DN(0,1)*crLHS3 + DN(0,2)*crLHS5 + crLHS12*crLHS6 + crLHS23);
    rLHS(0,1)+=crLHS0*(DN(0,0)*crLHS24 + DN(0,1)*crLHS26 + DN(0,2)*crLHS29 + crLHS31);
    rLHS(0,2)+=crLHS0*(DN(0,0)*crLHS32 + DN(0,1)*crLHS34 + DN(0,2)*crLHS36 + crLHS37);
    rLHS(0,3)+=DN(0,0)*crLHS44;
    rLHS(0,4)+=crLHS0*(DN(0,0)*crLHS45 + DN(0,1)*crLHS47 + DN(0,2)*crLHS49 + crLHS50 + crLHS55);
    rLHS(0,5)+=crLHS0*(DN(0,0)*crLHS56 + DN(0,1)*crLHS58 + DN(0,2)*crLHS61 + crLHS62);
    rLHS(0,6)+=crLHS0*(DN(0,0)*crLHS63 + DN(0,1)*crLHS65 + DN(0,2)*crLHS67 + crLHS68);
    rLHS(0,7)+=crLHS0*(crLHS16*crLHS73 + crLHS69*crLHS70 - crLHS69 + crLHS71*crLHS72);
    rLHS(0,8)+=crLHS0*(DN(0,0)*crLHS74 + DN(0,1)*crLHS76 + DN(0,2)*crLHS78 + crLHS79 + crLHS84);
    rLHS(0,9)+=crLHS0*(DN(0,0)*crLHS85 + DN(0,1)*crLHS87 + DN(0,2)*crLHS90 + crLHS91);
    rLHS(0,10)+=crLHS0*(DN(0,0)*crLHS92 + DN(0,1)*crLHS94 + DN(0,2)*crLHS96 + crLHS97);
    rLHS(0,11)+=crLHS0*(crLHS100*crLHS16 + crLHS70*crLHS98 + crLHS72*crLHS99 - crLHS98);
    rLHS(0,12)+=crLHS0*(DN(0,0)*crLHS101 + DN(0,1)*crLHS103 + DN(0,2)*crLHS105 + crLHS106 + crLHS111);
    rLHS(0,13)+=crLHS0*(DN(0,0)*crLHS112 + DN(0,1)*crLHS114 + DN(0,2)*crLHS117 + crLHS118);
    rLHS(0,14)+=crLHS0*(DN(0,0)*crLHS119 + DN(0,1)*crLHS121 + DN(0,2)*crLHS123 + crLHS124);
    rLHS(0,15)+=crLHS0*(DN(3,0)*crLHS43 + crLHS125*crLHS70 - crLHS125 + crLHS126*crLHS72);
    rLHS(1,0)+=crLHS0*(DN(0,0)*crLHS3 + DN(0,1)*crLHS127 + DN(0,2)*crLHS128 + crLHS31);
    rLHS(1,1)+=crLHS0*(DN(0,0)*crLHS26 + DN(0,1)*crLHS129 + DN(0,2)*crLHS131 + crLHS12*crLHS132 + crLHS23);
    rLHS(1,2)+=crLHS0*(DN(0,0)*crLHS34 + DN(0,1)*crLHS133 + DN(0,2)*crLHS135 + crLHS137);
    rLHS(1,3)+=DN(0,1)*crLHS44;
    rLHS(1,4)+=crLHS0*(DN(0,0)*crLHS47 + DN(0,1)*crLHS138 + DN(0,2)*crLHS139 + crLHS140);
    rLHS(1,5)+=crLHS0*(DN(0,0)*crLHS58 + DN(0,1)*crLHS141 + DN(0,2)*crLHS143 + crLHS144 + crLHS55);
    rLHS(1,6)+=crLHS0*(DN(0,0)*crLHS65 + DN(0,1)*crLHS145 + DN(0,2)*crLHS147 + crLHS148);
    rLHS(1,7)+=crLHS0*(crLHS149*crLHS70 - crLHS149 + crLHS150*crLHS72 + crLHS151*crLHS16);
    rLHS(1,8)+=crLHS0*(DN(0,0)*crLHS76 + DN(0,1)*crLHS152 + DN(0,2)*crLHS153 + crLHS154);
    rLHS(1,9)+=crLHS0*(DN(0,0)*crLHS87 + DN(0,1)*crLHS155 + DN(0,2)*crLHS157 + crLHS158 + crLHS84);
    rLHS(1,10)+=crLHS0*(DN(0,0)*crLHS94 + DN(0,1)*crLHS159 + DN(0,2)*crLHS161 + crLHS162);
    rLHS(1,11)+=crLHS0*(crLHS16*crLHS165 + crLHS163*crLHS70 - crLHS163 + crLHS164*crLHS72);
    rLHS(1,12)+=crLHS0*(DN(0,0)*crLHS103 + DN(0,1)*crLHS166 + DN(0,2)*crLHS167 + crLHS168);
    rLHS(1,13)+=crLHS0*(DN(0,0)*crLHS114 + DN(0,1)*crLHS169 + DN(0,2)*crLHS171 + crLHS111 + crLHS172);
    rLHS(1,14)+=crLHS0*(DN(0,0)*crLHS121 + DN(0,1)*crLHS173 + DN(0,2)*crLHS175 + crLHS176);
    rLHS(1,15)+=crLHS0*(DN(3,1)*crLHS43 + crLHS177*crLHS70 - crLHS177 + crLHS178*crLHS72);
    rLHS(2,0)+=crLHS0*(DN(0,0)*crLHS5 + DN(0,1)*crLHS128 + DN(0,2)*crLHS179 + crLHS37);
    rLHS(2,1)+=crLHS0*(DN(0,0)*crLHS29 + DN(0,1)*crLHS131 + DN(0,2)*crLHS180 + crLHS137);
    rLHS(2,2)+=crLHS0*(DN(0,0)*crLHS36 + DN(0,1)*crLHS135 + DN(0,2)*crLHS181 + crLHS12*crLHS182 + crLHS23);
    rLHS(2,3)+=DN(0,2)*crLHS44;
    rLHS(2,4)+=crLHS0*(DN(0,0)*crLHS49 + DN(0,1)*crLHS139 + DN(0,2)*crLHS183 + crLHS185);
    rLHS(2,5)+=crLHS0*(DN(0,0)*crLHS61 + DN(0,1)*crLHS143 + DN(0,2)*crLHS186 + crLHS187);
    rLHS(2,6)+=crLHS0*(DN(0,0)*crLHS67 + DN(0,1)*crLHS147 + DN(0,2)*crLHS188 + crLHS189 + crLHS55);
    rLHS(2,7)+=crLHS0*(crLHS16*crLHS192 + crLHS190*crLHS70 - crLHS190 + crLHS191*crLHS72);
    rLHS(2,8)+=crLHS0*(DN(0,0)*crLHS78 + DN(0,1)*crLHS153 + DN(0,2)*crLHS193 + crLHS194);
    rLHS(2,9)+=crLHS0*(DN(0,0)*crLHS90 + DN(0,1)*crLHS157 + DN(0,2)*crLHS195 + crLHS196);
    rLHS(2,10)+=crLHS0*(DN(0,0)*crLHS96 + DN(0,1)*crLHS161 + DN(0,2)*crLHS197 + crLHS198 + crLHS84);
    rLHS(2,11)+=crLHS0*(crLHS16*crLHS201 + crLHS199*crLHS70 - crLHS199 + crLHS200*crLHS72);
    rLHS(2,12)+=crLHS0*(DN(0,0)*crLHS105 + DN(0,1)*crLHS167 + DN(0,2)*crLHS202 + crLHS203);
    rLHS(2,13)+=crLHS0*(DN(0,0)*crLHS117 + DN(0,1)*crLHS171 + DN(0,2)*crLHS204 + crLHS205);
    rLHS(2,14)+=crLHS0*(DN(0,0)*crLHS123 + DN(0,1)*crLHS175 + DN(0,2)*crLHS206 + crLHS111 + crLHS207);
    rLHS(2,15)+=crLHS0*(DN(3,2)*crLHS43 + crLHS208*crLHS70 - crLHS208 + crLHS209*crLHS72);
    rLHS(3,0)+=DN(0,0)*crLHS211;
    rLHS(3,1)+=DN(0,1)*crLHS211;
    rLHS(3,2)+=DN(0,2)*crLHS211;
    rLHS(3,3)+=crLHS0*(crLHS13*crLHS212 + crLHS132*crLHS213 + crLHS182*crLHS213 + crLHS213*crLHS6);
    rLHS(3,4)+=crLHS0*(crLHS214*crLHS54 + crLHS71);
    rLHS(3,5)+=crLHS0*(crLHS150 + crLHS215*crLHS54);
    rLHS(3,6)+=crLHS0*(crLHS191 + crLHS216*crLHS54);
    rLHS(3,7)+=crLHS220;
    rLHS(3,8)+=crLHS0*(crLHS214*crLHS83 + crLHS99);
    rLHS(3,9)+=crLHS0*(crLHS164 + crLHS215*crLHS83);
    rLHS(3,10)+=crLHS0*(crLHS200 + crLHS216*crLHS83);
    rLHS(3,11)+=crLHS221;
    rLHS(3,12)+=crLHS0*(crLHS110*crLHS214 + crLHS126);
    rLHS(3,13)+=crLHS0*(crLHS110*crLHS215 + crLHS178);
    rLHS(3,14)+=crLHS0*(crLHS110*crLHS216 + crLHS209);
    rLHS(3,15)+=crLHS222;
    rLHS(4,0)+=crLHS0*(DN(1,0)*crLHS1 + DN(1,1)*crLHS3 + DN(1,2)*crLHS5 + crLHS226 + crLHS50);
    rLHS(4,1)+=crLHS0*(DN(1,0)*crLHS24 + DN(1,1)*crLHS26 + DN(1,2)*crLHS29 + crLHS140);
    rLHS(4,2)+=crLHS0*(DN(1,0)*crLHS32 + DN(1,1)*crLHS34 + DN(1,2)*crLHS36 + crLHS185);
    rLHS(4,3)+=crLHS0*(crLHS214*crLHS53 + crLHS69*crLHS72 + crLHS70*crLHS71 - crLHS71);
    rLHS(4,4)+=crLHS0*(DN(1,0)*crLHS45 + DN(1,1)*crLHS47 + DN(1,2)*crLHS49 + crLHS12*crLHS227 + crLHS229);
    rLHS(4,5)+=crLHS0*(DN(1,0)*crLHS56 + DN(1,1)*crLHS58 + DN(1,2)*crLHS61 + crLHS231);
    rLHS(4,6)+=crLHS0*(DN(1,0)*crLHS63 + DN(1,1)*crLHS65 + DN(1,2)*crLHS67 + crLHS232);
    rLHS(4,7)+=DN(1,0)*crLHS235;
    rLHS(4,8)+=crLHS0*(DN(1,0)*crLHS74 + DN(1,1)*crLHS76 + DN(1,2)*crLHS78 + crLHS236 + crLHS239);
    rLHS(4,9)+=crLHS0*(DN(1,0)*crLHS85 + DN(1,1)*crLHS87 + DN(1,2)*crLHS90 + crLHS240);
    rLHS(4,10)+=crLHS0*(DN(1,0)*crLHS92 + DN(1,1)*crLHS94 + DN(1,2)*crLHS96 + crLHS241);
    rLHS(4,11)+=crLHS0*(crLHS100*crLHS53 + crLHS242*crLHS70 - crLHS242 + crLHS243*crLHS72);
    rLHS(4,12)+=crLHS0*(DN(1,0)*crLHS101 + DN(1,1)*crLHS103 + DN(1,2)*crLHS105 + crLHS244 + crLHS246);
    rLHS(4,13)+=crLHS0*(DN(1,0)*crLHS112 + DN(1,1)*crLHS114 + DN(1,2)*crLHS117 + crLHS247);
    rLHS(4,14)+=crLHS0*(DN(1,0)*crLHS119 + DN(1,1)*crLHS121 + DN(1,2)*crLHS123 + crLHS248);
    rLHS(4,15)+=crLHS0*(DN(3,0)*crLHS234 + crLHS249*crLHS70 - crLHS249 + crLHS250*crLHS72);
    rLHS(5,0)+=crLHS0*(DN(1,0)*crLHS3 + DN(1,1)*crLHS127 + DN(1,2)*crLHS128 + crLHS62);
    rLHS(5,1)+=crLHS0*(DN(1,0)*crLHS26 + DN(1,1)*crLHS129 + DN(1,2)*crLHS131 + crLHS144 + crLHS226);
    rLHS(5,2)+=crLHS0*(DN(1,0)*crLHS34 + DN(1,1)*crLHS133 + DN(1,2)*crLHS135 + crLHS187);
    rLHS(5,3)+=crLHS0*(crLHS149*crLHS72 + crLHS150*crLHS70 - crLHS150 + crLHS215*crLHS53);
    rLHS(5,4)+=crLHS0*(DN(1,0)*crLHS47 + DN(1,1)*crLHS138 + DN(1,2)*crLHS139 + crLHS231);
    rLHS(5,5)+=crLHS0*(DN(1,0)*crLHS58 + DN(1,1)*crLHS141 + DN(1,2)*crLHS143 + crLHS12*crLHS251 + crLHS229);
    rLHS(5,6)+=crLHS0*(DN(1,0)*crLHS65 + DN(1,1)*crLHS145 + DN(1,2)*crLHS147 + crLHS253);
    rLHS(5,7)+=DN(1,1)*crLHS235;
    rLHS(5,8)+=crLHS0*(DN(1,0)*crLHS76 + DN(1,1)*crLHS152 + DN(1,2)*crLHS153 + crLHS254);
    rLHS(5,9)+=crLHS0*(DN(1,0)*crLHS87 + DN(1,1)*crLHS155 + DN(1,2)*crLHS157 + crLHS239 + crLHS255);
    rLHS(5,10)+=crLHS0*(DN(1,0)*crLHS94 + DN(1,1)*crLHS159 + DN(1,2)*crLHS161 + crLHS256);
    rLHS(5,11)+=crLHS0*(crLHS165*crLHS53 + crLHS257*crLHS70 - crLHS257 + crLHS258*crLHS72);
    rLHS(5,12)+=crLHS0*(DN(1,0)*crLHS103 + DN(1,1)*crLHS166 + DN(1,2)*crLHS167 + crLHS259);
    rLHS(5,13)+=crLHS0*(DN(1,0)*crLHS114 + DN(1,1)*crLHS169 + DN(1,2)*crLHS171 + crLHS246 + crLHS260);
    rLHS(5,14)+=crLHS0*(DN(1,0)*crLHS121 + DN(1,1)*crLHS173 + DN(1,2)*crLHS175 + crLHS261);
    rLHS(5,15)+=crLHS0*(DN(3,1)*crLHS234 + crLHS262*crLHS70 - crLHS262 + crLHS263*crLHS72);
    rLHS(6,0)+=crLHS0*(DN(1,0)*crLHS5 + DN(1,1)*crLHS128 + DN(1,2)*crLHS179 + crLHS68);
    rLHS(6,1)+=crLHS0*(DN(1,0)*crLHS29 + DN(1,1)*crLHS131 + DN(1,2)*crLHS180 + crLHS148);
    rLHS(6,2)+=crLHS0*(DN(1,0)*crLHS36 + DN(1,1)*crLHS135 + DN(1,2)*crLHS181 + crLHS189 + crLHS226);
    rLHS(6,3)+=crLHS0*(crLHS190*crLHS72 + crLHS191*crLHS70 - crLHS191 + crLHS216*crLHS53);
    rLHS(6,4)+=crLHS0*(DN(1,0)*crLHS49 + DN(1,1)*crLHS139 + DN(1,2)*crLHS183 + crLHS232);
    rLHS(6,5)+=crLHS0*(DN(1,0)*crLHS61 + DN(1,1)*crLHS143 + DN(1,2)*crLHS186 + crLHS253);
    rLHS(6,6)+=crLHS0*(DN(1,0)*crLHS67 + DN(1,1)*crLHS147 + DN(1,2)*crLHS188 + crLHS12*crLHS264 + crLHS229);
    rLHS(6,7)+=DN(1,2)*crLHS235;
    rLHS(6,8)+=crLHS0*(DN(1,0)*crLHS78 + DN(1,1)*crLHS153 + DN(1,2)*crLHS193 + crLHS266);
    rLHS(6,9)+=crLHS0*(DN(1,0)*crLHS90 + DN(1,1)*crLHS157 + DN(1,2)*crLHS195 + crLHS267);
    rLHS(6,10)+=crLHS0*(DN(1,0)*crLHS96 + DN(1,1)*crLHS161 + DN(1,2)*crLHS197 + crLHS239 + crLHS268);
    rLHS(6,11)+=crLHS0*(crLHS201*crLHS53 + crLHS269*crLHS70 - crLHS269 + crLHS270*crLHS72);
    rLHS(6,12)+=crLHS0*(DN(1,0)*crLHS105 + DN(1,1)*crLHS167 + DN(1,2)*crLHS202 + crLHS271);
    rLHS(6,13)+=crLHS0*(DN(1,0)*crLHS117 + DN(1,1)*crLHS171 + DN(1,2)*crLHS204 + crLHS272);
    rLHS(6,14)+=crLHS0*(DN(1,0)*crLHS123 + DN(1,1)*crLHS175 + DN(1,2)*crLHS206 + crLHS246 + crLHS273);
    rLHS(6,15)+=crLHS0*(DN(3,2)*crLHS234 + crLHS274*crLHS70 - crLHS274 + crLHS275*crLHS72);
    rLHS(7,0)+=crLHS0*(crLHS17*crLHS73 + crLHS69);
    rLHS(7,1)+=crLHS0*(crLHS149 + crLHS151*crLHS17);
    rLHS(7,2)+=crLHS0*(crLHS17*crLHS192 + crLHS190);
    rLHS(7,3)+=crLHS220;
    rLHS(7,4)+=DN(1,0)*crLHS277;
    rLHS(7,5)+=DN(1,1)*crLHS277;
    rLHS(7,6)+=DN(1,2)*crLHS277;
    rLHS(7,7)+=crLHS0*(crLHS212*crLHS228 + crLHS213*crLHS227 + crLHS213*crLHS251 + crLHS213*crLHS264);
    rLHS(7,8)+=crLHS0*(crLHS243 + crLHS73*crLHS83);
    rLHS(7,9)+=crLHS0*(crLHS151*crLHS83 + crLHS258);
    rLHS(7,10)+=crLHS0*(crLHS192*crLHS83 + crLHS270);
    rLHS(7,11)+=crLHS281;
    rLHS(7,12)+=crLHS0*(crLHS110*crLHS73 + crLHS250);
    rLHS(7,13)+=crLHS0*(crLHS110*crLHS151 + crLHS263);
    rLHS(7,14)+=crLHS0*(crLHS110*crLHS192 + crLHS275);
    rLHS(7,15)+=crLHS282;
    rLHS(8,0)+=crLHS0*(DN(2,0)*crLHS1 + DN(2,1)*crLHS3 + DN(2,2)*crLHS5 + crLHS286 + crLHS79);
    rLHS(8,1)+=crLHS0*(DN(2,0)*crLHS24 + DN(2,1)*crLHS26 + DN(2,2)*crLHS29 + crLHS154);
    rLHS(8,2)+=crLHS0*(DN(2,0)*crLHS32 + DN(2,1)*crLHS34 + DN(2,2)*crLHS36 + crLHS194);
    rLHS(8,3)+=crLHS0*(crLHS214*crLHS81 + crLHS70*crLHS99 + crLHS72*crLHS98 - crLHS99);
    rLHS(8,4)+=crLHS0*(DN(2,0)*crLHS45 + DN(2,1)*crLHS47 + DN(2,2)*crLHS49 + crLHS236 + crLHS287);
    rLHS(8,5)+=crLHS0*(DN(2,0)*crLHS56 + DN(2,1)*crLHS58 + DN(2,2)*crLHS61 + crLHS254);
    rLHS(8,6)+=crLHS0*(DN(2,0)*crLHS63 + DN(2,1)*crLHS65 + DN(2,2)*crLHS67 + crLHS266);
    rLHS(8,7)+=crLHS0*(crLHS242*crLHS72 + crLHS243*crLHS70 - crLHS243 + crLHS73*crLHS81);
    rLHS(8,8)+=crLHS0*(DN(2,0)*crLHS74 + DN(2,1)*crLHS76 + DN(2,2)*crLHS78 + crLHS12*crLHS288 + crLHS290);
    rLHS(8,9)+=crLHS0*(DN(2,0)*crLHS85 + DN(2,1)*crLHS87 + DN(2,2)*crLHS90 + crLHS292);
    rLHS(8,10)+=crLHS0*(DN(2,0)*crLHS92 + DN(2,1)*crLHS94 + DN(2,2)*crLHS96 + crLHS293);
    rLHS(8,11)+=DN(2,0)*crLHS296;
    rLHS(8,12)+=crLHS0*(DN(2,0)*crLHS101 + DN(2,1)*crLHS103 + DN(2,2)*crLHS105 + crLHS297 + crLHS300);
    rLHS(8,13)+=crLHS0*(DN(2,0)*crLHS112 + DN(2,1)*crLHS114 + DN(2,2)*crLHS117 + crLHS301);
    rLHS(8,14)+=crLHS0*(DN(2,0)*crLHS119 + DN(2,1)*crLHS121 + DN(2,2)*crLHS123 + crLHS302);
    rLHS(8,15)+=crLHS0*(DN(3,0)*crLHS295 + crLHS303*crLHS70 - crLHS303 + crLHS304*crLHS72);
    rLHS(9,0)+=crLHS0*(DN(2,0)*crLHS3 + DN(2,1)*crLHS127 + DN(2,2)*crLHS128 + crLHS91);
    rLHS(9,1)+=crLHS0*(DN(2,0)*crLHS26 + DN(2,1)*crLHS129 + DN(2,2)*crLHS131 + crLHS158 + crLHS286);
    rLHS(9,2)+=crLHS0*(DN(2,0)*crLHS34 + DN(2,1)*crLHS133 + DN(2,2)*crLHS135 + crLHS196);
    rLHS(9,3)+=crLHS0*(crLHS163*crLHS72 + crLHS164*crLHS70 - crLHS164 + crLHS215*crLHS81);
    rLHS(9,4)+=crLHS0*(DN(2,0)*crLHS47 + DN(2,1)*crLHS138 + DN(2,2)*crLHS139 + crLHS240);
    rLHS(9,5)+=crLHS0*(DN(2,0)*crLHS58 + DN(2,1)*crLHS141 + DN(2,2)*crLHS143 + crLHS255 + crLHS287);
    rLHS(9,6)+=crLHS0*(DN(2,0)*crLHS65 + DN(2,1)*crLHS145 + DN(2,2)*crLHS147 + crLHS267);
    rLHS(9,7)+=crLHS0*(crLHS151*crLHS81 + crLHS257*crLHS72 + crLHS258*crLHS70 - crLHS258);
    rLHS(9,8)+=crLHS0*(DN(2,0)*crLHS76 + DN(2,1)*crLHS152 + DN(2,2)*crLHS153 + crLHS292);
    rLHS(9,9)+=crLHS0*(DN(2,0)*crLHS87 + DN(2,1)*crLHS155 + DN(2,2)*crLHS157 + crLHS12*crLHS305 + crLHS290);
    rLHS(9,10)+=crLHS0*(DN(2,0)*crLHS94 + DN(2,1)*crLHS159 + DN(2,2)*crLHS161 + crLHS307);
    rLHS(9,11)+=DN(2,1)*crLHS296;
    rLHS(9,12)+=crLHS0*(DN(2,0)*crLHS103 + DN(2,1)*crLHS166 + DN(2,2)*crLHS167 + crLHS308);
    rLHS(9,13)+=crLHS0*(DN(2,0)*crLHS114 + DN(2,1)*crLHS169 + DN(2,2)*crLHS171 + crLHS300 + crLHS309);
    rLHS(9,14)+=crLHS0*(DN(2,0)*crLHS121 + DN(2,1)*crLHS173 + DN(2,2)*crLHS175 + crLHS310);
    rLHS(9,15)+=crLHS0*(DN(3,1)*crLHS295 + crLHS311*crLHS70 - crLHS311 + crLHS312*crLHS72);
    rLHS(10,0)+=crLHS0*(DN(2,0)*crLHS5 + DN(2,1)*crLHS128 + DN(2,2)*crLHS179 + crLHS97);
    rLHS(10,1)+=crLHS0*(DN(2,0)*crLHS29 + DN(2,1)*crLHS131 + DN(2,2)*crLHS180 + crLHS162);
    rLHS(10,2)+=crLHS0*(DN(2,0)*crLHS36 + DN(2,1)*crLHS135 + DN(2,2)*crLHS181 + crLHS198 + crLHS286);
    rLHS(10,3)+=crLHS0*(crLHS199*crLHS72 + crLHS200*crLHS70 - crLHS200 + crLHS216*crLHS81);
    rLHS(10,4)+=crLHS0*(DN(2,0)*crLHS49 + DN(2,1)*crLHS139 + DN(2,2)*crLHS183 + crLHS241);
    rLHS(10,5)+=crLHS0*(DN(2,0)*crLHS61 + DN(2,1)*crLHS143 + DN(2,2)*crLHS186 + crLHS256);
    rLHS(10,6)+=crLHS0*(DN(2,0)*crLHS67 + DN(2,1)*crLHS147 + DN(2,2)*crLHS188 + crLHS268 + crLHS287);
    rLHS(10,7)+=crLHS0*(crLHS192*crLHS81 + crLHS269*crLHS72 + crLHS270*crLHS70 - crLHS270);
    rLHS(10,8)+=crLHS0*(DN(2,0)*crLHS78 + DN(2,1)*crLHS153 + DN(2,2)*crLHS193 + crLHS293);
    rLHS(10,9)+=crLHS0*(DN(2,0)*crLHS90 + DN(2,1)*crLHS157 + DN(2,2)*crLHS195 + crLHS307);
    rLHS(10,10)+=crLHS0*(DN(2,0)*crLHS96 + DN(2,1)*crLHS161 + DN(2,2)*crLHS197 + crLHS12*crLHS313 + crLHS290);
    rLHS(10,11)+=DN(2,2)*crLHS296;
    rLHS(10,12)+=crLHS0*(DN(2,0)*crLHS105 + DN(2,1)*crLHS167 + DN(2,2)*crLHS202 + crLHS315);
    rLHS(10,13)+=crLHS0*(DN(2,0)*crLHS117 + DN(2,1)*crLHS171 + DN(2,2)*crLHS204 + crLHS316);
    rLHS(10,14)+=crLHS0*(DN(2,0)*crLHS123 + DN(2,1)*crLHS175 + DN(2,2)*crLHS206 + crLHS300 + crLHS317);
    rLHS(10,15)+=crLHS0*(DN(3,2)*crLHS295 + crLHS318*crLHS70 - crLHS318 + crLHS319*crLHS72);
    rLHS(11,0)+=crLHS0*(crLHS100*crLHS17 + crLHS98);
    rLHS(11,1)+=crLHS0*(crLHS163 + crLHS165*crLHS17);
    rLHS(11,2)+=crLHS0*(crLHS17*crLHS201 + crLHS199);
    rLHS(11,3)+=crLHS221;
    rLHS(11,4)+=crLHS0*(crLHS100*crLHS54 + crLHS242);
    rLHS(11,5)+=crLHS0*(crLHS165*crLHS54 + crLHS257);
    rLHS(11,6)+=crLHS0*(crLHS201*crLHS54 + crLHS269);
    rLHS(11,7)+=crLHS281;
    rLHS(11,8)+=DN(2,0)*crLHS321;
    rLHS(11,9)+=DN(2,1)*crLHS321;
    rLHS(11,10)+=DN(2,2)*crLHS321;
    rLHS(11,11)+=crLHS0*(crLHS212*crLHS289 + crLHS213*crLHS288 + crLHS213*crLHS305 + crLHS213*crLHS313);
    rLHS(11,12)+=crLHS0*(crLHS100*crLHS110 + crLHS304);
    rLHS(11,13)+=crLHS0*(crLHS110*crLHS165 + crLHS312);
    rLHS(11,14)+=crLHS0*(crLHS110*crLHS201 + crLHS319);
    rLHS(11,15)+=crLHS322;
    rLHS(12,0)+=crLHS0*(DN(3,0)*crLHS1 + DN(3,1)*crLHS3 + DN(3,2)*crLHS5 + crLHS106 + crLHS326);
    rLHS(12,1)+=crLHS0*(DN(3,0)*crLHS24 + DN(3,1)*crLHS26 + DN(3,2)*crLHS29 + crLHS168);
    rLHS(12,2)+=crLHS0*(DN(3,0)*crLHS32 + DN(3,1)*crLHS34 + DN(3,2)*crLHS36 + crLHS203);
    rLHS(12,3)+=crLHS0*(crLHS108*crLHS214 + crLHS125*crLHS72 + crLHS126*crLHS70 - crLHS126);
    rLHS(12,4)+=crLHS0*(DN(3,0)*crLHS45 + DN(3,1)*crLHS47 + DN(3,2)*crLHS49 + crLHS244 + crLHS327);
    rLHS(12,5)+=crLHS0*(DN(3,0)*crLHS56 + DN(3,1)*crLHS58 + DN(3,2)*crLHS61 + crLHS259);
    rLHS(12,6)+=crLHS0*(DN(3,0)*crLHS63 + DN(3,1)*crLHS65 + DN(3,2)*crLHS67 + crLHS271);
    rLHS(12,7)+=crLHS0*(crLHS108*crLHS73 + crLHS249*crLHS72 + crLHS250*crLHS70 - crLHS250);
    rLHS(12,8)+=crLHS0*(DN(3,0)*crLHS74 + DN(3,1)*crLHS76 + DN(3,2)*crLHS78 + crLHS297 + crLHS328);
    rLHS(12,9)+=crLHS0*(DN(3,0)*crLHS85 + DN(3,1)*crLHS87 + DN(3,2)*crLHS90 + crLHS308);
    rLHS(12,10)+=crLHS0*(DN(3,0)*crLHS92 + DN(3,1)*crLHS94 + DN(3,2)*crLHS96 + crLHS315);
    rLHS(12,11)+=crLHS0*(crLHS100*crLHS108 + crLHS303*crLHS72 + crLHS304*crLHS70 - crLHS304);
    rLHS(12,12)+=crLHS0*(DN(3,0)*crLHS101 + DN(3,1)*crLHS103 + DN(3,2)*crLHS105 + crLHS12*crLHS329 + crLHS331);
    rLHS(12,13)+=crLHS0*(DN(3,0)*crLHS112 + DN(3,1)*crLHS114 + DN(3,2)*crLHS117 + crLHS333);
    rLHS(12,14)+=crLHS0*(DN(3,0)*crLHS119 + DN(3,1)*crLHS121 + DN(3,2)*crLHS123 + crLHS334);
    rLHS(12,15)+=DN(3,0)*crLHS335;
    rLHS(13,0)+=crLHS0*(DN(3,0)*crLHS3 + DN(3,1)*crLHS127 + DN(3,2)*crLHS128 + crLHS118);
    rLHS(13,1)+=crLHS0*(DN(3,0)*crLHS26 + DN(3,1)*crLHS129 + DN(3,2)*crLHS131 + crLHS172 + crLHS326);
    rLHS(13,2)+=crLHS0*(DN(3,0)*crLHS34 + DN(3,1)*crLHS133 + DN(3,2)*crLHS135 + crLHS205);
    rLHS(13,3)+=crLHS0*(crLHS108*crLHS215 + crLHS177*crLHS72 + crLHS178*crLHS70 - crLHS178);
    rLHS(13,4)+=crLHS0*(DN(3,0)*crLHS47 + DN(3,1)*crLHS138 + DN(3,2)*crLHS139 + crLHS247);
    rLHS(13,5)+=crLHS0*(DN(3,0)*crLHS58 + DN(3,1)*crLHS141 + DN(3,2)*crLHS143 + crLHS260 + crLHS327);
    rLHS(13,6)+=crLHS0*(DN(3,0)*crLHS65 + DN(3,1)*crLHS145 + DN(3,2)*crLHS147 + crLHS272);
    rLHS(13,7)+=crLHS0*(crLHS108*crLHS151 + crLHS262*crLHS72 + crLHS263*crLHS70 - crLHS263);
    rLHS(13,8)+=crLHS0*(DN(3,0)*crLHS76 + DN(3,1)*crLHS152 + DN(3,2)*crLHS153 + crLHS301);
    rLHS(13,9)+=crLHS0*(DN(3,0)*crLHS87 + DN(3,1)*crLHS155 + DN(3,2)*crLHS157 + crLHS309 + crLHS328);
    rLHS(13,10)+=crLHS0*(DN(3,0)*crLHS94 + DN(3,1)*crLHS159 + DN(3,2)*crLHS161 + crLHS316);
    rLHS(13,11)+=crLHS0*(crLHS108*crLHS165 + crLHS311*crLHS72 + crLHS312*crLHS70 - crLHS312);
    rLHS(13,12)+=crLHS0*(DN(3,0)*crLHS103 + DN(3,1)*crLHS166 + DN(3,2)*crLHS167 + crLHS333);
    rLHS(13,13)+=crLHS0*(DN(3,0)*crLHS114 + DN(3,1)*crLHS169 + DN(3,2)*crLHS171 + crLHS12*crLHS336 + crLHS331);
    rLHS(13,14)+=crLHS0*(DN(3,0)*crLHS121 + DN(3,1)*crLHS173 + DN(3,2)*crLHS175 + crLHS337);
    rLHS(13,15)+=DN(3,1)*crLHS335;
    rLHS(14,0)+=crLHS0*(DN(3,0)*crLHS5 + DN(3,1)*crLHS128 + DN(3,2)*crLHS179 + crLHS124);
    rLHS(14,1)+=crLHS0*(DN(3,0)*crLHS29 + DN(3,1)*crLHS131 + DN(3,2)*crLHS180 + crLHS176);
    rLHS(14,2)+=crLHS0*(DN(3,0)*crLHS36 + DN(3,1)*crLHS135 + DN(3,2)*crLHS181 + crLHS207 + crLHS326);
    rLHS(14,3)+=crLHS0*(crLHS108*crLHS216 + crLHS208*crLHS72 + crLHS209*crLHS70 - crLHS209);
    rLHS(14,4)+=crLHS0*(DN(3,0)*crLHS49 + DN(3,1)*crLHS139 + DN(3,2)*crLHS183 + crLHS248);
    rLHS(14,5)+=crLHS0*(DN(3,0)*crLHS61 + DN(3,1)*crLHS143 + DN(3,2)*crLHS186 + crLHS261);
    rLHS(14,6)+=crLHS0*(DN(3,0)*crLHS67 + DN(3,1)*crLHS147 + DN(3,2)*crLHS188 + crLHS273 + crLHS327);
    rLHS(14,7)+=crLHS0*(crLHS108*crLHS192 + crLHS274*crLHS72 + crLHS275*crLHS70 - crLHS275);
    rLHS(14,8)+=crLHS0*(DN(3,0)*crLHS78 + DN(3,1)*crLHS153 + DN(3,2)*crLHS193 + crLHS302);
    rLHS(14,9)+=crLHS0*(DN(3,0)*crLHS90 + DN(3,1)*crLHS157 + DN(3,2)*crLHS195 + crLHS310);
    rLHS(14,10)+=crLHS0*(DN(3,0)*crLHS96 + DN(3,1)*crLHS161 + DN(3,2)*crLHS197 + crLHS317 + crLHS328);
    rLHS(14,11)+=crLHS0*(crLHS108*crLHS201 + crLHS318*crLHS72 + crLHS319*crLHS70 - crLHS319);
    rLHS(14,12)+=crLHS0*(DN(3,0)*crLHS105 + DN(3,1)*crLHS167 + DN(3,2)*crLHS202 + crLHS334);
    rLHS(14,13)+=crLHS0*(DN(3,0)*crLHS117 + DN(3,1)*crLHS171 + DN(3,2)*crLHS204 + crLHS337);
    rLHS(14,14)+=crLHS0*(DN(3,0)*crLHS123 + DN(3,1)*crLHS175 + DN(3,2)*crLHS206 + crLHS12*crLHS338 + crLHS331);
    rLHS(14,15)+=DN(3,2)*crLHS335;
    rLHS(15,0)+=crLHS0*(DN(3,0)*crLHS210 + crLHS125);
    rLHS(15,1)+=crLHS0*(DN(3,1)*crLHS210 + crLHS177);
    rLHS(15,2)+=crLHS0*(DN(3,2)*crLHS210 + crLHS208);
    rLHS(15,3)+=crLHS222;
    rLHS(15,4)+=crLHS0*(DN(3,0)*crLHS276 + crLHS249);
    rLHS(15,5)+=crLHS0*(DN(3,1)*crLHS276 + crLHS262);
    rLHS(15,6)+=crLHS0*(DN(3,2)*crLHS276 + crLHS274);
    rLHS(15,7)+=crLHS282;
    rLHS(15,8)+=crLHS0*(DN(3,0)*crLHS320 + crLHS303);
    rLHS(15,9)+=crLHS0*(DN(3,1)*crLHS320 + crLHS311);
    rLHS(15,10)+=crLHS0*(DN(3,2)*crLHS320 + crLHS318);
    rLHS(15,11)+=crLHS322;
    rLHS(15,12)+=DN(3,0)*crLHS339;
    rLHS(15,13)+=DN(3,1)*crLHS339;
    rLHS(15,14)+=DN(3,2)*crLHS339;
    rLHS(15,15)+=crLHS0*(crLHS212*crLHS330 + crLHS213*crLHS329 + crLHS213*crLHS336 + crLHS213*crLHS338);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<2,3>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<2,3>& rData,
    VectorType& rRHS)
{
    const array_1d<double,3>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const array_1d<double,3>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,2,3>& v = rData.Velocity;
    const BoundedMatrix<double,2,3>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,3>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = v - vmesh;
    const BoundedMatrix<double,2,3>& f = rData.BodyForce;
    const array_1d<double,3>& p = rData.Pressure;
    const array_1d<double,3>& pn = rData.Pressure_OldStep1;
    const array_1d<double,3>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,3>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // RHS Gauss point contribution
    const double weight = rData.Weight;

    const double crRHS0 =     N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crRHS1 =     N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2];
const double crRHS2 =     crRHS1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crRHS3 =     crRHS1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crRHS4 =     DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS5 =     N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS6 =     N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS7 =     crRHS1*(crRHS4*crRHS5 + crRHS6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crRHS8 =     crRHS1*stab_c2*sqrt(pow(crRHS5, 2) + pow(crRHS6, 2));
const double crRHS9 =     DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS10 =     crRHS4 + crRHS9;
const double crRHS11 =     1.0/crRHS1;
const double crRHS12 =     crRHS11*(crRHS5*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2]) + crRHS6*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2]));
const double crRHS13 =     crRHS11*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2], 2);
const double crRHS14 =     (crRHS8*h/stab_c1 + mu)*(crRHS10 + crRHS12 + crRHS13);
const double crRHS15 =     DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crRHS16 =     N[0]*crRHS1*crRHS15;
const double crRHS17 =     1.0/(crRHS1*dyn_tau/dt + crRHS8/h + mu*stab_c1/pow(h, 2));
const double crRHS18 =     1.0*crRHS17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS2 + crRHS3 + crRHS7);
const double crRHS19 =     crRHS1*(DN(0,0)*crRHS5 + DN(0,1)*crRHS6);
const double crRHS20 =     crRHS1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS21 =     crRHS1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crRHS22 =     crRHS1*(crRHS5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crRHS6*crRHS9);
const double crRHS23 =     1.0*crRHS17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS20 + crRHS21 + crRHS22);
const double crRHS24 =     N[1]*crRHS1*crRHS15;
const double crRHS25 =     crRHS1*(DN(1,0)*crRHS5 + DN(1,1)*crRHS6);
const double crRHS26 =     N[2]*crRHS1*crRHS15;
const double crRHS27 =     crRHS1*(DN(2,0)*crRHS5 + DN(2,1)*crRHS6);
    rRHS[0]+=-weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS14 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS2 + N[0]*crRHS3 + N[0]*crRHS7 + crRHS16*crRHS18 + crRHS18*crRHS19);
    rRHS[1]+=-weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS14 + DN(0,1)*stress[1] - N[0]*crRHS20 + N[0]*crRHS21 + N[0]*crRHS22 + crRHS16*crRHS23 + crRHS19*crRHS23);
    rRHS[2]+=-weight*(DN(0,0)*crRHS18 + DN(0,1)*crRHS23 + N[0]*crRHS10 + N[0]*crRHS12 + N[0]*crRHS13);
    rRHS[3]+=-weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS14 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS2 + N[1]*crRHS3 + N[1]*crRHS7 + crRHS18*crRHS24 + crRHS18*crRHS25);
    rRHS[4]+=-weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS14 + DN(1,1)*stress[1] - N[1]*crRHS20 + N[1]*crRHS21 + N[1]*crRHS22 + crRHS23*crRHS24 + crRHS23*crRHS25);
    rRHS[5]+=-weight*(DN(1,0)*crRHS18 + DN(1,1)*crRHS23 + N[1]*crRHS10 + N[1]*crRHS12 + N[1]*crRHS13);
    rRHS[6]+=-weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS14 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS2 + N[2]*crRHS3 + N[2]*crRHS7 + crRHS18*crRHS26 + crRHS18*crRHS27);
    rRHS[7]+=-weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS14 + DN(2,1)*stress[1] - N[2]*crRHS20 + N[2]*crRHS21 + N[2]*crRHS22 + crRHS23*crRHS26 + crRHS23*crRHS27);
    rRHS[8]+=-weight*(DN(2,0)*crRHS18 + DN(2,1)*crRHS23 + N[2]*crRHS10 + N[2]*crRHS12 + N[2]*crRHS13);

}

template <>
void WeaklyCompressibleNavierStokes<WeaklyCompressibleNavierStokesData<3,4>>::ComputeGaussPointRHSContribution(
    WeaklyCompressibleNavierStokesData<3,4>& rData,
    VectorType& rRHS)
{
    const array_1d<double,4>& rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const array_1d<double,4>& c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    const BoundedMatrix<double,3,4>& v = rData.Velocity;
    const BoundedMatrix<double,3,4>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,3,4>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = v - vmesh;
    const BoundedMatrix<double,3,4>& f = rData.BodyForce;
    const array_1d<double,4>& p = rData.Pressure;
    const array_1d<double,4>& pn = rData.Pressure_OldStep1;
    const array_1d<double,4>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,6>& stress = rData.ShearStress;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // RHS Gauss point contribution
    const double weight = rData.Weight;

    const double crRHS0 =     N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crRHS1 =     N[0]*rho[0] + N[1]*rho[1] + N[2]*rho[2] + N[3]*rho[3];
const double crRHS2 =     crRHS1*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS3 =     crRHS1*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crRHS4 =     DN(0,0)*v(0,0);
const double crRHS5 =     DN(1,0)*v(1,0);
const double crRHS6 =     DN(2,0)*v(2,0);
const double crRHS7 =     DN(3,0)*v(3,0);
const double crRHS8 =     N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS9 =     N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS10 =     N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS11 =     crRHS1*(crRHS10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crRHS8*(crRHS4 + crRHS5 + crRHS6 + crRHS7) + crRHS9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)));
const double crRHS12 =     crRHS1*stab_c2*sqrt(pow(crRHS10, 2) + pow(crRHS8, 2) + pow(crRHS9, 2));
const double crRHS13 =     DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS14 =     DN(0,1)*v(0,1);
const double crRHS15 =     DN(1,1)*v(1,1);
const double crRHS16 =     DN(2,1)*v(2,1);
const double crRHS17 =     DN(3,1)*v(3,1);
const double crRHS18 =     crRHS13 + crRHS14 + crRHS15 + crRHS16 + crRHS17 + crRHS4 + crRHS5 + crRHS6 + crRHS7;
const double crRHS19 =     1.0/crRHS1;
const double crRHS20 =     crRHS19*(N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(N[0]*c[0] + N[1]*c[1] + N[2]*c[2] + N[3]*c[3], 2);
const double crRHS21 =     crRHS19*(crRHS10*(DN(0,2)*rho[0] + DN(1,2)*rho[1] + DN(2,2)*rho[2] + DN(3,2)*rho[3]) + crRHS8*(DN(0,0)*rho[0] + DN(1,0)*rho[1] + DN(2,0)*rho[2] + DN(3,0)*rho[3]) + crRHS9*(DN(0,1)*rho[0] + DN(1,1)*rho[1] + DN(2,1)*rho[2] + DN(3,1)*rho[3]));
const double crRHS22 =     (crRHS12*h/stab_c1 + mu)*(crRHS18 + crRHS20 + crRHS21);
const double crRHS23 =     DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crRHS24 =     N[0]*crRHS1*crRHS23;
const double crRHS25 =     1.0/(crRHS1*dyn_tau/dt + crRHS12/h + mu*stab_c1/pow(h, 2));
const double crRHS26 =     1.0*crRHS25*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + crRHS11 - crRHS2 + crRHS3);
const double crRHS27 =     crRHS1*(DN(0,0)*crRHS8 + DN(0,1)*crRHS9 + DN(0,2)*crRHS10);
const double crRHS28 =     crRHS1*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS29 =     crRHS1*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crRHS30 =     crRHS1*(crRHS10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crRHS8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crRHS9*(crRHS14 + crRHS15 + crRHS16 + crRHS17));
const double crRHS31 =     1.0*crRHS25*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS28 + crRHS29 + crRHS30);
const double crRHS32 =     crRHS1*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS33 =     crRHS1*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crRHS34 =     crRHS1*(crRHS10*crRHS13 + crRHS8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crRHS9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crRHS35 =     1.0*crRHS25*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS32 + crRHS33 + crRHS34);
const double crRHS36 =     N[1]*crRHS1*crRHS23;
const double crRHS37 =     crRHS1*(DN(1,0)*crRHS8 + DN(1,1)*crRHS9 + DN(1,2)*crRHS10);
const double crRHS38 =     N[2]*crRHS1*crRHS23;
const double crRHS39 =     crRHS1*(DN(2,0)*crRHS8 + DN(2,1)*crRHS9 + DN(2,2)*crRHS10);
const double crRHS40 =     N[3]*crRHS1*crRHS23;
const double crRHS41 =     crRHS1*(DN(3,0)*crRHS8 + DN(3,1)*crRHS9 + DN(3,2)*crRHS10);
    rRHS[0]+=-weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS22 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] + N[0]*crRHS11 - N[0]*crRHS2 + N[0]*crRHS3 + crRHS24*crRHS26 + crRHS26*crRHS27);
    rRHS[1]+=-weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS22 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS28 + N[0]*crRHS29 + N[0]*crRHS30 + crRHS24*crRHS31 + crRHS27*crRHS31);
    rRHS[2]+=-weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS22 + DN(0,2)*stress[2] - N[0]*crRHS32 + N[0]*crRHS33 + N[0]*crRHS34 + crRHS24*crRHS35 + crRHS27*crRHS35);
    rRHS[3]+=-weight*(DN(0,0)*crRHS26 + DN(0,1)*crRHS31 + DN(0,2)*crRHS35 + N[0]*crRHS18 + N[0]*crRHS20 + N[0]*crRHS21);
    rRHS[4]+=-weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS22 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] + N[1]*crRHS11 - N[1]*crRHS2 + N[1]*crRHS3 + crRHS26*crRHS36 + crRHS26*crRHS37);
    rRHS[5]+=-weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS22 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS28 + N[1]*crRHS29 + N[1]*crRHS30 + crRHS31*crRHS36 + crRHS31*crRHS37);
    rRHS[6]+=-weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS22 + DN(1,2)*stress[2] - N[1]*crRHS32 + N[1]*crRHS33 + N[1]*crRHS34 + crRHS35*crRHS36 + crRHS35*crRHS37);
    rRHS[7]+=-weight*(DN(1,0)*crRHS26 + DN(1,1)*crRHS31 + DN(1,2)*crRHS35 + N[1]*crRHS18 + N[1]*crRHS20 + N[1]*crRHS21);
    rRHS[8]+=-weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS22 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] + N[2]*crRHS11 - N[2]*crRHS2 + N[2]*crRHS3 + crRHS26*crRHS38 + crRHS26*crRHS39);
    rRHS[9]+=-weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS22 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS28 + N[2]*crRHS29 + N[2]*crRHS30 + crRHS31*crRHS38 + crRHS31*crRHS39);
    rRHS[10]+=-weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS22 + DN(2,2)*stress[2] - N[2]*crRHS32 + N[2]*crRHS33 + N[2]*crRHS34 + crRHS35*crRHS38 + crRHS35*crRHS39);
    rRHS[11]+=-weight*(DN(2,0)*crRHS26 + DN(2,1)*crRHS31 + DN(2,2)*crRHS35 + N[2]*crRHS18 + N[2]*crRHS20 + N[2]*crRHS21);
    rRHS[12]+=-weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS22 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] + N[3]*crRHS11 - N[3]*crRHS2 + N[3]*crRHS3 + crRHS26*crRHS40 + crRHS26*crRHS41);
    rRHS[13]+=-weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS22 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS28 + N[3]*crRHS29 + N[3]*crRHS30 + crRHS31*crRHS40 + crRHS31*crRHS41);
    rRHS[14]+=-weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS22 + DN(3,2)*stress[2] - N[3]*crRHS32 + N[3]*crRHS33 + N[3]*crRHS34 + crRHS35*crRHS40 + crRHS35*crRHS41);
    rRHS[15]+=-weight*(DN(3,0)*crRHS26 + DN(3,1)*crRHS31 + DN(3,2)*crRHS35 + N[3]*crRHS18 + N[3]*crRHS20 + N[3]*crRHS21);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void WeaklyCompressibleNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
}


template< class TElementData >
void WeaklyCompressibleNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> >;
template class WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> >;

}