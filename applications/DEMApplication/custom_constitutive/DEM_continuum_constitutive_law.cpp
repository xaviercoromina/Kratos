#include "DEM_continuum_constitutive_law.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw() {}

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw) {}

    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw() {}

    std::string DEMContinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) shouldn't be accessed, use derived class instead"<<std::endl;
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMContinuumConstitutiveLaw::Initialize(SphericContinuumParticle* element1,
                                                 SphericContinuumParticle* element2,
                                                 Properties::Pointer pProps) {
        mpProperties = pProps;
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if (verbose) KRATOS_INFO("DEM")  << "Assigning " << pProp->GetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning " << pProp->GetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());

        TransferParametersToProperties(parameters, pProp);
        this->Check(pProp);
    }

    void DEMContinuumConstitutiveLaw::TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp)  {
    }

    void DEMContinuumConstitutiveLaw::Check(Properties::Pointer pProp) const {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::Check) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }

    void DEMContinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
                                                            double ViscoDampingLocalContactForce[3],
                                                            double indentation,
                                                            double equiv_visco_damp_coeff_normal,
                                                            double equiv_visco_damp_coeff_tangential,
                                                            bool& sliding,
                                                            int failure_id) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateViscoDamping) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                                       SphericContinuumParticle* neighbor,
                                                                       double equiv_young,
                                                                       double distance,
                                                                       double calculation_area,
                                                                       double LocalCoordSystem[3][3],
                                                                       double ElasticLocalRotationalMoment[3],
                                                                       double ViscoLocalRotationalMoment[3],
                                                                       double equiv_poisson,
                                                                       double indentation) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMContinuumConstitutiveLaw::AddPoissonContribution(const double equiv_poisson,
                                                             double LocalCoordSystem[3][3],
                                                             double& normal_force,
                                                             double calculation_area,
                                                             BoundedMatrix<double, 3, 3>* mSymmStressTensor,
                                                             SphericContinuumParticle* element1, SphericContinuumParticle* element2,
                                                             const ProcessInfo& r_process_info,
                                                             const int i_neighbor_count,
                                                             const double indentation) {
    }


    double DEMContinuumConstitutiveLaw::LocalMaxSearchDistance(const int i,
                                          SphericContinuumParticle* element1,
                                          SphericContinuumParticle* element2) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::LocalMaxSearchDistance) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    double DEMContinuumConstitutiveLaw::LocalPeriod(const int i, SphericContinuumParticle* element1,
                                                                 SphericContinuumParticle* element2) {

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        double initial_delta = element1->GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants
        double kn_el = equiv_young * calculation_area / initial_dist;

        const double mRealMass = element1->GetMass();  // { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
        const double &other_real_mass = element2->GetMass();

        // calculation of damping gamma

        const double& damping_gamma = (*mpProperties)[DAMPING_GAMMA];

        double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);
        const double viscous_damping_coeff     = 2.0 * damping_gamma * sqrt(equiv_mass * kn_el);
        //double viscous_damping_coeff = (1-equiv_coefficientOfRestitution) * 2.0 * sqrt(kn_el * equiv_mass);

        double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);
        //double a = 1.4142-equiv_gamma*equiv_gamma;

        //double sqr_period = kn_el / equiv_mass - rescaled_damping*rescaled_damping;
        double sqr_period = sqrt(2.0) * kn_el / equiv_mass - rescaled_damping*rescaled_damping;   //esta es la correcta en continuu suponiendo un maximo de Kt= Kn
        return sqr_period;
    }
    
    bool DEMContinuumConstitutiveLaw::CheckRequirementsOfStressTensor() {
        return false;
    }

} //kratos
