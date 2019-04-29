//Authors: M.A. Celigueta, G. Casas and M. Chung To Sang (CIMNE)
//Date: February 2019

#include "DEM_electromagnetic_CL.h"
#include "custom_elements/ion_particle.h"




namespace Kratos {

    void DEM_electromagnetic::Initialize(const ProcessInfo& r_process_info) {}

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_electromagnetic::Clone() const
    {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_electromagnetic(*this));
        return p_clone;
    }

    void DEM_electromagnetic::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        KRATOS_INFO("DEM") << "Assigning DEM_electromagnetic to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    std::string DEM_electromagnetic::GetTypeOfLaw() {
        std::string type_of_law = "Linear";
        return type_of_law;
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_electromagnetic::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation)
    {}

    void DEM_electromagnetic::CalculateForces(const ProcessInfo& r_process_info,
                                                  const double OldLocalContactForce[3],
                                                  double LocalElasticContactForce[3],
                                                  double LocalDeltDisp[3],
                                                  double LocalRelVel[3],
                                                  double indentation,
                                                  double previous_indentation,
                                                  double ViscoDampingLocalContactForce[3],
                                                  double& cohesive_force,
                                                  SphericParticle* element1,
                                                  SphericParticle* element2,
                                                  bool& sliding,
                                                  double LocalCoordSystem[3][3])
    {
        //InitializeContact(element1, element2, indentation);
        //LocalElasticContactForce[2]  = CalculateNormalForce(indentation);
        //KRATOS_INFO("DEM:  inside DEM_electromag CalculateForces")<< std::endl;

        const double distance = element1->GetRadius() + element2->GetRadius() - indentation;
        //KRATOS_INFO("DEM: distance =")<< distance << std::endl;
        //const array_1d<double, 3>& global_force = element2->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
        //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, global_force, LocalElasticContactForce);
        LocalElasticContactForce[0] = 0.0;
        LocalElasticContactForce[1] = 0.0;
        LocalElasticContactForce[2] = CalculateNormalForce(distance);

        cohesive_force = 0.0;
    }



    /////////////////////////
    // DEM-FEM INTERACTION //  
    /////////////////////////

    void DEM_electromagnetic::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta){}

    void DEM_electromagnetic::CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                                         const double OldLocalContactForce[3],
                                                         double LocalElasticContactForce[3],
                                                         double LocalDeltDisp[3],
                                                         double LocalRelVel[3],
                                                         double indentation,
                                                         double previous_indentation,
                                                         double ViscoDampingLocalContactForce[3],
                                                         double& cohesive_force,
                                                         SphericParticle* const element,
                                                         Condition* const wall,
                                                         bool& sliding) 
    {
        //InitializeContactWithFEM(element, wall, indentation);
        const double distance = element->GetInteractionRadius() - indentation;

        const double smoother = 1.0;//std::max(1.0, 9.0 * indentation / (element1->GetInteractionRadius() + element2->GetInteractionRadius()));
        LocalElasticContactForce[2] = smoother * CalculateNormalForce(distance);
    }

    double DEM_electromagnetic::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation)
    {
        return 0.0;
    }

    double DEM_electromagnetic::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation)
    {
        return 0.0;
    }

    double DEM_electromagnetic::CalculateNormalForce(const double distance)
    {
          //const double F_elec = CalculateElectromagneticForce(distance);
          return 0.0;
          //return F_elec;
    }

    double DEM_electromagnetic::CalculateElectromagneticForce(const double distance)
    {
        return mCoulombParameterCharge / (distance * distance);
        //KRATOS_INFO("inside CalculateElectromagneticForces");
    }

} // namespace Kratos
