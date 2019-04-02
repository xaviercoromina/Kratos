// Authors: Miguel Angel Celigueta, maceli@cimne.upc.edu
//          Guillermo Casas, gcasas@cimne.upc.edu
//          Marc Chung To Sang, mchungtosang@cimne.upc.edu

#if !defined(DEM_ELECTROMAGNETIC_CL_H_INCLUDED)
#define  DEM_ELECTROMAGNETIC_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "../DEMApplication/custom_constitutive/DEM_discontinuum_constitutive_law.h"


namespace Kratos {
    
    class SphericParticle;

    class KRATOS_API(PLASMA_DYNAMICS_APPLICATION) DEM_electromagnetic : public DEMDiscontinuumConstitutiveLaw 
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_electromagnetic);

        DEM_electromagnetic()
        {
            mCoulombParameterCharge = 2.3e-28; //Charge times charge times 1/(4 * Pi * ElectricPermittivity) in N.m2
        } 

        void Initialize(const ProcessInfo& r_process_info) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;
        
        std::string GetTypeOfLaw() override;

        ~DEM_electromagnetic() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;  

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

        void CalculateForces(const ProcessInfo& r_process_info,
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
                            bool& sliding, double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(ProcessInfo& r_process_info,
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
                                    bool& sliding) override;

        double CalculateCohesiveNormalForce(SphericParticle* const element1,
                                            SphericParticle* const element2,
                                            const double indentation) override;

        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                            Condition* const wall,
                                            const double indentation) override;

        double CalculateNormalForce(const double distance) override;

        double CalculateElectromagneticForce(const double distance);
        
        double mCoulombParameterCharge;
        
    private:


        
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override 
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override 
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} /* namespace Kratos.*/
#endif /* DEM_electromagnetic_CL_H_INCLUDED  defined */
