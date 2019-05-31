//
// Author: Marc Chung To Sang mchungtosang@cimne.upc.edu
//

#if !defined(KRATOS_ION_PARTICLE_H_INCLUDED )
#define  KRATOS_ION_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "../DEMApplication/custom_elements/spheric_particle.h"


namespace Kratos
{
class KRATOS_API(PLASMA_DYNAMICS_APPLICATION) IonParticle : public SphericParticle
{
public:

    /// Pointer definition of IonParticle
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IonParticle);

    using SphericParticle::GetGeometry;
    using SphericParticle::GetDensity;
    using SphericParticle::mRadius;

    IonParticle():SphericParticle()
    {
    mSingleIonCharge = 1.60e-19; // in Coulomb, single charged ion, can go into node
    mDoubleIonCharge = 3.20e-19; // in Coulomb, Double charged ion
    mXenonMass = 2.18e-25; // in kg, Xenon is the most common gas used in plasma thrusters
    mExternalElectricField[0]=0.0 ; // External Electric Field initialized, should be improved to compute complexe external fields
    mExternalElectricField[1]=0.0 ; 
    mExternalElectricField[2]=0.0 ;
    mExternalMagneticField[0]=0.0 ; // External Magnetic Field initialized
    mExternalMagneticField[1]=0.0 ;
    mExternalMagneticField[2]=0.0 ;
    }

    IonParticle( IndexType NewId, GeometryType::Pointer pGeometry ):SphericParticle(NewId, pGeometry)
    {   
                        mSingleIonCharge = 1.60e-19;
                        mDoubleIonCharge = 3.20e-19;
                        mXenonMass = 2.18e-25;
                        mExternalElectricField[0]=0.0 ; // mExternalElectricField[2]= {0.0};
                        mExternalElectricField[1]=0.0 ; 
                        mExternalElectricField[2]=0.0 ;
                        mExternalMagneticField[0]=0.0 ;
                        mExternalMagneticField[1]=0.0 ;
                        mExternalMagneticField[2]=0.0 ;                       
    }
    IonParticle( IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes)
    {
                        mSingleIonCharge = 1.60e-19;
                        mDoubleIonCharge = 3.20e-19;
                        mXenonMass = 2.18e-25;
                        mExternalElectricField[0]=0.0 ;
                        mExternalElectricField[1]=0.0 ;
                        mExternalElectricField[2]=0.0 ;
                        mExternalMagneticField[0]=0.0 ;
                        mExternalMagneticField[1]=0.0 ;
                        mExternalMagneticField[2]=0.0 ;
    }
    IonParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):SphericParticle(NewId, pGeometry, pProperties)
    {
                        mSingleIonCharge = 1.60e-19;
                        mDoubleIonCharge = 3.20e-19;
                        mXenonMass = 2.18e-25;
                        mExternalElectricField[0]=0.0 ;
                        mExternalElectricField[1]=0.0 ;
                        mExternalElectricField[2]=0.0 ;
                        mExternalMagneticField[0]=0.0 ;
                        mExternalMagneticField[1]=0.0 ;
                        mExternalMagneticField[2]=0.0 ;
    }

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return SphericParticle::Pointer(new IonParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor.
    virtual ~IonParticle();


    std::vector<Node<3>::Pointer> mNeighbourNodes;
    std::vector<double> mNeighbourNodesDistances;

    
    /// Assignment operator.
    IonParticle& operator=(IonParticle const& rOther); 

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IonParticle" ;
        return buffer.str();
    }
    
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "IonParticle";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}    

    void Initialize(const ProcessInfo& r_process_info) override;

    void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                 array_1d<double, 3>& additionally_applied_moment,
                                 const ProcessInfo& r_current_process_info,
                                 const array_1d<double,3>& gravity) override;

    void MemberDeclarationFirstStep(const ProcessInfo& r_process_info) override;
    
    virtual void CalculateCoulombForce(array_1d<double, 3>& Coulomb_force);
    virtual void CalculateLaplaceForce(array_1d<double, 3>& Laplace_force);


    double GetSingleIonCharge();
    double GetDoubleIonCharge();
    double GetXenonMass();
    array_1d<double, 3> GetExternalElectricField();
    array_1d<double, 3> GetExternalMagneticField();


protected:

    double mSingleIonCharge;
    double mDoubleIonCharge;
    double mXenonMass;
    array_1d<double, 3> mExternalElectricField;
    array_1d<double, 3> mExternalMagneticField;


private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
        rSerializer.save("mSingleIonCharge",mSingleIonCharge);
        rSerializer.save("mDoubleIonCharge",mDoubleIonCharge);
        rSerializer.save("mXenonMass",mXenonMass);
        rSerializer.save("mExternalElectricField",mExternalElectricField);
        rSerializer.save("mExternalMagneticField",mExternalMagneticField);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
        rSerializer.load("mSingleIonCharge",mSingleIonCharge);
        rSerializer.load("mDoubleIonCharge",mDoubleIonCharge);
        rSerializer.load("mXenonMass",mXenonMass);  
        rSerializer.load("mExternalElectricField",mExternalElectricField);
        rSerializer.load("mExternalMagneticField",mExternalMagneticField);      
    }    

    /// Copy constructor.
    IonParticle(IonParticle const& rOther)
    {
    *this = rOther;
    }


}; // Class IonParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, IonParticle& rThis) {return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const IonParticle& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_ION_PARTICLE_H_INCLUDED  defined
