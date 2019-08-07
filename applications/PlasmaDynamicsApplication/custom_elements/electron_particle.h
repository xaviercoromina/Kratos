//
// Author: Marc Chung To Sang mchungtosang@cimne.upc.edu
//

#if !defined(KRATOS_ELECTRON_PARTICLE_H_INCLUDED )
#define  KRATOS_ELECTRON_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "ion_particle.h"


namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) ElectronParticle : public IonParticle
{
public:

    /// Pointer definition of ElectronParticle
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ElectronParticle);


    ElectronParticle():IonParticle()
    {
    mElectronCharge = -1.60e-19; // in Coulomb, Hard-coded but should go into node
    }

    ElectronParticle( IndexType NewId, GeometryType::Pointer pGeometry ):IonParticle(NewId, pGeometry)
    {   
                        mElectronCharge = -1.60e-19;                     
    }
    ElectronParticle( IndexType NewId, NodesArrayType const& ThisNodes):IonParticle(NewId, ThisNodes)
    {
                        mElectronCharge = -1.60e-19;
    }
    ElectronParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):IonParticle(NewId, pGeometry, pProperties)
    {
                        mElectronCharge = -1.60e-19;
    }

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return SphericParticle::Pointer(new ElectronParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor.
    virtual ~ElectronParticle();

    /// Assignment operator.
    ElectronParticle& operator=(const ElectronParticle& rOther); 

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ElectronParticle" ;
        return buffer.str();
    }
    
    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ElectronParticle";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}    

    void Initialize(const ProcessInfo& r_process_info) override;

    void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;

    void MemberDeclarationFirstStep(const ProcessInfo& r_process_info) override;
    
    void CalculateCoulombForce(array_1d<double, 3>& Coulomb_force) override;
    void CalculateLaplaceForce(array_1d<double, 3>& Laplace_force) override;



protected:

    double mElectronCharge;



private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
        rSerializer.save("mElectronCharge",mElectronCharge);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
        rSerializer.load("mElectronCharge",mElectronCharge);     
    }    

    /// Copy constructor.
    ElectronParticle(ElectronParticle const& rOther)
    {
    *this = rOther;
    }


}; // Class ElectronParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, ElectronParticle& rThis) {return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const ElectronParticle& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_ELECTRON_PARTICLE_H_INCLUDED  defined
