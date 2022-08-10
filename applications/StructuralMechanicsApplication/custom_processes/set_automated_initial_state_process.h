// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#if !defined(KRATOS_SET_AUTOMATED_INITIAL_STATE_PROCESS )
#define  KRATOS_SET_AUTOMATED_INITIAL_STATE_PROCESS

#include "processes/process.h"

namespace Kratos
{

/**
 * @class SetAutomatedInitialStateProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class set the local axes of the elements according to a cylindrical coordinates
 * @author Alejandro Cornejo
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SetAutomatedInitialStateProcess
    : public Process
{

public:


    KRATOS_CLASS_POINTER_DEFINITION(SetAutomatedInitialStateProcess);


    /// Constructor
    SetAutomatedInitialStateProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );


    /// Destructor
    ~SetAutomatedInitialStateProcess() override = default;

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function is designed for being called at the beginning each time step
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SetAutomatedInitialStateProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetAutomatedInitialStateProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


protected:

    /// Member Variables

    ModelPart& mrThisModelPart;
    Parameters mThisParameters;

private:

    /// Assignment operator.
    SetAutomatedInitialStateProcess& operator=(SetAutomatedInitialStateProcess const& rOther);

    /// Copy constructor.
    //SetAutomatedInitialStateProcess(SetAutomatedInitialStateProcess const& rOther);

}; // Class SetAutomatedInitialStateProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetAutomatedInitialStateProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetAutomatedInitialStateProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_SET_CYLINDRICAL_LOCAL_AXES_PROCESS defined */
