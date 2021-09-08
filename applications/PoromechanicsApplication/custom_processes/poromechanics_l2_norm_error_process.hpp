//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_POROMECHANICS_L2_NORM_ERROR_PROCESS )
#define  KRATOS_POROMECHANICS_L2_NORM_ERROR_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"
#include "utilities/math_utils.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class PoromechanicsL2NormErrorProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsL2NormErrorProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    PoromechanicsL2NormErrorProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME"
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~PoromechanicsL2NormErrorProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the PoromechanicsL2NormErrorProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        // KRATOS_TRY;

        // int NCons = static_cast<int>(mr_model_part.Conditions().size());
        // ModelPart::ConditionsContainerType::iterator con_begin = mr_model_part.ConditionsBegin();

        // #pragma omp parallel for
        // for(int i = 0; i < NCons; i++)
        // {
        //     ModelPart::ConditionsContainerType::iterator itCond = con_begin + i;
        //     Condition::GeometryType& rGeom = itCond->GetGeometry();

        //     itCond->Set(PERIODIC,true);

        //     rGeom[0].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[1].Id();
        //     rGeom[1].FastGetSolutionStepValue(PERIODIC_PAIR_INDEX) = rGeom[0].Id();
        // }

        // int NElems = static_cast<int>(mr_model_part.Elements().size());
        // ModelPart::ElementsContainerType::iterator el_begin = mr_model_part.ElementsBegin();

        // #pragma omp parallel for
        // for(int i = 0; i < NElems; i++)
        // {
        //     ModelPart::ElementsContainerType::iterator itElem = el_begin + i;
        //     itElem->Set(ACTIVE,false);
        // }

        // KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        // TODO: calculate L2 Norm of Error and write it in a file

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PoromechanicsL2NormErrorProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PoromechanicsL2NormErrorProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    PoromechanicsL2NormErrorProcess& operator=(PoromechanicsL2NormErrorProcess const& rOther);

    /// Copy constructor.
    //PoromechanicsL2NormErrorProcess(PoromechanicsL2NormErrorProcess const& rOther);

}; // Class PoromechanicsL2NormErrorProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  PoromechanicsL2NormErrorProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PoromechanicsL2NormErrorProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_POROMECHANICS_L2_NORM_ERROR_PROCESS defined */
