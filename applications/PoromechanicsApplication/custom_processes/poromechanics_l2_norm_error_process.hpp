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

#include <fstream>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

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
        KRATOS_TRY;

        std::fstream l2_error_file;
        l2_error_file.open ("time_l2-rel-error_l2-abs-error.txt", std::fstream::out | std::fstream::app);
        l2_error_file.precision(12);
        l2_error_file << 0.0 << " " << 1.0 << " " << 1.0 << std::endl;
        l2_error_file.close();

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteAfterOutputStep() override
    {
        KRATOS_TRY;

        const ProcessInfo& r_current_process_info = mr_model_part.GetProcessInfo();
        const int NNodes = static_cast<int>(mr_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();


        double l2_abs_error = 1.0;
        double l2_rel_error = 1.0;

        double l2_numerator = 0.0;
        double l2_denominator = 0.0;

        #pragma omp parallel for reduction(+:l2_numerator,l2_denominator)
        for(int i = 0; i<NNodes; i++) {
            ModelPart::NodesContainerType::iterator itCurrentNode = it_begin + i;
            const array_1d<double, 3>& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& r_reference_displacement = itCurrentNode->FastGetSolutionStepValue(REFERENCE_DISPLACEMENT);
            const double& r_current_fluid_pressure = itCurrentNode->FastGetSolutionStepValue(WATER_PRESSURE);
            const double& r_reference_fluid_pressure = itCurrentNode->FastGetSolutionStepValue(REFERENCE_FLUID_PRESSURE);
            array_1d<double, 3> delta_displacement;
            noalias(delta_displacement) = r_current_displacement - r_reference_displacement;
            const double delta_fluid_pressure = r_current_fluid_pressure - r_reference_fluid_pressure;
            const double norm_2_du = inner_prod(delta_displacement,delta_displacement) + delta_fluid_pressure*delta_fluid_pressure;
            const double norm_2_u_ref = inner_prod(r_reference_displacement,r_reference_displacement) + r_reference_fluid_pressure*r_reference_fluid_pressure;

            l2_numerator += norm_2_du;
            l2_denominator += norm_2_u_ref;
        }
        if (l2_denominator > 1.0e-12) {
            l2_abs_error = std::sqrt(l2_numerator);
            l2_rel_error = l2_abs_error/std::sqrt(l2_denominator);
        }

        std::fstream l2_error_file;
        l2_error_file.open ("time_l2-rel-error_l2-abs-error.txt", std::fstream::out | std::fstream::app);
        l2_error_file.precision(12);
        l2_error_file << r_current_process_info[TIME] << " " << l2_rel_error << " " << l2_abs_error << std::endl;
        l2_error_file.close();

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
