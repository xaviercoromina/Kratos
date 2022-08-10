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

#include "includes/model_part.h"
#include "custom_processes/set_automated_initial_state_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetAutomatedInitialStateProcess::SetAutomatedInitialStateProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialStateProcess::ExecuteInitialize()
{
    KRATOS_TRY
    array_1d<double, 6> strain_vector;
    array_1d<double, 6> initial_stress_vector;
    array_1d<double, 3> element_centroid;
    double element_radial_cordinate;

    // KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        
        int TableId = mThisParameters["table_id"].GetInt();
        element_centroid = rElement.GetGeometry().Center();
        element_radial_cordinate = sqrt(element_centroid[0] * element_centroid[0] + element_centroid[1] * element_centroid[1]);
        initial_stress_vector[0] = mrThisModelPart.GetTable(TableId).GetValue(element_radial_cordinate);
        initial_stress_vector[1] = 0;
        initial_stress_vector[2] = 0;
        initial_stress_vector[3] = 0;
        initial_stress_vector[4] = 0;
        initial_stress_vector[5] = 0;

        rElement.SetValue(INITIAL_STRESS_VECTOR, initial_stress_vector );
        
    });


    KRATOS_CATCH("")
}

} // namespace Kratos.