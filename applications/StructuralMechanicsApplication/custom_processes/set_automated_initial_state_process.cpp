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

    // const array_1d<double, 3> element_centroid;
    //  element_radial_cordinate;

    // KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;

    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {
        
        
        int TableId;
        int ElemId;
        array_1d<double, 3> element_centroid;
        double element_radial_cordinate;
        array_1d<double, 6> initial_stress_vector;

        TableId = mThisParameters["initial_state_table"]["table_id"].GetInt();
        KRATOS_WATCH(TableId)
        
        ElemId = rElement.Id();
        KRATOS_WATCH(ElemId)    
        
        element_centroid = rElement.GetGeometry().Center();
        KRATOS_WATCH(element_centroid)
        
        element_radial_cordinate = sqrt(element_centroid[0] * element_centroid[0] + element_centroid[1] * element_centroid[1]);
        KRATOS_WATCH(element_radial_cordinate)    

        initial_stress_vector[0] = mrThisModelPart.GetTable(TableId).GetValue(element_radial_cordinate);
        initial_stress_vector[1] = 0;
        initial_stress_vector[2] = 0;
        initial_stress_vector[3] = 0;
        initial_stress_vector[4] = 0;
        initial_stress_vector[5] = 0;
        KRATOS_WATCH(initial_stress_vector)

        rElement.SetValue(INITIAL_STRESS_VECTOR, initial_stress_vector);
        KRATOS_WATCH(rElement.GetValue(INITIAL_STRESS_VECTOR))
        
    });
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetAutomatedInitialStateProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "help"                     : "This sets the initial conditions in terms of imposed strain, stress or deformation gradient",
        "mesh_id"                  : 0,
        "model_part_name"          : "please_specify_model_part_name",
        "dimension"                : 3,
        "initial_state_table"      : {
                    "name"             : "csv_table",
                    "filename"         : "sample.csv",
                    "delimiter"        : ",",
                    "skiprows"         : 1,
                    "first_column_id"  : 0,
                    "second_column_id" : 1,
                    "table_id"         : -1,
                    "na_replace"       : 0.0
                }
    })");

    return default_parameters;
}

} // namespace Kratos.