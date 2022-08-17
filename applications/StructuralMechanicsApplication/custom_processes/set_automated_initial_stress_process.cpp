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
#include "custom_processes/set_automated_initial_stress_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
SetAutomatedInitialStressProcess::SetAutomatedInitialStressProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
}

/***********************************************************************************/
/***********************************************************************************/

void SetAutomatedInitialStressProcess::ExecuteInitialize()
{
    // KRATOS_WATCH(mrThisModelPart.Tables())
    KRATOS_TRY
        
    // KRATOS_ERROR_IF(MathUtils<double>::Norm3(r_generatrix_axis) < std::numeric_limits<double>::epsilon()) << "The r_generatrix_axis has norm zero" << std::endl;
    const double tol = 1.0e-6;
    
    block_for_each(mrThisModelPart.Elements(), [&](Element &rElement) {

        int TableFirstId;       
        int TableId;
        int ElemId;
        int index;
        array_1d<double, 3> hole_generatrix_axis;   
        array_1d<double, 3> normalized_generatrix_vector;   
        array_1d<double, 3> element_centroid;
        array_1d<double, 3> hole_generatrix_point;
        array_1d<double, 3> relative_position_vector;
        double vector_scaler;
        array_1d<double, 3> intersection_point;
        double hole_radius_offset;
        array_1d<double, 3> radial_position_vector;
        double centroid_relative_distance;
        array_1d<double, 6> initial_stress_vector;
        double radial_vector_norm;

        ElemId = rElement.Id();
        // KRATOS_WATCH(ElemId)    
        
        noalias(hole_generatrix_axis) = mThisParameters["hole_generatrix_axis"].GetVector();
        // KRATOS_WATCH(hole_generatrix_axis)  

        noalias(hole_generatrix_point) = mThisParameters["hole_generatrix_point"].GetVector();

        normalized_generatrix_vector[0] = (hole_generatrix_axis[0] * hole_generatrix_axis[0]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);
        normalized_generatrix_vector[1] = (hole_generatrix_axis[1] * hole_generatrix_axis[1]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);
        normalized_generatrix_vector[2] = (hole_generatrix_axis[2] * hole_generatrix_axis[2]) / std::sqrt(hole_generatrix_axis[0] * hole_generatrix_axis[0] + hole_generatrix_axis[1] * hole_generatrix_axis[1] + hole_generatrix_axis[2] * hole_generatrix_axis[2]);
        // KRATOS_WATCH(normalized_generatrix_vector);

        noalias(element_centroid) = rElement.GetGeometry().Center();
        // KRATOS_WATCH(element_centroid)

        relative_position_vector[0] = element_centroid[0] - hole_generatrix_point[0];
        relative_position_vector[1] = element_centroid[1] - hole_generatrix_point[1];
        relative_position_vector[2] = element_centroid[2] - hole_generatrix_point[2];
        // KRATOS_WATCH(relative_position_vector);

        vector_scaler = relative_position_vector[0] * normalized_generatrix_vector[0] + relative_position_vector[1] * normalized_generatrix_vector[1] + relative_position_vector[2] * normalized_generatrix_vector[2];
        // KRATOS_WATCH(vector_scaler)
        
        noalias(intersection_point) = hole_generatrix_point + vector_scaler * normalized_generatrix_vector;
        // KRATOS_WATCH(intersection_point)

        hole_radius_offset = mThisParameters["hole_radius_offset"].GetDouble();
        // KRATOS_WATCH(hole_radius_offset)  

        noalias(radial_position_vector) = element_centroid - intersection_point;
        // KRATOS_WATCH(radial_position_vector)

        // radial_vector_norm = std::sqrt(radial_position_vector[0] * radial_position_vector[0] + radial_position_vector[1] * radial_position_vector[1] + radial_position_vector[2] * radial_position_vector[2]);
        // KRATOS_WATCH(radial_vector_norm)

        centroid_relative_distance = std::sqrt(radial_position_vector[0] * radial_position_vector[0] + radial_position_vector[1] * radial_position_vector[1] + radial_position_vector[2] * radial_position_vector[2]) - hole_radius_offset;
        // KRATOS_ERROR_IF(centroid_relative_distance < 0.0) << "Thickness of element " << ElemId << " is too small." << std::endl;
        
        if (centroid_relative_distance < 0.0 && abs(centroid_relative_distance) <= tol) {
            centroid_relative_distance = 0.0;
        }
        else if(centroid_relative_distance < 0.0 && abs(centroid_relative_distance) >= tol){
            KRATOS_ERROR << "Thickness of element " << ElemId << " is too small." << std::endl;
        }
            
        
        // KRATOS_WATCH(centroid_relative_distance)    

        TableFirstId = mThisParameters["initial_stress_table"]["table_id"].GetInt()-mrThisModelPart.Tables().size()+1;
        
        // KRATOS_WATCH(TableFirstId)

        index=0;    
        for (IndexType TableId = TableFirstId; TableId < TableFirstId + mrThisModelPart.Tables().size(); ++TableId) {
                
                // KRATOS_WATCH(TableId)

                initial_stress_vector[index] = mrThisModelPart.GetTable(TableId).GetValue(centroid_relative_distance);

                index+=1;
        }

        rElement.SetValue(INITIAL_STRESS_VECTOR, initial_stress_vector);
        // KRATOS_WATCH(rElement.GetValue(INITIAL_STRESS_VECTOR))
    });
    
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SetAutomatedInitialStressProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "help"                     : "This automates the application of initial conditions in terms of imposed stress",
        "model_part_name"          : "please_specify_model_part_name",
        "hole_generatrix_axis"     : [0.0,0.0,1.0],
        "hole_generatrix_point"    : [0.0,0.0,0.0],
        "hole_radius_offset"       : 0.0,
        "last_layer"               : false,
        "initial_stress_table"     : {
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