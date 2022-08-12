# Importing the Kratos Library
from genericpath import isfile
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KM.Parameters(
        """{
            "help"                     : "This sets the initial conditions in terms of imposed strain, stress or deformation gradient",
            "model_part_name"          : "please_specify_model_part_name",
            "hole_generatrix_axis"     : [0.0,0.0,1.0],
            "hole_generatrix_point"    : [0.0,0.0,0.0],
            "hole_radius_offset"       : 0.0,
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
        }""")
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    layer_string = process_settings["initial_stress_table"]["filename"].GetString().split("_")[0]
    stress_string = process_settings["initial_stress_table"]["filename"].GetString().split("_")[1].split(".")[0]
      
    # for i in range(0,6):

    #     process_settings["initial_stress_table"]["filename"].SetString(layer_string + "_" + stress_string[:-1] + str(i+1)+".csv")
    #     process_settings["initial_stress_table"]["table_id"].SetInt(i)

    #     # print(process_settings["initial_stress_table"]["table_id"].SetInt(i))

    #     ReadCsvTableUtility(process_settings["initial_stress_table"]).Read(computing_model_part)
    #     print(computing_model_part.GetTable(i))
        
    #     if not isfile(layer_string + "_" + stress_string[:-1] + str(i+2) + ".csv"):
    #         break
    #     # i=+1
    
    # print(process_settings["model_part_name"].GetString())

    i=0

    while isfile(layer_string + "_" + stress_string[:-1] + str(i+1) + ".csv") and i<6:

        # print(layer_string + "_" + stress_string[:-1] + str(i+1) + ".csv")

        process_settings["initial_stress_table"]["filename"].SetString(layer_string + "_" + stress_string[:-1] + str(i+1)+".csv")
        process_settings["initial_stress_table"]["table_id"].SetInt(i)

        ReadCsvTableUtility(process_settings["initial_stress_table"]).Read(computing_model_part)
        # print(computing_model_part.GetTable(i))
        
        i += 1

        # if not isfile(layer_string[:-1] + str(i+1) + "_" + stress_string + ".csv"):
        #     Logger.PrintWarning("SetAutomatedInitialStressProcess:: ","Initial stress tables of " + layer_string[:-1] + str(i+1)+ " were not imported")
        #     break
            
    Logger.PrintInfo("SetAutomatedInitialStressProcess:: ","Initial stress tables of " + layer_string + " were successfully imported")

    return SMA.SetAutomatedInitialStressProcess(computing_model_part, process_settings)