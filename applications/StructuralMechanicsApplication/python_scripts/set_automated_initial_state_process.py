# Importing the Kratos Library
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
        }""")
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]
    process_settings.RemoveValue("model_part_name")

    initial_state_table = ReadCsvTableUtility(process_settings["initial_state_table"]).Read(computing_model_part)
    Logger.PrintInfo("SetAutomatedInitialStateProcess:: ","Initial state table was sucessfully imported")
    #print(computing_model_part.GetTable(0))


    #print("How much is the value interpolated of 3.77? --> ", initial_state_table.GetValue(3.77))

    return SMA.SetAutomatedInitialStateProcess(computing_model_part, process_settings)