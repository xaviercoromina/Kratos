import KratosMultiphysics
from importlib import import_module

def CreateSolverByParameters(model, solver_settings, parallelism):

    solver_type_raw = solver_settings["solver_type"].GetString()
    
    if solver_settings.Has("time_integration_method"):
        time_integration_method = solver_settings["time_integration_method"].GetString()
    else:
        time_integration_method = "implicit" # defaulting to implicit time-integration
        
    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        solver_type = solver_type_raw.lower()
        if solver_type in ("u_pw", "geomechanics_u_pw_solver", "twophase"):
            #solver_settings["solver_settings"]["time_stepping"].AddValue("end_time", solver_settings["problem_data"]["end_time"])
            solver_module_name = "geomechanics_U_Pw_solver"

        elif solver_type in ("pw", "geomechanics_pw_solver"):
            #solver_settings["solver_settings"]["time_stepping"].AddValue("end_time", solver_settings["problem_data"]["end_time"])
            solver_module_name = "geomechanics_Pw_solver"

        elif (solver_type.lower() == "t" or solver_type.lower() == "geomechanics_t_solver" or solver_type.lower() == "twophase"):
            #solver_settings["solver_settings"]["time_stepping"].AddValue("end_time", solver_settings["problem_data"]["end_time"])
            solver_module_name = "geomechanics_T_solver"
            
        elif (solver_type.lower() == "tpw" or solver_type.lower() == "thermalpressurecoupled" or solver_type.lower() == "twophase"):
            #solver_settings["solver_settings"]["time_stepping"].AddValue("end_time", solver_settings["problem_data"]["end_time"])
            solver_module_name = "coupled_thermal_pressure_solver"

        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"geomechanics_U_Pw_solver\""
            raise Exception(err_msg)

    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\""
        raise Exception(err_msg)

    module_full_name = 'KratosMultiphysics.GeoMechanicsApplication.' + solver_module_name
    solver = import_module(module_full_name).CreateSolver(model, solver_settings)

    return solver
    
    
def CreateSolver(model, custom_settings):

    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("input is expected to be provided as a Kratos Model object")#

    if not isinstance(custom_settings, KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    solver_settings = custom_settings["solver_settings"]
    parallelism = custom_settings["problem_data"]["parallel_type"].GetString()

    return CreateSolverByParameters(model, solver_settings, parallelism)