import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion

class ContainerVariableDataHolderOutputDecorator(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optmization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "optimization_info_data_location": "PLEASE/PROVIDE/DATA/LOCATION/SEPERATED/BY/BACK/SLASH",
            "output_variable_name"           : "PLEASE_PROVIDE_MATCHING_VARIABLE_NAME_FOR_DATA",
            "output_processes"               : []
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__optimization_info = optmization_info
        self.__optimization_info_location = parameters["optimization_info_data_location"].GetString()
        self.__output_variable = Kratos.KratosGlobals.GetVariable(parameters["output_variable_name"].GetString())
        self.__output_processes_list: 'list[Kratos.OutputProcess]' = KratosProcessFactory(model).ConstructListOfProcesses(parameters["output_processes"])

        for output_process in self.__output_processes_list:
            if not isinstance(output_process, Kratos.OutputProcess):
                raise RuntimeError(f"The provided output process is not of the type Kratos.OutputProcess. Details of the process: " + str(output_process))

    def ExecuteInitialize(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteInitialize)

    def ExecuteBeforeSolutionLoop(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteBeforeSolutionLoop)

    def ExecuteInitializeSolutionStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteInitializeSolutionStep)

    def ExecuteBeforeOutputStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteBeforeOutputStep)

    def ExecuteAfterOutputStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteAfterOutputStep)

    def ExecuteFinalizeSolutionStep(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteFinalizeSolutionStep)

    def ExecuteFinalize(self):
        CallOnAll(self.__output_processes_list, Kratos.OutputProcess.ExecuteFinalize)

    def IsOutputStep(self) -> bool:
        is_output_step = False
        for output_process in self.__output_processes_list:
            is_output_step = is_output_step | output_process.IsOutputStep()
        return is_output_step

    def PrintOutput(self):
        # obtain the container variable data holder
        current_pos = self.__optimization_info.GetSolutionStepData(0)
        for location_sub_key in self.__optimization_info_location.split("/"):
            try:
                current_pos = current_pos[location_sub_key]
            except KeyError:
                raise KeyError(f"The requested subkey \"{location_sub_key}\" for location \"{self.__optimization_info_location}\" is not found. Followings are available subkeys: " + "\n\t".join(current_pos.keys()))
            except Exception:
                raise Exception(f"The requested subkey \"{location_sub_key}\" for location \"{self.__optimization_info_location}\" is not found. Please check the location if it is correct.")

        current_pos: ContainerVariableDataHolderUnion = current_pos
        if not isinstance(current_pos, ContainerVariableDataHolderUnion):
            raise TypeError(f"The provided location at \"{self.__optimization_info_location}\" does not contain a container variable data holder.")

        if current_pos.GetDataDimension() != KratosOA.OptimizationUtils.GetVariableDimension(self.__output_variable):
            raise RuntimeError(f"Provided {self.__output_variable.Name()} is in compatible with the data container dimension of {current_pos.GetDataDimension()} for location at \"{self.__optimization_info_location}\".")

        temp_values = current_pos.Clone()
        temp_values.ReadDataFromContainerVariable(self.__output_variable)

        current_pos.AssignDataToContainerVariable(self.__output_variable)
        for output_process in self.__output_processes_list:
            if output_process.IsOutputStep():
                output_process.PrintOutput()

        temp_values.AssignDataToContainerVariable(self.__output_variable)






