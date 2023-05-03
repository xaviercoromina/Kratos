import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseRoutine:
    def __init__(self, master_control: MasterControl, response: ResponseFunction) -> None:
        # set the master control
        self.__master_control = master_control

        # set the response
        self.__response = response
        self.__response_value = None

        self.__contributing_controls_list: 'list[Control]' = []
        self.__required_physical_gradients: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]' = {}

    def Initialize(self):
        """Initializes the response routine and underlying response.

        Raises:
            RuntimeError: If control domain and response domain does not have any intersection.
        """
        # initialize the response
        self.__response.Initialize()

        # create the required physical control fields to compute gradients
        self.__required_physical_gradients = self.__master_control.GetPhysicalKratosVariableCollectiveExpressionsMap()

        # now check which are the dependent physical space variables for the response, if not then remove
        # that variable
        list_of_independent_variables = []
        for required_physical_variable in self.__required_physical_gradients.keys():
            if required_physical_variable not in self.__response.GetDependentPhysicalKratosVariables():
                list_of_independent_variables.append(required_physical_variable)

        # now remove this independent collective expression from the require collective expressions map.
        for independent_variable in list_of_independent_variables:
            del self.__required_physical_gradients[independent_variable]

        for control in self.__master_control.GetListOfControls():
            # check whether control has keys given by required gradients
            if set(control.GetPhysicalKratosVariables()).intersection(self.__required_physical_gradients.keys()):
                # check whether there is an intersection of model parts between respones domain and control domain.
                #   1. in the case where response does not require an analysis, then intersection between evaluated and control domain is checked.
                #   2. in the case where response require an analysis, then intersection between analysis and control domain is checked.
                if self.__response.GetAnalysisModelPart() is None:
                    checked_model_part: Kratos.ModelPart = self.__response.GetEvaluatedModelPart()
                else:
                    checked_model_part: Kratos.ModelPart = self.__response.GetAnalysisModelPart()

                if Kratos.ModelPartOperationUtilities.HasIntersection([checked_model_part, control.GetEmptyField().GetModelPart()]):
                    self.__contributing_controls_list.append(control)

        if not self.__contributing_controls_list:
            raise RuntimeError(f"The controls does not have any influence over the response {self.GetReponseName()}.")

    def Check(self):
        self.__response.Check()

    def Finalize(self):
        self.__response.Finalize()

    def GetReponse(self) -> ResponseFunction:
        return self.__response

    def CalculateValue(self, control_field: KratosOA.ContainerExpression.CollectiveExpressions) -> float:
        """Calculates the value of the response.

        This method updates the design with the provided control field. If a control field is updated
        which affects the this response value, then a new value is computed. Otherwise, the design
        used in this reponse has not changed, hence the old value is returned without new calculations.

        Args:
            control_field (KratosOA.ContainerExpression.CollectiveExpressions): Control field of the new design.

        Returns:
            float: Respone value.
        """
        # update using the master control and get updated states.
        update_states = self.__master_control.Update(control_field)

        compute_response_value_flag = False
        for control, is_updated in update_states.items():
            if is_updated and control in self.__contributing_controls_list:
                compute_response_value_flag = True

        compute_response_value_flag = compute_response_value_flag or self.__response_value is None

        # TODO: In the case of having two analysis with the same mesh (model parts) for two different
        #       responses, we need to flag all the anayses which are affected by the control update_state
        #       from the first call, otherwise the second call will not update anything, hence no execution
        #       policies will be executed.

        #       Example: Drag 1 with 10m/s, Drag 2 with 20 m/s using the same mesh and model parts with two
        #                different analyses.
        # # now set the response value in response data container
        # if len(modified_model_parts) > 0:
        #     # if there are any modified model parts, then this changes state of is_executed to false
        #     # in the execution policies which has some intersection with the modified model parts.
        #     ChangeExecutionPolicyStates(modified_model_parts, False, self.__optimization_problem)

        if compute_response_value_flag:
            self.__response_value = self.__response.CalculateValue()

        return self.__response_value

    def CalculateGradient(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        """Returns Collective expression containing all the control space gradients for all control variable types (fields).

        Notes:
            1. It expects that the CalcaulteValue is called.
            2. The gradients are computed with respect to updates from master control.

        Returns:
            KratosOA.ContainerExpression.CollectiveExpressions: Returns mapped gradients collective expression.
        """
        # fills the proper physical gradients from the response
        self.__response.CalculateGradient(self.__required_physical_gradients)

        # calculate and return the control space gradients from respective controls
        return self.__master_control.MapGradient(self.__required_physical_gradients)
