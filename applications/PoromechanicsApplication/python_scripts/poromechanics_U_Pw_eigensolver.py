# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.PoromechanicsApplication.poromechanics_U_Pw_solver import UPwSolver

from KratosMultiphysics import eigen_solver_factory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

def CreateSolver(model, custom_settings):
    return UPwEigenSolver(model, custom_settings)

class UPwEigenSolver(UPwSolver):
    """The Poromechanics eigen solver.

    This class creates the solvers for eigenvalue analysis.
    """
    def __init__(self, model, custom_settings):
        if custom_settings.Has("linear_solver_settings"):
            IssueDeprecationWarning('UPwEigenSolver', '"linear_solver_settings" was specified which is not used in the UPwEigenSolver. Use "eigensolver_settings"!')
            custom_settings.RemoveValue("linear_solver_settings")

        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[UPwEigenSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "compute_modal_decomposition": false,
            "eigensolver_settings" : {
                "solver_type"           : "spectra_sym_g_eigs_shift",
                "max_iteration"         : 1000,
                "number_of_eigenvalues" : 5,
                "echo_level"            : 1
            },
            "eigensolver_diagonal_values" : { }
        }""")
        base_parameters = super().GetDefaultParameters()
        base_parameters.RemoveValue("linear_solver_settings")
        this_defaults.AddMissingParameters(base_parameters)
        return this_defaults

    #### Private functions ####

    def _ConstructScheme(self, scheme_type, solution_type):

        scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()

        return scheme

    def _ConstructLinearSolver(self):
        """Create the eigensolver.

        This overrides the base class method and replaces the usual linear solver
        with an eigenvalue problem solver.
        """
        return eigen_solver_factory.ConstructSolver(self.settings["eigensolver_settings"])

    def _ConstructSolver(self, builder_and_solver, strategy_type):

        self.linear_solver = _ConstructLinearSolver()
        eigen_scheme = self.scheme # The scheme defines the matrices of the eigenvalue problem.
        builder_and_solver = self._ConstructBuilderAndSolver(True) # The eigensolver is created here.

        solver_type = self.settings["eigensolver_settings"]["solver_type"].GetString()
        if solver_type in ["eigen_eigensystem", "spectra_sym_g_eigs_shift"]: # TODO evaluate what has to be used for spectra
            mass_matrix_diagonal_value = 0.0
            stiffness_matrix_diagonal_value = 1.0
        elif solver_type == "feast":
            mass_matrix_diagonal_value = 1.0
            stiffness_matrix_diagonal_value = -1.0
        else:
            diag_values = self.settings["eigensolver_diagonal_values"]
            if not diag_values.Has("mass_matrix_diagonal_value") or not diag_values.Has("stiffness_matrix_diagonal_value"):
                err_msg  = 'For the used eigensolver "{}" no defaults for '.format(solver_type)
                err_msg += '"mass_matrix_diagonal_value" and "stiffness_matrix_diagonal_value" exist, '
                err_msg += 'please specify them under "eigensolver_diagonal_values"'
                raise Exception(err_msg)

            mass_matrix_diagonal_value = diag_values["mass_matrix_diagonal_value"].GetDouble()
            stiffness_matrix_diagonal_value = diag_values["stiffness_matrix_diagonal_value"].GetDouble()

        return StructuralMechanicsApplication.EigensolverStrategy(self.computing_model_part,
                                                                  self.scheme,
                                                                  builder_and_solver,
                                                                  mass_matrix_diagonal_value,
                                                                  stiffness_matrix_diagonal_value,
                                                                  self.settings["compute_modal_decomposition"].GetBool())
