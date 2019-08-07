from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.PlasmaDynamicsApplication import *

from sphere_strategy import ExplicitStrategy
BaseStrategy = ExplicitStrategy
class PlasmaStrategy(BaseStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, parameters, procedures):
        self.project_parameters = parameters
        super(PlasmaStrategy, self).__init__(all_model_parts, creator_destructor, dem_fem_search, parameters['dem_parameters'], procedures)

    def TranslationalIntegrationSchemeTranslator(self, name):
        class_name = BaseStrategy.TranslationalIntegrationSchemeTranslator(self, name)

        return class_name

    def RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational):
        class_name = BaseStrategy.RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational)

        return class_name

    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()

        if (self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, self.solver_settings)
        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.search_strategy, self.solver_settings)

    def GetTranslationalSchemeInstance(self, class_name):
        try:
            translational_scheme = super(PlasmaStrategy, self).GetTranslationalSchemeInstance(class_name)
        except Exception:
            translational_scheme = PlasmaStrategy.class_name()
        return translational_scheme

    def GetRotationalSchemeInstance(self, class_name):
        try:
            rotational_scheme = super(PlasmaStrategy, self).GetRotationalSchemeInstance(class_name)
        except Exception:
            rotational_scheme = PlasmaStrategy.class_name()
        return rotational_scheme

    def GetPlasmaConstitutiveLawParametersIfItExists(self, properties):
        if self.project_parameters.Has('properties'):
            for p in self.project_parameters["properties"]:
                if p['properties_id'].GetInt() == int(properties.Id) and p.Has('plasma_dynamics_law_parameters'):
                    return p['plasma_dynamics_law_parameters']
        return None

    @staticmethod
    def CreatePlasmaDynamicsLaw(properties, plasma_dynamics_law_parameters):

        #plasma_dynamics_name = plasma_dynamics_law_parameters['name'].GetString() 
        # TODO: remove this when decided on which way to call the name of the law
        #PlasmaDynamicsInteractionLaw = globals().get(properties['DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME'])()
        #PlasmaDynamicsInteractionLaw.SetConstitutiveLawInProperties(properties, True)
        pass

    def ModifyProperties(self, properties, param = 0):
        #print(properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME])
        #globals()['DEM_electromagnetic'] = DEM_electromagnetic()
        super(PlasmaStrategy,self).ModifyProperties(properties, param)

        #plasma_dynamics_law_parameters = self.GetPlasmaConstitutiveLawParametersIfItExists(properties)
        #if plasma_dynamics_law_parameters:
        #    PlasmaStrategy.CreatePlasmaDynamicsLaw(properties, plasma_dynamics_law_parameters)

        if not param:
            if not properties.Has(PARTICLE_SPHERICITY):
                properties[PARTICLE_SPHERICITY] = 1.0