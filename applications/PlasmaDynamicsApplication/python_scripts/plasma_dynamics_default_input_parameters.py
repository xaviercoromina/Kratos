from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics

def GetDefaultInputParameters():

    default_settings = KratosMultiphysics.Parameters(
        """{
            "Dimension" : 3,
            "GravityX" : 0.0,
            "GravityY" : 0.0,
            "GravityZ" : 0.0,
            "RotationOption" : false,

            "strategy_parameters"            : {
                "strategy"                   : "plasma_sphere_strategy",
                "RemoveBallsInitiallyTouchingWalls" : false
            },
            "do_print_results_option" : true,
            "full_particle_history_watcher" : "Empty",

            "OutputFileType" : "Binary",
            "Multifile" : "multiple_files",

            "TranslationalIntegrationScheme" : "Symplectic_Euler",
            "RotationalIntegrationScheme"    : "Direct_Integration",
            "MaxTimeStep" : 0.005,
            "FinalTime" : 1.0,
            "ControlTime" : 4.0,
            "OutputTimeStep" : 0.5,

            "coupling_level_type" : 0,
            "NeighbourSearchFrequency" : 1,
            "interaction_start_time" : 0.0,
            "do_search_neighbours" : false,
            "do_solve_dem" : true,

            "TestType" : "None",

            "ElementType" : "IonParticle3D",
            "echo_level" : 1,
            "problem_data" : {
                "problem_name" : "dummy_name2.Provide_a_real_one",
                "parallel_type" : "OpenMP",
                "echo_level" : 1,
                "start_time" : 0.0,
                "end_time" : 1
            },


            "processes" : {
                "auxiliar_process_list": []
            },
            "json_output_process" : [],
            "plasma_dynamics_output_processes" : {},
            "fluid_already_calculated" : false,

            "store_full_gradient_option" : false,

            "coupling_weighing_type" : 2,
            "coupling_weighing_type_comment" : "{fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)",
            "fluid_model_type" : 1,
            "fluid_model_type_comment" : " untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)",
            "print_particles_results_option" : false,

            "stationary_problem_option" : false,
            "stationary_problem_option_comment" : " stationary, stop calculating the fluid after it reaches the stationary state",


            "embedded_option" : false,
            "embedded_option_comment" : "the embedded domain tools are to be used",
            "make_results_directories_option" : true,
            "make_results_directories_option_comment": "results are written into a folder (../results) inside the problem folder",

            "print_debug_info_option" : false,
            "print_debug_info_option_comment" : " print a summary of global physical measures",
            "print_particles_results_cycle" : 1,
            "print_particles_results_cycle_comment" : " number of 'ticks' per printing cycle",
            "debug_tool_cycle" : 10,
            "debug_tool_cycle_comment" : " number of 'ticks' per debug computations cycle",
            "similarity_transformation_type" : 0,
            "similarity_transformation_type_comment" : " no transformation (0), Tsuji (1)",

            "min_fluid_fraction" : 0.2,

            "model_over_real_diameter_factor" : 1.0,
            "model_over_real_diameter_factor_comment": " not active if similarity_transformation_type = 0",


            "time_steps_per_analytic_processing_step": 1,
            "do_process_analytic_data" : true,
            "meso_scale_length" : 0.2,
            "meso_scale_length_comment" : " the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)",
            "shape_factor" : 0.5,


            "material_acceleration_calculation_type" : 0,
            "laplacian_calculation_type" : 0,
            "vorticity_calculation_type" : 0,
            "pressure_grad_recovery_type" : 0,

            "basset_force_type" : 0,
            "basset_force_integration_type" : 2,
            "n_init_basset_steps" : 0,
            "time_window" : 0.04,
            "number_of_exponentials" : 2,
            "frame_of_reference_type" : 0,
            "angular_velocity_of_frame_Z" : 0.0,

            "filter_velocity_option" : false,
            "apply_time_filter_to_fluid_fraction_option" : false,

            "type_of_dem_inlet" : "VelocityImposed",
            "type_of_dem_inlet_comment" : "VelocityImposed or ForceImposed",



            "dem_nodal_results" : {
                                   "SLIP_VELOCITY" : false,
                                   "RADIUS" : false,
                                   "ANGULAR_VELOCITY" : false,
                                   "ELASTIC_FORCES" : false,
                                   "CONTACT_FORCES" : false,
                                   
                                   "EXTERNAL_APPLIED_FORCE" : false,
                                   
                                   "FLUID_VEL_LAPL_PROJECTED" : false,
                                   
                                   "FLUID_ACCEL_PROJECTED" : false,
                                   "FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED" : false,
                                   
                                   
                                   "FLUID_FRACTION_PROJECTED" : false,

                                   "BASSET_FORCE" : false,

                                   "PRESSURE" : false,
                                   
                                   "ELECTRIC_FIELD_PROJECTED_TO_PARTICLE": false,
                                   "NUMBER_OF_PARTICLES_IN_A_ION_MACROPARTICLE": false,
                                   "PARTICLE_ION_VELOCITY" : false,
                                   "MACROPARTICLE_ION_DENSITY": false,       
                                   "FLUID_VEL_PROJECTED" : false},

            "fluid_nodal_results" : {
                                     "VELOCITY_GRADIENT" : false,
                                     
                                     "AVERAGED_FLUID_VELOCITY" : false,
                                     "FLUID_FRACTION" : false,
                                     "FLUID_FRACTION_OLD" : false,
                                     "PARTICLE_VEL_FILTERED" : false,
                                
                                     "FLUID_FRACTION_GRADIENT" : false,

                                     "DISTANCE" : false,
                                     "SLIP_VELOCITY" : false,
                                     
                                     "VELOCITY_LAPLACIAN" : false,
                                     "BODY_FORCE" : false,
                                     
                                     "VELOCITY_LAPLACIAN_RATE" : false,
                                     "ELECTRIC_POTENTIAL": false,
                                     "FLUID_ION_DENSITY": false,
                                     "FLUID_ELECTRON_DENSITY": false,
                                     "FLUID_NEUTRAL_DENSITY": false,
                                     "ELECTRIC_FIELD": false,
                                     "MAGNETIC_FIELD": false
                                     },

            "print_distance_option" : false,

            "print_steps_per_plot_step" : 1,

            "problem_name" : "dummy_name1.Provide_a_real_one",

            "PredefinedSkinOption" : false,
            "MeanRadius" : 0.0001,

            "properties": [{
                "model_part_name": "dummy_name.Provide_a_real_one",
                "properties_id": 1,
                "plasma_dynamics_law_parameters": {
                    "name": "DEM_electromagnetic"}
                },
                {"electromagnetic_field_parameters":{
                    "external_electric_field_X": 0.0,
                    "external_electric_field_Y": 0.0,
                    "external_electric_field_Z": 0.0,
                    "external_magnetic_field_X": 0.0,
                    "external_magnetic_field_Y": 0.0,
                    "external_magnetic_field_Z": 0.0
                    }
                },
                {"density_parameters":{
                    "number_of_particles_in_a_ion_macroparticle": 1,
                    "particle_ion_density": 4.16e-22,
                    "fluid_neutral_density": 1.0,
                    "fluid_electron_density": 1.0
                    }
            }],

            "fluid_domain_volume"            : 1.0, 

            "fluid_parameters" : {},

            "dem_parameters" : {

                "problem_data"     : {
                    "problem_name"  : "Plasma_Dynamics_test_particles_in_a_virtual_cylinder_DEM",
                    "parallel_type" : "OpenMP",
                    "echo_level"    : 0,
                    "start_time"    : 0.0,
                    "end_time"      : 1
                },
                "solver_settings"   :{
                    "strategy" : "plasma_sphere_strategy",
                    "RemoveBallsInitiallyTouchingWalls": false

                },
                "do_print_results_option"          : false,
                "Dimension"                        : 3,
                "PeriodicDomainOption"             : false,
                "BoundingBoxOption"                : true,
                "AutomaticBoundingBoxOption"       : false,
                "BoundingBoxEnlargementFactor"     : 1.0,
                "BoundingBoxStartTime"             : 0.0,
                "BoundingBoxStopTime"              : 1000.0,
                "BoundingBoxMaxX"                  : 1000.0,
                "BoundingBoxMaxY"                  : 1000.0,
                "BoundingBoxMaxZ"                  : 1000.0,
                "BoundingBoxMinX"                  : -1000.0,
                "BoundingBoxMinY"                  : -1000.0,
                "BoundingBoxMinZ"                  : -1000.0,

                "dem_inlet_option"                 : false,
                "GravityX"                         : 0.0,
                "GravityY"                         : 0.0,
                "GravityZ"                         : 0.0,

                "VelocityTrapOption"               : false,
                "RotationOption"                   : false,
                "RemoveBallsInEmbeddedOption"      : true,

                "DeltaOption"                      : "Absolute",
                "SearchTolerance"                  : 0.0001,
                "CoordinationNumber"               : 10,
                "AmplifiedSearchRadiusExtension"   : 0.0,
                "ModelDataInfo"                    : false,
                "VirtualMassCoefficient"           : 1.0,
                "RollingFrictionOption"            : false,
                "DontSearchUntilFailure"           : false,
                "ContactMeshOption"                : false,
                "OutputFileType"                   : "Binary",
                "Multifile"                        : "multiple_files",

                "TranslationalIntegrationScheme"   : "Symplectic_Euler",
                "RotationalIntegrationScheme"      : "Direct_Integration",
                "DeltaTimeSafetyFactor"            : 1.0,
                "MaxTimeStep"                      : 1e-6,
                "FinalTime"                        : 1.0,
                "ControlTime"                      : 5.0,
                "NeighbourSearchFrequency"         : 1,
                "TestType"                         : "None",
                "ElementType"                      : "IonParticle3D",
                "problem_name"                     : "Plasma_Dynamics_test_particles_in_a_virtual_cylinder_DEM",
                "GraphExportFreq"                  : 1e-3,
                "VelTrapGraphExportFreq"           : 1e-3,
                "OutputTimeStep"                   : 0.01,
                "PostDisplacement"                 : true,
                "PostVelocity"                     : true,
                "PostElasticForces"                : false,
                "PostContactForces"                : false,
                "PostRigidElementForces"           : false,
                "PostTangentialElasticForces"      : false,
                "PostTotalForces"                  : false,
                "PostShearStress"                  : false,
                "PostNonDimensionalVolumeWear"     : false,
                "PostNodalArea"                    : false,
                "PostRHS"                          : false,
                "PostDampForces"                   : false,
                "PostAppliedForces"                : false,
                "PostRadius"                       : true,
                "PostGroupId"                      : false,
                "PostExportId"                     : false,
                "PostAngularVelocity"              : false,
                "PostParticleMoment"               : false,
                "PostEulerAngles"                  : false,
                "PostBoundingBox"                  : false                
            },


            "frame_of_reference" : {
                "frame_type": 0,
                "angular_velocity_of_frame_X" : 0.0,
                "angular_velocity_of_frame_Y" : 0.0,
                "angular_velocity_of_frame_Z" : 0.0,
                "angular_velocity_of_frame_old_X" : 0.0,
                "angular_velocity_of_frame_old_Y" : 0.0,
                "angular_velocity_of_frame_old_Z" : 0.0,
                "acceleration_of_frame_origin_X" : 0.0,
                "acceleration_of_frame_origin_Y" : 0.0,
                "acceleration_of_frame_origin_Z" : 0.0,
                "angular_acceleration_of_frame_X" : 0.0,
                "angular_acceleration_of_frame_Y" : 0.0,
                "angular_acceleration_of_frame_Z" : 0.0,
                "frame_rotation_axis_initial_point" : [0.0,0.0,0.0],
                "frame_rotation_axis_final_point" : [0.0,0.0,1.0],
                "angular_velocity_magnitude" : 0.0
            },

            "custom_dem" : {
                "do_solve_dem" : true,
                "do_search_neighbours" : false,
                "type_of_dem_inlet" : "VelocityImposed",
                "type_of_dem_inlet_comment" : "VelocityImposed or ForceImposed",
                "translational_integration_scheme" : "Symplectic_Euler"
            },            

            "custom_fluid" : {
                "fluid_already_calculated" : false,
                "embedded_option" : false,
                "embedded_option_comment" : "the embedded domain tools are to be used",
                "do_impose_flow_from_field_option" : false,
                "body_force_on_fluid_option" : true,
                "ALE_option" : false,
                "fluid_model_type" : 1,
                "fluid_model_type_comment" : " untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)"
            },

            "stationarity" : {
                "stationary_problem_option" : false,
                "stationary_problem_option_comment" : " stationary, stop calculating the fluid after it reaches the stationary state",
                "max_pressure_variation_rate_tol" : 1e-3,
                "max_pressure_variation_rate_tol_comment": " for stationary problems, criterion to stop the fluid calculations",
                "time_steps_per_stationarity_step" : 15,
                "time_steps_per_stationarity_step_comment": " number of fluid time steps between consecutive assessment of stationarity steps",
                "time_steps_per_analytic_processing_step": 1
            },

            "output_interval" : 0.5,

            "coupling" : {
                "coupling_level_type" : 1,
                "coupling_weighing_type" : 2,
                "coupling_weighing_type_comment" : "{fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)",
                "interaction_start_time" : 0.0,
                "time_averaging_type" : 0,

                "backward_coupling" : {
                    "meso_scale_length" : 0.2,
                    "meso_scale_length_comment" : " the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)",
                    "shape_factor" : 0.5,
                    "filter_velocity_option" : false,
                    "apply_time_filter_to_fluid_fraction_option" : false,
                    "min_fluid_fraction" : 0.2,
                    "fluid_fraction_grad_type" : 0,
                    "calculate_diffusivity_option" : false,
                    "viscosity_modification_type" : 0
                }
            },


            "time_stepping" : {
                "automatic_time_step" : false,
                "time_step"           : 1e-5
            }

            }""")

    return default_settings
