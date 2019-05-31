set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPlasmaDynamicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SwimmingDEMApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidTransportApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/PlasmaDynamicsApplication )

## Define application core sources (pure C++ parts)
set( KRATOS_PLASMA_DYNAMICS_APPLICATION_CORE_SOURCES

    ## Application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/plasma_dynamics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/plasma_dynamics_application_variables.cpp

    # Utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/binbased_DEM_fluid_coupled_mapping_for_plasma_dynamics.cpp
    # Processes

    # Constitutive Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_electromagnetic_CL.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ion_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/electron_particle.cpp
)

## Define python interface sources
set( KRATOS_PLASMA_DYNAMICS_APPLICATION_PYTHON_INTERFACE_SOURCES

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/plasma_dynamics_python_application.cpp
)

## PlasmaDynamicsApplication core library (C++ parts)
add_library( KratosPlasmaDynamicsCore SHARED ${KRATOS_PLASMA_DYNAMICS_APPLICATION_CORE_SOURCES})
target_link_libraries( KratosPlasmaDynamicsCore PUBLIC KratosCore KratosSwimmingDEMCore KratosDEMCore KratosFluidDynamicsCore KratosFluidTransportCore)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## PlasmaDynamicsApplication core library (C++ parts)
set_target_properties( KratosPlasmaDynamicsCore PROPERTIES COMPILE_DEFINITIONS "PLASMA_DYNAMICS_APPLICATION=EXPORT,API")

## PlasmaDynamicsApplication python module
pybind11_add_module(KratosPlasmaDynamicsApplication MODULE THIN_LTO ${KRATOS_PLASMA_DYNAMICS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosPlasmaDynamicsApplication PRIVATE KratosPlasmaDynamicsCore)
set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosPlasmaDynamicsCore)
    cotire(KratosPlasmaDynamicsApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosPlasmaDynamicsCore DESTINATION libs )
install(TARGETS KratosPlasmaDynamicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PlasmaDynamicsApplication.py" DESTINATION "KratosMultiphysics/PlasmaDynamicsApplication" RENAME "__init__.py")


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
