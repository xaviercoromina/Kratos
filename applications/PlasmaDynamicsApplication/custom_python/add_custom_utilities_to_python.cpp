//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Chung To Sang
//


// External includes

// Project includes
#include "add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"
#include "custom_utilities/renumbering_nodes_utility_for_plasma_dynamics.h"
#include "custom_utilities/binbased_DEM_fluid_coupled_mapping_for_plasma_dynamics.h"
#include "custom_elements/ion_particle.h"
#include "custom_elements/electron_particle.h"


namespace Kratos
{

namespace Python
{




void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<RenumberingNodesUtilityForPlasmaDynamics> (m, "RenumberingNodesUtilityForPlasmaDynamics")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def(py::init<ModelPart&,ModelPart&,ModelPart&,ModelPart&,ModelPart&>())
        .def("Renumber", &RenumberingNodesUtilityForPlasmaDynamics::Renumber)
        .def("UndoRenumber", &RenumberingNodesUtilityForPlasmaDynamics::UndoRenumber)
        ;  

    py::class_<BinBasedDEMFluidCoupledMappingForPlasmaDynamics> (m, "BinBasedDEMFluidCoupledMappingForPlasmaDynamics3D")
        .def(py::init<Parameters&>())
        .def("InterpolateVelocityOnAuxVelocity", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::InterpolateVelocityOnAuxVelocity)
        .def("ImposeVelocityOnDEMFromFieldToAuxVelocity", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::ImposeVelocityOnDEMFromFieldToAuxVelocity)
        .def("InterpolateFromFluidMesh", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::InterpolateFromFluidMesh)
        .def("ImposeFlowOnDEMFromField", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::ImposeFlowOnDEMFromField)
        .def("InterpolateFromDEMMesh", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::InterpolateFromDEMMesh)
        .def("HomogenizeFromDEMMesh", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::HomogenizeFromDEMMesh)
        .def("ComputePostProcessResults", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::ComputePostProcessResults)
        .def("AddDEMCouplingVariable", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::AddDEMCouplingVariable)
        .def("AddFluidCouplingVariable", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::AddFluidCouplingVariable)
        .def("AddDEMVariablesToImpose", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::AddDEMVariablesToImpose)
        .def("AddFluidVariableToBeTimeFiltered", &BinBasedDEMFluidCoupledMappingForPlasmaDynamics ::AddFluidVariableToBeTimeFiltered)
        ;

    py::class_<IonParticle> (m, "IonParticle")
        .def("Calculate", &IonParticle ::Calculate)
        ;

    py::class_<ElectronParticle> (m, "ElectronParticle")
        .def("Calculate", &ElectronParticle ::Calculate)
        ;
    


}

}  // namespace Python.
} // Namespace Kratos
