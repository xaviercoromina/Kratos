//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez
//

#ifndef KRATOS_COMPUTE_DISTANCE_SENSITIVITIES_PROCESS_H
#define KRATOS_COMPUTE_DISTANCE_SENSITIVITIES_PROCESS_H

#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

	class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ComputeDistanceSensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ComputeDistanceSensitivitiesProcess);

    // Constructor for ComputeDistanceSensitivitiesProcess Process
    ComputeDistanceSensitivitiesProcess(ModelPart& rModelPart, Parameters ThisParameters);

    /// Destructor.
    ~ComputeDistanceSensitivitiesProcess() = default;

    /// Assignment operator.
    ComputeDistanceSensitivitiesProcess& operator=(ComputeDistanceSensitivitiesProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeDistanceSensitivitiesProcess(ComputeDistanceSensitivitiesProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    void Execute() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ComputeDistanceSensitivitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeDistanceSensitivitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;

}; // Class Process
} // namespace Kratos


#endif // KRATOS_CALCULATE_DISTANCE_SENSITIVITIES_PROCESS_H
