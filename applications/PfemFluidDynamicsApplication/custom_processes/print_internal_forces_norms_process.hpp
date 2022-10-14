//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics PfemFluidDynamics Application
//
//  License:         BSD License
//    Kratos default license:
//  kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include <fstream>
#include <iostream>

namespace Kratos
{

  ///@name Kratos Classes
  ///@{

  /** This class retrieves internal forces vector norms and prints it in a file
   */
  class PrintInternalForcesNormProcess : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PrintInternalForcesNormProcess
    KRATOS_CLASS_POINTER_DEFINITION(PrintInternalForcesNormProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    PrintInternalForcesNormProcess(ModelPart &rModelPart, const std::string OutputFileName = "WaveHeight", const double TimeInterval = 0.0) : mrModelPart(rModelPart), mOutputFileName(OutputFileName), mTimeInterval(TimeInterval)
    {
      std::ofstream my_file;
      const std::string file_name = mOutputFileName + ".txt";
      my_file.open(file_name, std::ios_base::trunc);
      my_file << "    TIME     INERTIAL_FORCES_NORM     VISCOUS_FORCES_NORM     VOLUMETRIC_FORCES_NORM      EXTERNAL_FORCES_NORM" << std::endl;
      my_file.close();
    }

    /// Destructor.
    virtual ~PrintInternalForcesNormProcess() {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
      Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void ExecuteFinalizeSolutionStep() override
    {
      KRATOS_TRY
      const double time = mrModelPart.GetProcessInfo()[TIME];
      const int step = mrModelPart.GetProcessInfo()[STEP];
      const std::string file_name = mOutputFileName + ".txt";

      // Initialize values to be printed
      double inertial_forces = 0.0;
      double viscous_forces = 0.0;
      double vol_forces = 0.0;
      double ext_forces = 0.0;

      if (time - mPreviousPlotTime > mTimeInterval || step == 1) {
        // We loop over the elements...
        const auto it_elem_begin = mrModelPart.ElementsBegin();

        //#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); ++i) {
          auto it_elem = it_elem_begin + i;

          inertial_forces += it_elem->GetValue(INERTIAL_FORCES_NORM);
          viscous_forces += it_elem->GetValue(VISCOUS_FORCES_NORM);
          vol_forces += it_elem->GetValue(VOLUMETRIC_FORCES_NORM);
          ext_forces += it_elem->GetValue(EXTERNAL_FORCES_NORM);
        }

        // We open the file where we print the wave height values
        std::ofstream my_file;
        my_file.open(file_name, std::ios_base::app);
        my_file << "  " + std::to_string(time) + "    " + std::to_string(inertial_forces)+ "    " + std::to_string(viscous_forces)+ "    " + std::to_string(vol_forces)+ "    " + std::to_string(ext_forces) << std::endl;
        mPreviousPlotTime = time;

      }
      KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
      return "PrintInternalForcesNormProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "PrintInternalForcesNormProcess";
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    PrintInternalForcesNormProcess(PrintInternalForcesNormProcess const &rOther);

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

  private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;

    std::string mOutputFileName;
    double mTimeInterval;
    double mPreviousPlotTime = 0.0;


    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    PrintInternalForcesNormProcess &operator=(PrintInternalForcesNormProcess const &rOther);

    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

  }; // Class PrintInternalForcesNormProcess

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  inline std::istream &operator>>(std::istream &rIStream,
                                  PrintInternalForcesNormProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const PrintInternalForcesNormProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.
