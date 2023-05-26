//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jonathan Nuttall
//
#include "set_multiple_moving_loads.h"

#include <custom_conditions/moving_load_condition.h>
#include "includes/variables.h"



namespace Kratos
{

    SetMultipleMovingLoadsProcess::SetMultipleMovingLoadsProcess(ModelPart& rModelPart,
        Parameters Settings)
        : mrModelPart(rModelPart),
        mParameters(Settings)
    {
    	Parameters default_parameters(R"(
        {
            "help"                    : "This process applies multiple coupled moving load conditions belonging to a model part. The loads move over line elements.",
            "model_part_name"         : "please_specify_model_part_name",
			"compute_model_part_name" : "please_specify_compute_body_part_name",
            "variable_name"           : "POINT_LOAD",
            "load"                    : [0.0, 1.0, 0.0],
            "direction"               : [1,1,1],
            "velocity"                : 1,
			"origin"                  : [0.0, 0.0, 0.0],
			"configuration"           : [0.0],
			"function_path"           : "please specify a string to the UVEC function path",
			"uvec_parameters"         : {"parameters":{}, "state":{}}
        }  )"
        );


    	// Set default velocity as a string, if the input velocity is a string
        if (mParameters.Has("velocity")) {
            if (mParameters["velocity"].IsString()) {
                default_parameters["velocity"].SetString("1");
            }
        }

        // ToDo JDN: Needs dealing with (also in python)
        //mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // check if load parameter has size 3
        KRATOS_ERROR_IF(mParameters["load"].size() != 3) <<
            "'load' has to have size 3!" << std::endl;

        // check if all elements in load parameter are either string or a number
        bool is_all_string = true;
        bool is_all_number = true;
        for (IndexType i = 0; i < mParameters["load"].size(); i++)
        {
            if (!mParameters["load"][i].IsString()) {
                is_all_string = false;
            }
            if (!mParameters["load"][i].IsNumber()) {
                is_all_number = false;
            }
        }

        KRATOS_ERROR_IF(!is_all_string && !is_all_number) << "'load' has to be a vector of numbers, or an array with strings" << std::endl;

        
        int count = 0;
    	for (double offset : mParameters["configuration"].GetVector())
        {

            auto parameters_moving_load = mParameters.Clone();

    		count++;
            mConfigurationCallback.push_back(false);
            const std::string newModelPartName = mrModelPart.Name() + "_cloned_" + std::to_string(count);
            auto& new_cloned_model_part = CloneMovingConditionInComputeModelPart(newModelPartName, count);

            parameters_moving_load.RemoveValue("configuration");
            parameters_moving_load.RemoveValue("compute_model_part_name");
            parameters_moving_load.RemoveValue("function_path");
            parameters_moving_load.RemoveValue("uvec_parameters");
    		parameters_moving_load.AddDouble("offset", offset);
    		auto r_moving_point_process = SetMovingLoadProcess(new_cloned_model_part, parameters_moving_load);
            mMovingPointLoadsProcesses.push_back(r_moving_point_process);
           
        }

        RemoveClonedConditions();

    }

    ModelPart& SetMultipleMovingLoadsProcess::CloneMovingConditionInComputeModelPart(std::string NewBodyPartName, int ConfigurationIndex)
    {
        auto& compute_model_part = mrModelPart.GetRootModelPart().GetSubModelPart(mParameters["compute_model_part_name"].GetString());
        auto& new_model_part = compute_model_part.CreateSubModelPart(NewBodyPartName);
        new_model_part.AddNodes(mrModelPart.NodesBegin(), mrModelPart.NodesEnd());

    	int index = GetMaxConditionsIndex();
    	for(auto& moving_load_condition: mrModelPart.Conditions())
        {
            index++;
    		const Condition::Pointer clone_condition = moving_load_condition.Clone(index, moving_load_condition.GetGeometry());
    		auto clone_condition_ptr = dynamic_cast<MovingLoadCondition<2, 2>*>(clone_condition.get());
            if (clone_condition_ptr == nullptr) KRATOS_ERROR << "SetMultipleMovingLoadsProcess:" << "Could not cast condition to MovingLoadCondition" << std::endl;

    		if (mPythonUvecFunction == nullptr) SetPythonUvecFunction();
    		clone_condition_ptr->SetUvecFunction(
                std::bind(&SetMultipleMovingLoadsProcess::PythonInitializeNonLinearFunction, this,
                    std::placeholders::_1, std::placeholders::_2, std::placeholders::_3));
            clone_condition_ptr->SetConfigurationIndex(ConfigurationIndex-1);
    		new_model_part.AddCondition(clone_condition);
        }
        return new_model_part;
    }

    void SetMultipleMovingLoadsProcess::SetPythonUvecFunction()
    {
        KRATOS_INFO("SetPythonUvecFunction") << "Setting Up" << std::endl;
    	PyObject* python_uvec_function_module = PyImport_ImportModule(mParameters["function_path"].GetString().c_str());
        mPythonUvecFunction = PyObject_GetAttrString(python_uvec_function_module, "UVEC");
        KRATOS_INFO("SetPythonUvecFunction") << "Complete" << std::endl;
    }


	void SetMultipleMovingLoadsProcess::PythonInitializeNonLinearFunction(int configurationIndex, Vector u, Vector theta)
    {
        mConfigurationCallback[configurationIndex] = true;
		
		//Update_json

        // Displacements per axle
    	if(!mParameters["uvec_parameters"]["u"].Has(std::to_string(configurationIndex)))
        {
            mParameters["uvec_parameters"]["u"].AddVector(std::to_string(configurationIndex), u);
        }
        else
        {
            mParameters["uvec_parameters"]["u"][std::to_string(configurationIndex)].SetVector(u);
        }

        // Rotations per axle
        if (!mParameters["uvec_parameters"]["theta"].Has(std::to_string(configurationIndex)))
        {
            mParameters["uvec_parameters"]["theta"].AddVector(std::to_string(configurationIndex), theta);
        }
        else
        {
            mParameters["uvec_parameters"]["theta"][std::to_string(configurationIndex)].SetVector(theta);
        }

    	for (auto configuration : mConfigurationCallback)
        {
        	if (!configuration)
            {
                return;
            }
        }

        mParameters["uvec_parameters"]["dt"].SetDouble(mrModelPart.GetProcessInfo()[DELTA_TIME]);

        PyObject* py_string = PyUnicode_FromString(mParameters["uvec_parameters"].WriteJsonString().c_str());
    	PyObject* result = PyObject_CallFunction(mPythonUvecFunction, "O", py_string);

    	if (result)
        {
            const char* jsonStr = PyUnicode_AsUTF8(result);
            Parameters kp = Parameters(jsonStr);
            mParameters["uvec_parameters"] = kp;
        }

        // update the moving load on all conditions
        int configurationNo = 0;
        for (double offset : mParameters["configuration"].GetVector()) // Improve this loop
        {
            configurationNo++;
            const std::string newModelPartName = mrModelPart.Name() + "_cloned_" + std::to_string(configurationNo);
            auto& compute_model_part = mrModelPart.GetRootModelPart().GetSubModelPart(mParameters["compute_model_part_name"].GetString());
            auto& relevantModelPart = compute_model_part.GetSubModelPart(newModelPartName);
            for (Condition condition : relevantModelPart.Conditions())
            {
                condition.SetValue(POINT_LOAD, mParameters["uvec_parameters"]["loads"][std::to_string((configurationNo - 1))].GetVector());
            }
        }
        for (auto configuration : mConfigurationCallback)
        {
            configuration = false;
        }

    }

    int SetMultipleMovingLoadsProcess::GetMaxConditionsIndex()
    {
        auto& root_model_part = mrModelPart.GetRootModelPart();
        int max_index = 0;
    	for(auto& condition: root_model_part.Conditions())
        {
            max_index = std::max(max_index, static_cast<int>(condition.Id()));
        }
        return max_index;
    }

    void SetMultipleMovingLoadsProcess::RemoveClonedConditions()
    {
        auto& compute_model_part = mrModelPart.GetRootModelPart().GetSubModelPart(mParameters["compute_model_part_name"].GetString());

        KRATOS_INFO("RemoveClonedModelPart: Compute") << compute_model_part.Conditions().size() << std::endl;
        for (auto& moving_load_condition : mrModelPart.Conditions())
        {
            compute_model_part.pGetCondition(moving_load_condition.Id())->Set(TO_ERASE, true);
        }

        // Call method
        compute_model_part.RemoveConditions(TO_ERASE);
        KRATOS_INFO("RemoveClonedModelPart: Compute") << compute_model_part.Conditions().size() << std::endl;

    }

    void SetMultipleMovingLoadsProcess::ExecuteInitialize()
    {
        for (Process& rMovingPointLoad : mMovingPointLoadsProcesses)
        {
            rMovingPointLoad.ExecuteInitialize();
        }
    }

    void SetMultipleMovingLoadsProcess::ExecuteInitializeSolutionStep()
    {
		for (Process& rMovingPointLoad: mMovingPointLoadsProcesses)
		{
            rMovingPointLoad.ExecuteInitializeSolutionStep();
		}
    };

    void SetMultipleMovingLoadsProcess::ExecuteFinalizeSolutionStep()
    {
        for (Process& rMovingPointLoad : mMovingPointLoadsProcesses)
        {
            rMovingPointLoad.ExecuteFinalizeSolutionStep();
        }

    };

}