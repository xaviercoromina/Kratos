#include <iostream>
#include <string>

#include "custom_workflows/dgeosettlement.h"

using namespace Kratos;


int main()
{
    std::cout << "About to run a single stage of a settlement calculation" << std::endl;

    const std::string working_directory{"C:/Users/graaf/tmp/settlement_test"};
    const std::string project_parameters_file_name{"ProjectParameters_stage1.json"};

    auto dummy_log_callback = [](char*){};
    auto dummy_report_progress = [](double){};
    auto dummy_report_textual_progress = [](char*){};
    auto never_cancel = [](){ return false; };

    KratosGeoSettlement geo_settlement;
    geo_settlement.RunStage(working_directory, project_parameters_file_name, dummy_log_callback,
                            dummy_report_progress, dummy_report_textual_progress, never_cancel);

    return 0;
}
