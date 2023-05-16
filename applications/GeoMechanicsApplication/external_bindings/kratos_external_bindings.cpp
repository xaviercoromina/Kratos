#define EXPORT __declspec(dllexport)
#include "kratos_external_bindings.h"

extern "C"
{

#if defined(KRATOS_COMPILED_IN_WINDOWS)

    EXPORT Kratos::KratosExecute *KratosExecute_CreateInstance()
    {
        return new Kratos::KratosExecute();
    }

    EXPORT int __stdcall execute_flow_analysis(Kratos::KratosExecute *instance,
                                               char *workingDirectory,
                                               char *projectFile,
                                               double minCriticalHead,
                                               double maxCriticalHead,
                                               double stepCriticalHead,
                                               char *criticalHeadBoundaryModelPartName,
                                               void __stdcall logCallback(char *),
                                               void __stdcall reportProgress(double),
                                               void __stdcall reportTextualProgress(char *),
                                               bool __stdcall shouldCancel())
    {
        int errorCode = instance->execute_flow_analysis(workingDirectory,
                                                        projectFile,
                                                        minCriticalHead,
                                                        maxCriticalHead,
                                                        stepCriticalHead,
                                                        criticalHeadBoundaryModelPartName,
                                                        logCallback,
                                                        reportProgress,
                                                        reportTextualProgress,
                                                        shouldCancel);
        return errorCode;
    }



    EXPORT Kratos::KratosGeoSettlement* KratosGeoSettlement_CreateInstance()
    {
        return new Kratos::KratosGeoSettlement{};
    }

    EXPORT int __stdcall runSettlementStage(Kratos::KratosGeoSettlement* instance,
                                            char*                        workingDirectory,
                                            char*                        projectFileName,
                                            void __stdcall               logCallback(char*),
                                            void __stdcall               reportProgress(double),
                                            void __stdcall               reportTextualProgress(char*),
                                            bool __stdcall               shouldCancel())
    {
        return instance->RunStage(workingDirectory,
                                  projectFileName,
                                  logCallback,
                                  reportProgress,
                                  reportTextualProgress,
                                  shouldCancel);
    }

#endif
}
