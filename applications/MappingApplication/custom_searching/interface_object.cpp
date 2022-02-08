//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// Project includes
#include "utilities/openmp_utils.h"
#include "interface_object.h"

namespace Kratos {

InterfaceObject::InterfaceObject(const CoordinatesArrayType& rCoordinates)
    : Point(rCoordinates) { }

void InterfaceObject::InitializeThreadLocalStorage()
{
    mThreadFlags.resize(ParallelUtilities::GetNumThreads(), true);
}

bool InterfaceObject::IsEnabledInThisThread() const
{
    return mThreadFlags[OpenMPUtils::ThisThread()];
}

void InterfaceObject::EnableInThisThread()
{
    mThreadFlags[OpenMPUtils::ThisThread()] = true;
}

void InterfaceObject::DisableInThisThread()
{
    mThreadFlags[OpenMPUtils::ThisThread()] = false;
}

}  // namespace Kratos.
