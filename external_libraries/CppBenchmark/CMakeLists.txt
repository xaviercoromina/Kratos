SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Include the header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)

## HDR

# Generate variables with the sources
FILE(GLOB_RECURSE HDR_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/hdr/*.c") #all .cpp

# Compile library
ADD_LIBRARY(hdr STATIC ${HDR_SOURCES})

# Specify the include directories for the library's users
TARGET_INCLUDE_DIRECTORIES(hdr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Define include and library
SET(HDR_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/hdr)
SET(HDR_LIBRARIES hdr)

# This compiles the "DLL" interface (C API)
TARGET_COMPILE_DEFINITIONS(hdr PRIVATE HDR_DLL)

# Shared libraries
TARGET_COMPILE_DEFINITIONS(hdr PRIVATE HDR_EXPORTS)

# Add definitions
ADD_DEFINITIONS( -DHDR_EXPORTS )

# Install targets
INSTALL(TARGETS hdr DESTINATION libs)

## CPPBENCHMARK

# Generate variables with the sources
FILE(GLOB_RECURSE CPPBENCHMARK_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/benchmark/*.cpp") #all .cpp

# Compile library
ADD_LIBRARY(cppbenchmark STATIC ${CPPBENCHMARK_SOURCES})
TARGET_LINK_LIBRARIES(cppbenchmark PRIVATE ${HDR_LIBRARIES})

# Specify the include directories for the library's users
TARGET_INCLUDE_DIRECTORIES(cppbenchmark PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Define include and library
SET(CPPBENCHMARK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include/benchmark)
SET(CPPBENCHMARK_LIBRARIES cppbenchmark)

# This compiles the "DLL" interface (C API)
TARGET_COMPILE_DEFINITIONS(cppbenchmark PRIVATE CPPBENCHMARK_DLL)

# Shared libraries
TARGET_COMPILE_DEFINITIONS(cppbenchmark PRIVATE CPPBENCHMARKLIB_EXPORTS)

# Add definitions
ADD_DEFINITIONS( -DCPPBENCHMARKLIB_EXPORTS )

# Install targets
INSTALL(TARGETS cppbenchmark DESTINATION libs)

# Printing a message
MESSAGE(STATUS "cppbenchmark : ${CPPBENCHMARK_INCLUDE_DIR}" )