SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Generate variables with the sources
FILE(GLOB_RECURSE CPPBENCHMARK_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp") #all .cpp

# Include the header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)

# Compile library
ADD_LIBRARY(cppbenchmark STATIC ${CPPBENCHMARK_SOURCES})

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