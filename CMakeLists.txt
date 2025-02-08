cmake_minimum_required(VERSION 3.10)

# Set the project name
project(DFHydro)
# Set the C++ standard to c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ad an include directory include
include_directories(include)

# List the sources 
set(SOURCES src/jsoncpp.cpp)

# Add an executable based on the sources
add_executable(test_json.exe src/test_json.cpp ${SOURCES})
