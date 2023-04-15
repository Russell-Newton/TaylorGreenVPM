cmake_minimum_required(VERSION 3.25)
project(TaylorGreenVPM CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

add_executable(TaylorGreenVPM main.cpp)
target_link_libraries(TaylorGreenVPM PUBLIC vpm)

add_executable(TaylorGreenVPMValidation validation.cpp)
target_link_libraries(TaylorGreenVPMValidation PUBLIC vpm)
