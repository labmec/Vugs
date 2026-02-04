# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)

project (MyProject)


find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)



add_subdirectory(Targets)

#add_subdirectory(TestNeoPz)

#add_subdirectory(TestOpenCV)
#add_subdirectory(TestBinaryImage)
#add_subdirectory(RAW3dVTK)


