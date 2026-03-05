# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 


find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

add_custom_target(GlobalFiles SOURCES ${PZ_SRC})

# Creates a shared library (dynamic library)
add_library(libsource SHARED "")

# Links the NeoPZ library to libsource publicly, meaning any project that links 
# against libsource will also link against NeoPZ.
target_link_libraries(libsource PUBLIC NeoPZ::pz)

# Adds the NeoPZ include directories to libsource so it can find the headers
target_include_directories(libsource PUBLIC ${PZ_INCLUDE_DIRS})

add_subdirectory(ComparacaoH1xHdiv)
add_subdirectory(ContornosVugs)
add_subdirectory(SimulacaoPconstH1)
add_subdirectory(SimulacaoPconstHdiv)
add_subdirectory(FormaResidual)
add_subdirectory(Sources)



