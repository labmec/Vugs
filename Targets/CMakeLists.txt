# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 




find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)



add_subdirectory(ComparacaoH1xHdiv)
add_subdirectory(ContornosVugs)
add_subdirectory(SimulacaoPconstH1)
add_subdirectory(SimulacaoPconstHdiv)
#add_subdirectory(FormaResidual)



