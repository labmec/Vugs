# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 



add_executable(TestH1 main.cpp)
target_link_libraries(TestH1 PRIVATE NeoPZ::pz)
target_link_libraries(TestH1 PRIVATE libsource)

add_executable(TestFracH1 mainFrac.cpp)
target_link_libraries(TestFracH1 PRIVATE NeoPZ::pz)
target_link_libraries(TestFracH1 PRIVATE libsource)

