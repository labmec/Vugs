# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 


add_executable(TestHDiv main.cpp)
add_executable(TestFrac mainFrac.cpp)
target_link_libraries(TestHDiv PRIVATE libsource)
target_link_libraries(TestHDiv PRIVATE NeoPZ::pz)

target_link_libraries(TestFrac PRIVATE libsource)
target_link_libraries(TestFrac PRIVATE NeoPZ::pz)