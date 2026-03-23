# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 



add_executable(Comparacao main.cpp)
#target_link_libraries(Comparacao PRIVATE NeoPZ::pz)

add_executable(Test3d main3D.cpp)
target_link_libraries(Comparacao PRIVATE libsource)
target_link_libraries(Comparacao PRIVATE NeoPZ::pz)

target_link_libraries(Test3d PRIVATE libsource)
target_link_libraries(Test3d PRIVATE NeoPZ::pz)
