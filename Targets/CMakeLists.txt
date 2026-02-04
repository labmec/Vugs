# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 



add_executable(Targets main.cpp)


set(All_Dirs  NeoPZ::pz  ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} )

target_include_directories(Targets PUBLIC All_Dirs)
