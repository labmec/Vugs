# Adds the current source directory to the include path for libsource.
target_include_directories(libsource PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Collects all C++ source files (*.cpp, .c) and header files (.hpp) in this directory into the SRC variable.
file(GLOB SRC CONFIGURE_DEPENDS "*.cpp" "*.c" "*.hpp")
# Collects all header files (*.h) into the HEADERS variable.
file(GLOB HEADERS CONFIGURE_DEPENDS "*.h")

# Add all the source and header files to the libsource library.
target_sources(libsource PRIVATE ${SRC} ${HEADERS})