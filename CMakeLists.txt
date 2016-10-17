cmake_minimum_required(VERSION 3.5)
project(advanced_multiloops)


SET(CMAKE_CXX_FLAGS_ASAN "-fsanitize=address -lasan -g" CACHE STRING
        "Flags used by the C++ compiler during ASAN builds"
        FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_ASAN
        "-fsanitize=address -lasan" CACHE STRING
        "Flags used for linking binaries during ASAN builds."
        FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
        "-fsanitize=address -lasan" CACHE STRING
        "Flags used by the shared libraries linker during ASAN builds."
        FORCE )
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_MAINTAINER
        CMAKE_EXE_LINKER_FLAGS_MAINTAINER
        CMAKE_SHARED_LINKER_FLAGS_MAINTAINER )
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel ASAN."
        FORCE )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pedantic -Wshadow")

add_subdirectory(librnary)
add_subdirectory(programs)