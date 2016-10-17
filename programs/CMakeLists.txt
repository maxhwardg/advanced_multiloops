cmake_minimum_required(VERSION 3.5)
project(programs)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")


add_executable(linear_fold src/linear_fold.cpp)
target_include_directories(linear_fold PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(linear_fold librnary)

add_executable(linear_efn src/linear_efn.cpp)
target_include_directories(linear_efn PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(linear_efn librnary)

add_executable(logarithmic_fold src/logarithmic_fold.cpp)
target_include_directories(logarithmic_fold PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(logarithmic_fold librnary)

add_executable(logarithmic_efn src/logarithmic_efn.cpp)
target_include_directories(logarithmic_efn PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(logarithmic_efn librnary)


add_executable(an_fold src/an_fold.cpp)
target_include_directories(an_fold PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(an_fold librnary)

add_executable(an_efn src/an_efn.cpp)
target_include_directories(an_efn PRIVATE "${INCLUDE_DIRS}")
target_link_libraries(an_efn librnary)