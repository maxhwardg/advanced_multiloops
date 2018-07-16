cmake_minimum_required(VERSION 3.5)
project(librnary)

find_package(Threads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pthread")


set(RNASTRUCTURE_SOURCES
        lib/RNAstructure/rna_library.cpp lib/RNAstructure/rna_library.h lib/RNAstructure/structure.cpp
        lib/RNAstructure/structure.h lib/RNAstructure/stackclass.cpp lib/RNAstructure/stackclass.h
        lib/RNAstructure/stackstruct.cpp lib/RNAstructure/stackstruct.h lib/RNAstructure/algorithm.cpp
        lib/RNAstructure/algorithm.h lib/RNAstructure/arrayclass.cpp lib/RNAstructure/arrayclass.h
        lib/RNAstructure/dotarray.cpp lib/RNAstructure/dotarray.h lib/RNAstructure/forceclass.cpp
        lib/RNAstructure/forceclass.h lib/RNAstructure/TProgressDialog.cpp lib/RNAstructure/TProgressDialog.h
        lib/RNAstructure/defines.h lib/RNAstructure/platform.h)


set(LIBRNARY_SOURCES
        include/primary_structure.hpp src/primary_structure.cpp include/secondary_structure.hpp
        src/secondary_structure.cpp include/ss_tree.hpp src/ss_tree.cpp src/rna_structure.cpp include/rna_structure.hpp
        src/ss_enumeration.cpp include/ss_enumeration.hpp include/statistics.hpp src/statistics.cpp include/models/nn_model.hpp
        src/models/nn_model.cpp include/scorers/nn_scorer.hpp include/models/nn_affine_model.hpp src/models/nn_affine_model.cpp
        include/models/nn_unpaired_model.hpp src/models/nn_unpaired_model.cpp include/folders/brute_folder.hpp
        include/folders/nn_unpaired_folder.hpp
        src/folders/nn_unpaired_folder.cpp include/folders/brute_folder.hpp include/read_cts.hpp src/read_cts.cpp
        include/models/aalberts_model.hpp src/models/aalberts_model.cpp include/folders/upbr_folder.hpp
        include/multi_array.hpp
        src/folders/upbr_folder.cpp include/scorers/aalberts_scorer.hpp include/folders/aalberts_folder.hpp
        src/folders/aalberts_folder.cpp
        include/vector_types.hpp include/folders/nussinov_folder.hpp
        src/folders/nussinov_folder.cpp include/pseudoknot_removal.hpp src/pseudoknot_removal.cpp
        include/models/average_asym_model.hpp include/folders/average_asym_folder.hpp
        src/models/average_asym_model.cpp src/folders/average_asym_folder.cpp src/multi_loop.cpp include/multi_loop.hpp
        src/energy.cpp
        src/scorers/aalberts_scorer.cpp include/folders/nn_affine_folder.hpp
        src/folders/nn_affine_folder.cpp include/training/IBF_multiloop.hpp
        include/parallel.hpp
        include/training/IBF_multiloop_aalberts.hpp)

set(SOURCE_FILES "${RNASTRUCTURE_SOURCES}" "${LIBRNARY_SOURCES}")

set(PROJECT_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/lib")

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_INCLUDE_DIRS}")

target_link_libraries (${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})
