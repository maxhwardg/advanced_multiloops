cmake_minimum_required(VERSION 3.5)
project(librnary)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")


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
        src/ss_enumeration.cpp include/ss_enumeration.hpp include/statistics.hpp src/statistics.cpp include/nn_model.hpp
        src/nn_model.cpp include/nn_scorer.hpp include/nn_affine_model.hpp src/nn_affine_model.cpp
        include/nn_unpaired_model.hpp src/nn_unpaired_model.cpp include/brute_folder.hpp include/nn_unpaired_folder.hpp
        src/nn_unpaired_folder.cpp include/brute_folder.hpp include/read_cts.hpp src/read_cts.cpp
        include/aalberts_model.hpp src/aalberts_model.cpp include/upbr_folder.hpp include/multi_array.hpp
        src/upbr_folder.cpp include/aalberts_scorer.hpp include/aalberts_folder.hpp src/aalberts_folder.cpp
        include/vector_types.hpp include/asymmetry_model.hpp src/asymmtery_model.cpp include/asymmetry_scorer.hpp
        include/asymmetry_folder.hpp include/nussinov_folder.hpp
        src/nussinov_folder.cpp include/pseudoknot_removal.hpp src/pseudoknot_removal.cpp
        include/average_asym_model.hpp include/average_asym_folder.hpp src/asymmetry_folder.cpp
        src/average_asym_model.cpp src/average_asym_folder.cpp src/multi_loop.cpp include/multi_loop.hpp src/energy.cpp
        include/energy.hpp include/average_asym_scorer.hpp)

set(SOURCE_FILES "${RNASTRUCTURE_SOURCES}" "${LIBRNARY_SOURCES}")

set(PROJECT_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/lib")

add_library(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_INCLUDE_DIRS}")