// Created by max on 5/30/16.
/*
 * This file contains a clean interface to useful RNAstructure functions.
 */

#ifndef RNARK_RNA_STRUCTURE_HPP
#define RNARK_RNA_STRUCTURE_HPP

#include <string>
#include <memory>

#include "RNAstructure/rna_library.h"
#include "RNAstructure/algorithm.h"
#include "primary_structure.hpp"
#include "secondary_structure.hpp"
#include "energy.hpp"

namespace librnary {

/**
 * Given a path to the "data_tables" folder, loads into a datatable object.
 * Returns the datatable object pointer.
 * A unique pointer is used to prevent shenanigans.
 */
std::unique_ptr<datatable> LoadDatatable(const std::string &path);

/**
 * Given an RNA defined by its primary sequence, returns a corresponding structure object.
 * A unique pointer is used to prevent shenanigans.
 */
std::unique_ptr<structure> LoadStructure(const PrimeStructure &primary_seq);

std::unique_ptr<structure> LoadStructure(const PrimeStructure &primary_seq, const Matching &match);

std::unique_ptr<structure> LoadStructure(const std::string &ct_file);

std::unique_ptr<structure> LoadSeqFile(const std::string &seq_file);

// TODO: Make StructureToMatchings for structures with many possible matchings.
librnary::Matching StructureToMatching(structure &struc);

librnary::PrimeStructure StructureToPrimary(structure &struc);

/**
 * Given a datatable and structure, evaluates and saves the energy of the frist valid secondary structure in struct.
 * Uses the RNAstructure's standard energy model for multi-loops.
 */
energy_t RunEFN2(datatable &dt, structure &struc);

/**
 * Given a datatable and structure, evaluates and saves the energy of the frist valid secondary structure in struct.
 * Uses the simple, affine multi-loop model.
 */
energy_t RunEFN2WithSimpleMulti(datatable &dt, structure &struc);

/**
 * Run RNAstructure's version of the Zuker-Stiegler-like dynamic programming algorithm.
 * TODO: Add an option to generate several suboptimal tracebacks.
 */
energy_t RunMFEFold(datatable &dt, structure &struc, int two_loop_max_nts = 999999);

}

#endif //RNARK_RNA_STRUCTURE_HPP
