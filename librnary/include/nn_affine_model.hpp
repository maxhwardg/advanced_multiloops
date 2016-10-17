//
// Created by max on 6/15/16.
//

#ifndef RNARK_NN_AFFINE_MODEL_HPP
#define RNARK_NN_AFFINE_MODEL_HPP

#include "nn_model.hpp"
#include "energy.hpp"

namespace librnary {
/**
 * Implements a typical nearest neighbour model with an affine function ofr multi-loop free energy.
 * Because RNAstructure is used as a back-end, the parameters from RNAstructure are used for multi-loops here.
 */
class NNAffineModel: public NNModel {

public:
	energy_t MLClosure(const librnary::SSTree &sstree, librnary::SSTreeNodeId node_id) const override;
	energy_t MLClosure(const librnary::Surface &surf) const override;
	energy_t MLInitCost() const;
	energy_t MLBranchCost() const;
	energy_t MLUnpairedCost() const;
	NNAffineModel(const std::string &data_path, const PrimeStructure &_rna)
		: NNModel(data_path, _rna) {}
	NNAffineModel(const std::string &data_path)
		: NNModel(data_path) {}
};
}

#endif //RNARK_NN_AFFINE_MODEL_HPP
